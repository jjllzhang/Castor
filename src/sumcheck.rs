//! Cubic spectral sumcheck with a separable inverse tensor kernel.

use std::borrow::Cow;

use p3_binary_qa::QaButterflyValue;
use p3_field::PrimeCharacteristicRing;
use p3_maybe_rayon::prelude::*;

use crate::tensor_kernel::inverse_factors;
use crate::{ChallengeField as F, Transcript};

/// Three transmitted evaluations. The running claim also determines the value at one.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CubicMessage {
    pub at_zero: F,
    pub at_omega: F,
    pub leading: F,
}

impl CubicMessage {
    pub(crate) fn fields(self) -> [F; 3] {
        [self.at_zero, self.at_omega, self.leading]
    }

    fn coefficients(self, sum: F) -> [F; 4] {
        let c0 = self.at_zero;
        let c3 = self.leading;
        let sigma = sum + c3;
        let c2 = sigma.mul_omega() + c0 + self.at_omega + c3;
        let c1 = sigma + c2;
        [c0, c1, c2, c3]
    }

    pub(crate) fn evaluate(self, sum: F, point: F) -> F {
        let [c0, c1, c2, c3] = self.coefficients(sum);
        c0 + point * (c1 + point * (c2 + point * c3))
    }
}

pub(crate) fn fold_table(table: &[F], point: F) -> Vec<F> {
    table
        .par_chunks_exact(2)
        .with_min_len(1024)
        .map(|pair| pair[0] + point * (pair[0] + pair[1]))
        .collect()
}

fn tensor(factors: &[[F; 2]]) -> Vec<F> {
    let mut table = Vec::with_capacity(1 << factors.len());
    table.push(F::ONE);
    for &[lo, hi] in factors {
        let width = table.len();
        table.resize(width * 2, F::ZERO);
        let (lower, upper) = table.split_at_mut(width);
        for (x, y) in lower.iter_mut().zip(upper) {
            *y = *x * hi;
            *x *= lo;
        }
    }
    table
}

fn add(a: [F; 3], b: [F; 3]) -> [F; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

struct SpectralTables<'a> {
    witness: Cow<'a, [F]>,
    selectors: Vec<Vec<F>>,
    spectra: Vec<Cow<'a, [F]>>,
}

impl SpectralTables<'_> {
    fn measure(&self, factors: &[[F; 2]], prefix: F) -> CubicMessage {
        let [lo, hi] = factors[0];
        let d0 = prefix * lo;
        let dw = prefix * (lo + (lo + hi).mul_omega());
        // Split the tensor tail into two small tables. Hoist the high-index
        // factor and multiply by the low-index factors within each tile.
        let (inner, outer) = factors[1..].split_at(factors.len() / 2);
        let inner = tensor(inner);
        let outer = tensor(outer);
        let sums = outer
            .par_iter()
            .enumerate()
            .with_min_len(4)
            .map(|(tile, &factor)| {
                let mut sums = [F::ZERO; 3];
                let mut target0 = F::ZERO;
                let mut targetw = F::ZERO;
                for (offset, &suffix) in inner.iter().enumerate() {
                    let i = 2 * (tile * inner.len() + offset);
                    let w0 = self.witness[i];
                    let delta_w = w0 + self.witness[i + 1];
                    let ww = w0 + delta_w.mul_omega();
                    let h0 = self.selectors[0][i];
                    let mut at0 = h0;
                    let mut atw = h0 + (h0 + self.selectors[0][i + 1]).mul_omega();
                    let mut leading = F::ZERO;
                    for (h, a) in self.selectors[1..].iter().zip(&self.spectra) {
                        let dh = h[i] + h[i + 1];
                        let da = a[i] + a[i + 1];
                        at0 += h[i] * a[i];
                        atw += (h[i] + dh.mul_omega()) * (a[i] + da.mul_omega());
                        leading += dh * da;
                    }
                    sums[0] += w0 * at0;
                    sums[1] += ww * atw;
                    sums[2] += delta_w * leading;
                    target0 += w0 * suffix;
                    targetw += ww * suffix;
                }
                sums[0] += target0 * (factor * d0);
                sums[1] += targetw * (factor * dw);
                sums
            })
            .par_fold_reduce(|| [F::ZERO; 3], add, add);
        CubicMessage {
            at_zero: sums[0],
            at_omega: sums[1],
            leading: sums[2],
        }
    }

    fn bind(&mut self, point: F) {
        self.witness = Cow::Owned(fold_table(&self.witness, point));
        for table in &mut self.selectors {
            *table = fold_table(table, point);
        }
        for table in &mut self.spectra {
            *table = Cow::Owned(fold_table(table, point));
        }
    }
}

pub(crate) struct CubicOpening {
    pub(crate) messages: Vec<CubicMessage>,
    pub(crate) point: Vec<F>,
    pub(crate) endpoints: Vec<F>,
}

pub(crate) fn prove(
    witness: &[F],
    selectors: Vec<Vec<F>>,
    spectra: &[F],
    target: &[F],
    power: F,
    transcript: &mut Transcript,
) -> CubicOpening {
    let mut tables = SpectralTables {
        witness: Cow::Borrowed(witness),
        selectors,
        spectra: spectra
            .chunks_exact(witness.len())
            .map(Cow::Borrowed)
            .collect(),
    };
    let factors = inverse_factors(target);
    let mut prefix = power;
    let mut messages = Vec::with_capacity(target.len());
    let mut point = Vec::with_capacity(target.len());
    for (axis, &[lo, hi]) in factors.iter().enumerate() {
        let message = tables.measure(&factors[axis..], prefix);
        let r = challenge(message, transcript);
        prefix *= lo + r * (lo + hi);
        tables.bind(r);
        point.push(r);
        messages.push(message);
    }
    let endpoints = std::iter::once(tables.witness[0])
        .chain(tables.spectra.iter().map(|a| a[0]))
        .collect();
    CubicOpening {
        messages,
        point,
        endpoints,
    }
}

pub(crate) fn reduce(
    messages: &[CubicMessage],
    mut claim: F,
    transcript: &mut Transcript,
) -> (F, Vec<F>) {
    let point = messages
        .iter()
        .map(|&message| {
            let r = challenge(message, transcript);
            claim = message.evaluate(claim, r);
            r
        })
        .collect();
    (claim, point)
}

fn challenge(message: CubicMessage, transcript: &mut Transcript) -> F {
    transcript.absorb_fields(b"cubic-message", &message.fields());
    transcript.field(b"cubic-challenge")
}

#[cfg(test)]
mod tests {
    use p3_binary_field::TowerLevel;
    use p3_field::Field;

    use super::*;
    use crate::tensor_kernel::{inverse_kernel, omega};

    #[test]
    fn all_cubic_rounds_match_full_table() {
        let bits = 5;
        let table = |seed: u128| {
            (0..1 << bits)
                .map(|i| F::from_repr((i * 307 + seed).pow(3)))
                .collect::<Vec<_>>()
        };
        let target: Vec<_> = (0..bits)
            .map(|i| F::from_repr(700 + 37 * i as u128))
            .collect();
        let power = F::from_repr(305);
        let mut tables = SpectralTables {
            witness: table(1).into(),
            selectors: vec![table(2), table(3), table(4)],
            spectra: vec![table(5).into(), table(6).into()],
        };
        let mut kernel: Vec<_> = (0..1 << bits)
            .map(|i| {
                let x: Vec<_> = (0..bits)
                    .map(|bit| F::from_bool(i & (1 << bit) != 0))
                    .collect();
                power * inverse_kernel(&target, &x)
            })
            .collect();
        let factors = inverse_factors(&target);
        let mut prefix = power;
        for (axis, &[lo, hi]) in factors.iter().enumerate() {
            let scalar = |point: F| {
                let interpolate = |table: &[F], i| table[i] + point * (table[i] + table[i + 1]);
                (0..tables.witness.len())
                    .step_by(2)
                    .map(|i| {
                        let mut weight =
                            interpolate(&kernel, i) + interpolate(&tables.selectors[0], i);
                        for (h, a) in tables.selectors[1..].iter().zip(&tables.spectra) {
                            weight += interpolate(h, i) * interpolate(a, i);
                        }
                        interpolate(&tables.witness, i) * weight
                    })
                    .sum::<F>()
            };
            let evaluations = [F::ZERO, F::ONE, omega(), omega() + F::ONE].map(scalar);
            // The sum over the four GF(4) points is the leading coefficient.
            let expected = CubicMessage {
                at_zero: evaluations[0],
                at_omega: evaluations[2],
                leading: evaluations.into_iter().sum(),
            };
            let actual = tables.measure(&factors[axis..], prefix);
            assert_eq!(actual, expected);
            let probe = F::from_repr(552 + 11 * axis as u128).try_inverse().unwrap();
            assert_eq!(
                actual.evaluate(evaluations[0] + evaluations[1], probe),
                scalar(probe)
            );
            tables.bind(probe);
            kernel = fold_table(&kernel, probe);
            prefix *= lo + probe * (lo + hi);
        }
    }
}
