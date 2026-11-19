//! Tensor kernels for the QA spectral transforms.

use p3_binary_field::TowerLevel;
use p3_binary_qa::QaButterflyValue;
use p3_field::PrimeCharacteristicRing;
use p3_maybe_rayon::prelude::*;

use crate::ChallengeField as F;

#[cfg(test)]
pub(crate) fn omega() -> F {
    F::from_repr(2)
}

#[inline]
fn omega2() -> F {
    F::from_repr(3)
}

/// Multilinear extension of the forward transform matrix.
pub fn forward_kernel(output: &[F], input: &[F]) -> F {
    assert_eq!(output.len(), input.len());
    output
        .iter()
        .zip(input)
        .map(|(&x, &y)| F::ONE + (omega2() + x) * y)
        .product()
}

/// Multilinear extension of the inverse transform matrix.
pub fn inverse_kernel(output: &[F], input: &[F]) -> F {
    assert_eq!(output.len(), input.len());
    output
        .iter()
        .zip(input)
        .map(|(&x, &y)| omega2() + x.mul_omega() + (F::ONE + x) * y)
        .product()
}

pub(crate) fn inverse_factors(output: &[F]) -> Vec<[F; 2]> {
    output
        .iter()
        .map(|&x| {
            let low = omega2() + x.mul_omega();
            [low, low + F::ONE + x]
        })
        .collect()
}

pub(crate) fn selector_tables(
    log_cols: usize,
    blocks: usize,
    indices: &[usize],
    tau: F,
) -> (Vec<Vec<F>>, F) {
    let width = 1 << log_cols;
    let mut entries = vec![Vec::new(); blocks];
    let mut power = F::ONE;
    for &index in indices {
        entries[index / width].push((index % width, power));
        power *= tau;
    }
    let sparse_layers = (indices.len() / blocks).max(1).ilog2() as usize;
    let tables = entries
        .into_iter()
        .map(|entries| transpose(log_cols, &entries, sparse_layers))
        .collect();
    (tables, power)
}

fn transpose(log_cols: usize, entries: &[(usize, F)], sparse_layers: usize) -> Vec<F> {
    let mut values = F::zero_vec(1 << log_cols);
    let mut active = Vec::with_capacity(entries.len());
    for &(index, weight) in entries {
        values[index] += weight;
        active.push(index);
    }
    let sparse_layers = sparse_layers.min(log_cols);
    for layer in 0..sparse_layers {
        let stride = 1 << layer;
        for index in &mut active {
            *index &= !stride;
        }
        active.sort_unstable();
        active.dedup();
        for &i in &active {
            let j = i | stride;
            let high = values[i].mul_omega() + values[j];
            values[i] += high;
            values[j] = high;
        }
        let low_count = active.len();
        for i in 0..low_count {
            active.push(active[i] | stride);
        }
    }
    drop(active);

    for layer in sparse_layers..log_cols {
        let stride = 1 << layer;
        // Within a layer, all tiles are disjoint even for strides exceeding one cache tile.
        let blocks = values.par_chunks_exact_mut(2 * stride);
        blocks
            .with_min_len((4096 / stride).max(1))
            .for_each(|block| {
                let (lo, hi) = block.split_at_mut(stride);
                lo.par_chunks_mut(2048)
                    .zip(hi.par_chunks_mut(2048))
                    .for_each(|(lo, hi)| {
                        transpose_dense(lo, hi);
                    });
            });
    }
    values
}

fn transpose_dense(lo: &mut [F], hi: &mut [F]) {
    #[cfg(all(target_arch = "x86_64", target_feature = "avx512f"))]
    type Simd = p3_binary_field::M512;
    #[cfg(all(
        target_arch = "x86_64",
        target_feature = "avx2",
        not(target_feature = "avx512f")
    ))]
    type Simd = p3_binary_field::M256;
    #[cfg(not(all(
        target_arch = "x86_64",
        any(target_feature = "avx2", target_feature = "avx512f")
    )))]
    type Simd = p3_binary_field::M128;

    let width = size_of::<Simd>() / size_of::<F>();
    let mut lo = lo.chunks_exact_mut(width);
    let mut hi = hi.chunks_exact_mut(width);
    for (left, right) in lo.by_ref().zip(hi.by_ref()) {
        let a = left.as_mut_ptr().cast::<Simd>();
        let b = right.as_mut_ptr().cast::<Simd>();
        // Each disjoint run has exactly one SIMD word's worth of transparent
        // u128 elements, for which all bit patterns represent field elements.
        unsafe {
            let low = a.read_unaligned();
            let high = low.mul_omega().bitxor(b.read_unaligned());
            a.write_unaligned(high.bitxor(low));
            b.write_unaligned(high);
        }
    }
    for (a, b) in lo.into_remainder().iter_mut().zip(hi.into_remainder()) {
        *b += a.mul_omega();
        *a += *b;
    }
}

pub(crate) fn selector_values(
    point: &[F],
    blocks: usize,
    indices: &[usize],
    tau: F,
) -> (Vec<F>, F) {
    let windows: Vec<_> = point
        .chunks(4)
        .map(|coordinates| {
            let mut table = [F::ONE; 16];
            for (bit, &x) in coordinates.iter().enumerate() {
                let high = 1 << bit;
                let factor = omega2() + x;
                table[high] = factor;
                for mask in 1..high {
                    table[mask | high] = table[mask] * factor;
                }
            }
            table
        })
        .collect();
    let width = 1 << point.len();
    let mut values = vec![F::ZERO; blocks];
    let mut power = F::ONE;
    for &index in indices {
        let absent = (width - 1) ^ (index % width);
        let mut factors = windows
            .iter()
            .enumerate()
            .map(|(i, window)| window[(absent >> (4 * i)) & 15]);
        let first = factors.next().unwrap_or(F::ONE);
        values[index / width] += power * factors.fold(first, |acc, factor| acc * factor);
        power *= tau;
    }
    (values, power)
}

#[cfg(test)]
mod tests {
    use p3_binary_qa::{forward_transform, inverse_transform};

    use super::*;
    use crate::poly::{dot, eq_table};

    fn binary_point(bits: usize, index: usize) -> Vec<F> {
        (0..bits)
            .map(|b| F::from_bool(index & (1 << b) != 0))
            .collect()
    }

    #[test]
    fn tensor_kernels_match_qa_transforms() {
        let point = [F::from_repr(15), F::from_repr(235), F::from_repr(7)];
        let input: Vec<_> = (0..8).map(|i| F::from_repr(i * i + 36)).collect();
        for (kernel, transform) in [
            (
                forward_kernel as fn(&[F], &[F]) -> F,
                forward_transform::<F> as fn(&mut [F], usize),
            ),
            (inverse_kernel, inverse_transform::<F>),
        ] {
            let mut output = input.clone();
            transform(&mut output, 1);
            let target: F = input
                .iter()
                .enumerate()
                .map(|(i, &v)| kernel(&point, &binary_point(3, i)) * v)
                .sum();
            assert_eq!(dot(&eq_table(&point), &output), target);
        }
    }

    #[test]
    fn pruned_transpose_and_windows_match_direct_kernels() {
        let indices: Vec<_> = (0..128).filter(|i| i % 3 == 0).collect();
        let point: Vec<_> = (0..5).map(|i| F::from_repr(100 + 117 * i)).collect();
        let tau = F::from_repr(96);
        let (tables, final_power) = selector_tables(5, 4, &indices, tau);
        let (evaluations, expected_power) = selector_values(&point, 4, &indices, tau);
        assert_eq!(final_power, expected_power);
        let mut direct = vec![vec![F::ZERO; 32]; 4];
        let mut direct_eval = [F::ZERO; 4];
        let mut power = F::ONE;
        for index in indices {
            let block = index / 32;
            let d = binary_point(5, index % 32);
            for (x, value) in direct[block].iter_mut().enumerate() {
                *value += power * inverse_kernel(&d, &binary_point(5, x));
            }
            direct_eval[block] += power * inverse_kernel(&d, &point);
            power *= tau;
        }
        assert_eq!(tables, direct);
        assert_eq!(evaluations, direct_eval);
    }
}
