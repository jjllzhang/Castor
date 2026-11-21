//! Additive RS folding with packet Merkle commitments and a fresh sumcheck on the claim.

use core::marker::PhantomData;
use p3_binary_dft::EncodableLevel;
use p3_binary_pcs::{FoldAlphabet, fold_pair};

use p3_binary_field::TowerLevel;
use p3_commit::Encoder;
use p3_field::{Algebra, PrimeCharacteristicRing};
use p3_matrix::dense::RowMajorMatrix;
use p3_maybe_rayon::prelude::*;

use crate::merkle::{MerkleTree, stream_fields, verify};
use crate::poly::{eq_points, eq_table};
use crate::sumcheck::fold_table;
use crate::{CastorError, ChallengeField as F, Digest, MultiOpening, Transcript};

/// Grouped FRI folding parameters. Queries count complete paths from the initial packets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FriConfig {
    pub log_blowup: usize,
    pub fold_step: usize,
    pub queries: usize,
    pub cap_log: usize,
}

/// Error terms for the analyzed interactive backend, expressed in bits.
#[derive(Clone, Copy, Debug)]
pub struct FriBudget {
    pub path_sampling_bits: f64,
    pub folding_bits: f64,
    pub evaluation_bits: f64,
    pub total_bits: f64,
}

impl Default for FriConfig {
    fn default() -> Self {
        Self {
            log_blowup: 1,
            fold_step: 2,
            queries: 320,
            cap_log: 0,
        }
    }
}

impl FriConfig {
    pub(crate) fn check(self, log_cols: usize) -> Result<(), CastorError> {
        if log_cols == 0
            || !(1..=8).contains(&self.log_blowup)
            || !(1..=8).contains(&self.fold_step)
            || self.queries == 0
            || self.cap_log > 8
            || log_cols
                .checked_add(self.log_blowup)
                .is_none_or(|bits| bits >= usize::BITS as usize - 4)
        {
            return Err(CastorError::InvalidParameters(
                "invalid FRI geometry or query count",
            ));
        }
        Ok(())
    }

    pub fn budget(self, log_cols: usize) -> Result<FriBudget, CastorError> {
        self.check(log_cols)?;
        let initial_packets = 1 << (log_cols + self.log_blowup - self.fold_step.min(log_cols));
        let queries = self.initial_queries(log_cols);
        let sampling = if queries == initial_packets {
            0.0
        } else {
            0.8_f64.powf(queries as f64)
        };
        // With the fixed agreement threshold a=4/5 and relative RS distance
        // at least 1/2, the exact-set MCA slack is at least 3/250.
        let mut remaining = log_cols + self.log_blowup;
        let mut folding_numerator = 0.0;
        for arity in self.arities(log_cols) {
            remaining -= arity;
            let child_length = (1_usize << remaining) as f64;
            folding_numerator += arity as f64 * (0.285 * child_length + 24.0) / (3.0 / 250.0);
        }
        let inv_q = (-128.0_f64).exp2();
        let folding = folding_numerator * inv_q;
        let evaluation = 2.0 * log_cols as f64 * inv_q;
        Ok(FriBudget {
            path_sampling_bits: -sampling.log2(),
            folding_bits: -folding.log2(),
            evaluation_bits: -evaluation.log2(),
            total_bits: -(sampling + folding + evaluation).log2(),
        })
    }

    pub fn initial_queries(self, log_cols: usize) -> usize {
        self.queries
            .min(1 << (log_cols + self.log_blowup - self.fold_step.min(log_cols)))
    }

    fn arities(self, log_cols: usize) -> impl Iterator<Item = usize> {
        (0..log_cols)
            .step_by(self.fold_step)
            .map(move |used| self.fold_step.min(log_cols - used))
    }

    fn absorb(self, transcript: &mut Transcript) {
        let config = [self.log_blowup, self.fold_step, self.queries, self.cap_log];
        transcript.absorb(
            b"fri-geometry",
            &config.map(|v| v as u64).map(u64::to_le_bytes).concat(),
        );
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FriRound {
    pub root: Digest,
    pub packets: MultiOpening<F>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FriProof<M = F> {
    pub sumcheck: Vec<[F; 2]>,
    pub final_value: F,
    pub witness_packets: MultiOpening<F>,
    pub parameter_packets: MultiOpening<M>,
    pub layers: Vec<FriRound>,
}

impl<M: TowerLevel> FriProof<M> {
    pub fn size_bytes(&self) -> usize {
        32 * self.sumcheck.len()
            + 16
            + self.witness_packets.values_bytes()
            + self.witness_packets.authentication_bytes()
            + self.parameter_packets.values_bytes()
            + self.parameter_packets.authentication_bytes()
            + self
                .layers
                .iter()
                .map(|round| {
                    32 + round.packets.values_bytes() + round.packets.authentication_bytes()
                })
                .sum::<usize>()
    }

    pub(crate) fn write_bytes(&self, out: &mut Vec<u8>) {
        stream_fields(self.sumcheck.as_flattened(), |bytes| {
            out.extend_from_slice(bytes)
        });
        for layer in &self.layers {
            out.extend_from_slice(&layer.root);
        }
        stream_fields(&[self.final_value], |bytes| out.extend_from_slice(bytes));
        self.witness_packets.write_bytes(out);
        self.parameter_packets.write_bytes(out);
        for layer in &self.layers {
            layer.packets.write_bytes(out);
        }
    }
}

pub(crate) struct Encoded<T = F> {
    values: Vec<T>,
    width: usize,
    tree: MerkleTree,
}

impl<T: TowerLevel> Encoded<T> {
    fn new(values: Vec<T>, width: usize) -> Self {
        let tree = MerkleTree::new(&values, width);
        Self {
            values,
            width,
            tree,
        }
    }

    pub(crate) fn root(&self) -> Digest {
        self.tree.root()
    }

    fn open(&self, indices: &[usize], cap_log: usize) -> MultiOpening<T> {
        self.tree.open(
            &self.values,
            self.width,
            indices,
            cap_log.min(p3_util::log2_strict_usize(self.values.len() / self.width)),
        )
    }
}

pub(crate) struct FriBackend<M> {
    pub(crate) verifier: FriVerifier<M>,
    parameters: Encoded<M>,
}

#[derive(Clone, Copy)]
pub(crate) struct FriVerifier<M> {
    config: FriConfig,
    log_cols: usize,
    spectra_count: usize,
    pub(crate) root: Digest,
    field: PhantomData<M>,
}

impl<M: EncodableLevel> FriBackend<M>
where
    F: Algebra<M>,
{
    pub(crate) fn setup(
        config: FriConfig,
        log_cols: usize,
        spectra: &[M],
    ) -> Result<Self, CastorError> {
        config.check(log_cols)?;
        if log_cols + config.log_blowup > 1 << M::LOG_BITS {
            return Err(CastorError::InvalidParameters(
                "FRI domain exceeds the multiplier field",
            ));
        }
        let k = 1 << log_cols;
        let arity = config.fold_step.min(log_cols);
        let words: Vec<_> = spectra
            .chunks_exact(k)
            .map(|table| encode(table.to_vec(), config.log_blowup))
            .collect();
        let count = words.len();
        let mut packed = M::zero_vec((k << config.log_blowup) * count);
        packed
            .par_chunks_exact_mut(count)
            .with_min_len(1024)
            .enumerate()
            .for_each(|(i, row)| {
                for (entry, word) in row.iter_mut().zip(&words) {
                    *entry = word[i];
                }
            });
        let parameters = Encoded::new(packed, count << arity);
        Ok(Self {
            verifier: FriVerifier {
                config,
                log_cols,
                spectra_count: count,
                root: parameters.root(),
                field: PhantomData,
            },
            parameters,
        })
    }

    pub(crate) fn commit(&self, spectrum: &[F]) -> Encoded {
        let arity = self.verifier.config.fold_step.min(self.verifier.log_cols);
        Encoded::new(
            encode(spectrum.to_vec(), self.verifier.config.log_blowup),
            1 << arity,
        )
    }

    pub(crate) fn prove(
        &self,
        encoded: &Encoded,
        mut message: Vec<F>,
        spectra: &[M],
        point: &[F],
        beta: F,
        transcript: &mut Transcript,
    ) -> FriProof<M> {
        let powers = beta_powers(beta, self.verifier.spectra_count);
        let cols = message.len();
        for (&power, table) in powers.iter().zip(spectra.chunks_exact(cols)) {
            message
                .par_iter_mut()
                .zip(table.par_iter())
                .with_min_len(1024)
                .for_each(|(dst, &x)| *dst += power * x);
        }
        let mut eq = eq_table(point);
        let mut messages = Vec::with_capacity(self.verifier.log_cols);
        let mut layers: Vec<Encoded> = Vec::new();
        let arities: Vec<_> = self
            .verifier
            .config
            .arities(self.verifier.log_cols)
            .collect();
        for (layer, &arity) in arities.iter().enumerate() {
            let challenges: Vec<_> = (0..arity)
                .map(|_| {
                    let poly = quadratic_message(&message, &eq);
                    let r = quadratic_challenge(poly, transcript);
                    messages.push(poly);
                    message = fold_table(&message, r);
                    eq = fold_table(&eq, r);
                    r
                })
                .collect();
            if let Some(&next) = arities.get(layer + 1) {
                // The intermediate commitment uses the folded virtual codeword.
                // A final batch can go straight to the terminal constant.
                let folded = if let Some(source) = layers.last() {
                    F::fold_rounds(&source.values, &challenges)
                } else {
                    let virtual_word = encoded
                        .values
                        .par_iter()
                        .zip(
                            self.parameters
                                .values
                                .par_chunks_exact(self.verifier.spectra_count),
                        )
                        .with_min_len(1024)
                        .map(|(&value, parameters)| combine(value, parameters, &powers))
                        .collect::<Vec<_>>();
                    F::fold_rounds(&virtual_word, &challenges)
                };
                let next_layer = Encoded::new(folded, 1 << next);
                transcript.absorb(b"fri-layer-root", &next_layer.root());
                layers.push(next_layer);
            }
        }
        let final_value = message[0];
        transcript.absorb_fields(b"fri-final", &[final_value]);
        let indices = self.verifier.queries(transcript);
        let witness_packets = encoded.open(&indices, self.verifier.config.cap_log);
        let parameter_packets = self.parameters.open(&indices, self.verifier.config.cap_log);
        let layers = layers
            .into_iter()
            .zip(arities.into_iter().skip(1))
            .scan(indices, |indices, (layer, next_arity)| {
                *indices = shifted(indices, next_arity);
                Some(FriRound {
                    root: layer.root(),
                    packets: layer.open(indices, self.verifier.config.cap_log),
                })
            })
            .collect();
        FriProof {
            sumcheck: messages,
            final_value,
            witness_packets,
            parameter_packets,
            layers,
        }
    }
}

impl<M: EncodableLevel> FriVerifier<M>
where
    F: Algebra<M>,
{
    pub(crate) fn bind_parameters(&self, transcript: &mut Transcript) {
        self.config.absorb(transcript);
        transcript.absorb(b"fri-parameters", &self.root);
    }

    pub(crate) fn verify(
        &self,
        root: &Digest,
        point: &[F],
        claims: &[F],
        beta: F,
        proof: &FriProof<M>,
        transcript: &mut Transcript,
    ) -> Result<(), CastorError> {
        use CastorError::InvalidProof;
        let arities: Vec<_> = self.config.arities(self.log_cols).collect();
        if proof.sumcheck.len() != self.log_cols
            || proof.layers.len() + 1 != arities.len()
            || claims.len() != self.spectra_count + 1
        {
            return Err(InvalidProof("FRI shape"));
        }
        let powers = beta_powers(beta, self.spectra_count);
        let mut sum = combine::<F>(claims[0], &claims[1..], &powers);
        let mut challenges = Vec::with_capacity(self.log_cols);
        let mut messages = proof.sumcheck.iter();
        for (layer, &arity) in arities.iter().enumerate() {
            for &[p0, c2] in messages.by_ref().take(arity) {
                let r = quadratic_challenge([p0, c2], transcript);
                sum = p0 + r * (sum + c2 + r * c2);
                challenges.push(r);
            }
            if let Some(layer) = proof.layers.get(layer) {
                transcript.absorb(b"fri-layer-root", &layer.root);
            }
        }
        if sum != proof.final_value * eq_points(point, &challenges) {
            return Err(InvalidProof("FRI evaluation sumcheck endpoint"));
        }
        transcript.absorb_fields(b"fri-final", &[proof.final_value]);

        let indices = self.queries(transcript);
        let initial_log = self.log_cols + self.config.log_blowup - arities[0];
        let cap = self.config.cap_log.min(initial_log);
        let packet_size = 1 << arities[0];
        if !verify(
            root,
            &indices,
            &proof.witness_packets,
            packet_size,
            initial_log,
            cap,
        ) || !verify(
            &self.root,
            &indices,
            &proof.parameter_packets,
            packet_size * self.spectra_count,
            initial_log,
            cap,
        ) {
            return Err(InvalidProof("FRI initial packet authentication"));
        }
        let initial_values = proof
            .witness_packets
            .values
            .iter()
            .zip(
                proof
                    .parameter_packets
                    .values
                    .chunks_exact(self.spectra_count),
            )
            .map(|(&w, a)| combine(w, a, &powers))
            .collect::<Vec<_>>();
        let mut previous: Vec<_> = indices
            .iter()
            .copied()
            .zip(initial_values.chunks_exact(packet_size))
            .map(|(q, packet)| (q, fold_packet(q, packet, &challenges[..arities[0]])))
            .collect();
        let mut used = arities[0];
        let mut current_indices = indices;
        for (layer, &arity) in proof.layers.iter().zip(&arities[1..]) {
            let width = 1 << arity;
            current_indices = shifted(&current_indices, arity);
            let leaf_log = self.log_cols + self.config.log_blowup - used - arity;
            if !verify(
                &layer.root,
                &current_indices,
                &layer.packets,
                width,
                leaf_log,
                self.config.cap_log.min(leaf_log),
            ) {
                return Err(InvalidProof("FRI folded packet authentication"));
            }
            for (index, value) in previous {
                let packet = current_indices.binary_search(&(index >> arity)).unwrap();
                if layer.packets.values[packet * width + (index & (width - 1))] != value {
                    return Err(InvalidProof("FRI fold consistency"));
                }
            }
            previous = current_indices
                .iter()
                .copied()
                .zip(layer.packets.values.chunks_exact(width))
                .map(|(q, packet)| (q, fold_packet(q, packet, &challenges[used..used + arity])))
                .collect();
            used += arity;
        }
        if previous
            .iter()
            .any(|&(_, value)| value != proof.final_value)
        {
            return Err(InvalidProof("FRI terminal constant"));
        }
        Ok(())
    }

    fn queries(&self, transcript: &mut Transcript) -> Vec<usize> {
        let width = 1
            << (self.log_cols + self.config.log_blowup - self.config.fold_step.min(self.log_cols));
        transcript.sorted_columns(
            b"fri-queries",
            self.config.initial_queries(self.log_cols),
            width,
        )
    }
}

fn encode<T: EncodableLevel>(table: Vec<T>, log_blowup: usize) -> Vec<T> {
    // Each Boolean MLE value is one novel-basis coefficient. Adjacent-pair RS
    // folding then performs the corresponding MLE interpolation on coefficients.
    T::Encoder::default()
        .encode_batch(RowMajorMatrix::new(table, 1), log_blowup)
        .values
}

fn beta_powers(beta: F, count: usize) -> Vec<F> {
    (0..count)
        .scan(beta, |power, _| {
            let current = *power;
            *power *= beta;
            Some(current)
        })
        .collect()
}

fn combine<M: Copy>(value: F, parameters: &[M], powers: &[F]) -> F
where
    F: Algebra<M>,
{
    parameters
        .iter()
        .zip(powers)
        .fold(value, |sum, (&parameter, &power)| sum + power * parameter)
}

fn shifted(indices: &[usize], arity: usize) -> Vec<usize> {
    let mut next: Vec<_> = indices.iter().map(|&index| index >> arity).collect();
    next.dedup();
    next
}

fn fold_packet(index: usize, packet: &[F], challenges: &[F]) -> F {
    let mut values = packet.to_vec();
    for (round, &alpha) in challenges.iter().enumerate() {
        let width = values.len() / 2;
        for pair in 0..width {
            values[pair] = fold_pair(
                index * width + pair,
                alpha,
                values[2 * pair],
                values[2 * pair + 1],
            );
        }
        debug_assert_eq!(width << (round + 1), packet.len());
        values.truncate(width);
    }
    values[0]
}

fn add(a: [F; 2], b: [F; 2]) -> [F; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

fn quadratic_message(table: &[F], eq: &[F]) -> [F; 2] {
    table
        .par_chunks_exact(2)
        .zip(eq.par_chunks_exact(2))
        .with_min_len(1024)
        .map(|(f, e)| [f[0] * e[0], (f[0] + f[1]) * (e[0] + e[1])])
        .par_fold_reduce(|| [F::ZERO; 2], add, add)
}

fn quadratic_challenge(message: [F; 2], transcript: &mut Transcript) -> F {
    transcript.absorb_fields(b"fri-quadratic", &message);
    transcript.field(b"fri-fold")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dot;
    use p3_binary_field::TowerLevel;

    fn values(seed: u128, count: usize) -> Vec<F> {
        (0..count as u128)
            .map(|i| F::from_repr((197 * i + seed).pow(3) << 47))
            .collect()
    }

    #[test]
    fn multiplier_domain_embeds_in_challenge_domain() {
        use p3_binary_field::BinaryField64;
        let narrow: Vec<_> = (0..1_u64 << 11)
            .map(|i| BinaryField64::from_repr((99 * i + 28).pow(2) << 27))
            .collect();
        for blowup in 1..=3 {
            let wide = encode(narrow.iter().copied().map(F::from).collect(), blowup);
            let lifted: Vec<_> = encode(narrow.clone(), blowup)
                .into_iter()
                .map(F::from)
                .collect();
            assert_eq!(wide, lifted);
        }
        let budget = FriConfig::default().budget(19).unwrap();
        assert!(budget.path_sampling_bits >= 103.0);
        assert!((budget.folding_bits - 104.01181).abs() < 0.00001);
    }

    #[test]
    fn rs_folding_matches_mle_binding_and_disjoint_packets() {
        let table = values(37, 128);
        let point = values(331, 7);
        for blowup in 1..=3 {
            let mut code = encode(table.clone(), blowup);
            let mut table = table.clone();
            for &r in &point {
                for (q, pair) in code.chunks_exact(2).enumerate() {
                    assert_eq!(
                        fold_packet(q, pair, &[r]),
                        fold_pair(q, r, pair[0], pair[1])
                    );
                }
                code = F::fold_rounds(&code, &[r]);
                table = fold_table(&table, r);
                assert_eq!(code, encode(table.clone(), blowup));
            }
        }
    }

    #[test]
    fn batch_opens_dynamic_and_static_root_with_ragged_folds() {
        let bits = 7;
        let witness = values(21, 1 << bits);
        let spectra = values(51, 3 << bits);
        let point = values(89, bits);
        let beta = F::from_repr(373);
        let chi = eq_table(&point);
        let claims: Vec<_> = std::iter::once(witness.as_slice())
            .chain(spectra.chunks_exact(1 << bits))
            .map(|table| dot(&chi, table))
            .collect();
        for step in [1, 2, 3, 8] {
            let backend = FriBackend::setup(
                FriConfig {
                    fold_step: step,
                    queries: 13,
                    cap_log: 2,
                    ..FriConfig::default()
                },
                bits,
                &spectra,
            )
            .unwrap();
            let data = backend.commit(&witness);
            let transcript = || Transcript::new(b"backend-test");
            let proof = backend.prove(
                &data,
                witness.clone(),
                &spectra,
                &point,
                beta,
                &mut transcript(),
            );
            backend
                .verifier
                .verify(
                    &data.root(),
                    &point,
                    &claims,
                    beta,
                    &proof,
                    &mut transcript(),
                )
                .unwrap();
            let mut bytes = Vec::new();
            proof.write_bytes(&mut bytes);
            assert_eq!(proof.size_bytes(), bytes.len());
            let mut changed = proof.clone();
            changed.parameter_packets.values[0] += F::ONE;
            assert!(
                backend
                    .verifier
                    .verify(
                        &data.root(),
                        &point,
                        &claims,
                        beta,
                        &changed,
                        &mut transcript()
                    )
                    .is_err()
            );
            let mut changed = proof.clone();
            changed.final_value += F::ONE;
            assert!(
                backend
                    .verifier
                    .verify(
                        &data.root(),
                        &point,
                        &claims,
                        beta,
                        &changed,
                        &mut transcript()
                    )
                    .is_err()
            );
            let mut changed = proof.clone();
            if let Some(layer) = changed.layers.first_mut() {
                layer.packets.values[0] += F::ONE;
                assert!(
                    backend
                        .verifier
                        .verify(
                            &data.root(),
                            &point,
                            &claims,
                            beta,
                            &changed,
                            &mut transcript()
                        )
                        .is_err()
                );
            }
            for i in 0..claims.len() {
                let mut changed = claims.clone();
                changed[i] += F::ONE;
                assert!(
                    backend
                        .verifier
                        .verify(
                            &data.root(),
                            &point,
                            &changed,
                            beta,
                            &proof,
                            &mut transcript()
                        )
                        .is_err()
                );
            }
        }
    }
}
