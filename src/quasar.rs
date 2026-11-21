//! Quasar opening by switching to one committed spectrum and a cubic sumcheck.

use core::any::Any;
use core::marker::PhantomData;
use std::time::Instant;

use p3_binary_dft::EncodableLevel;
use p3_binary_qa::{QaField, forward_transform};
use p3_field::{Algebra, ExtensionField, PrimeCharacteristicRing};
use p3_maybe_rayon::prelude::*;
use rand::Rng;

use crate::commit::{Code, ColumnVerifier, sample_spectra};
use crate::fri_binius::{FriBackend, FriVerifier};
use crate::merkle::stream_fields;
use crate::poly::{dot, eq_table, evaluate_rows, fold_row};
use crate::sumcheck;
use crate::tensor_kernel::{inverse_kernel, selector_tables, selector_values};
use crate::{
    CastorError, ChallengeField as F, CommitmentData, CubicMessage, Digest, EvalTimings, FriConfig,
    FriProof, MainOpening, PcsParams, ProofSizeBreakdown, Transcript, VerifyTimings,
};

pub struct QuasarQaPcs<E = F, M = E> {
    code: Code<E, M>,
    backend: FriBackend<M>,
    spectra: Vec<M>,
    promoted_spectra: Option<Vec<F>>,
}

/// Public Quasar verification key. It retains the code identifier and static FRI root.
#[derive(Clone, Copy)]
pub struct QuasarVerifier<E = F, M = E> {
    column: ColumnVerifier,
    fri: FriVerifier<M>,
    field: PhantomData<E>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QuasarProof<E = F, M = E, W = E> {
    pub row_evaluations: Vec<F>,
    pub spectrum_root: Digest,
    pub columns: MainOpening<W, E>,
    pub sumcheck: Vec<CubicMessage>,
    pub endpoints: Vec<F>,
    pub fri: FriProof<M>,
}

impl<E: QaField, M: QaField, W: QaField> QuasarProof<E, M, W> {
    pub fn size(&self) -> ProofSizeBreakdown {
        ProofSizeBreakdown {
            row_evaluations: 16 * self.row_evaluations.len(),
            spectrum_root: 32,
            opened_main_columns: self.columns.values_bytes(),
            main_multiproof: self.columns.authentication_bytes(),
            cubic_sumcheck: 48 * self.sumcheck.len(),
            endpoint_claims: 16 * self.endpoints.len(),
            fri_opening: self.fri.size_bytes(),
            ..ProofSizeBreakdown::default()
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.size().total_bytes());
        stream_fields(&self.row_evaluations, |bytes| out.extend_from_slice(bytes));
        out.extend_from_slice(&self.spectrum_root);
        self.columns.write_bytes(&mut out);
        for round in &self.sumcheck {
            stream_fields(&round.fields(), |bytes| out.extend_from_slice(bytes));
        }
        stream_fields(&self.endpoints, |bytes| out.extend_from_slice(bytes));
        self.fri.write_bytes(&mut out);
        out
    }
}

impl<E, M> QuasarQaPcs<E, M>
where
    E: QaField + ExtensionField<M>,
    M: QaField + EncodableLevel,
    F: Algebra<E> + Algebra<M>,
{
    pub fn setup(
        params: PcsParams,
        fri: FriConfig,
        rng: &mut impl Rng,
    ) -> Result<Self, CastorError> {
        Self::from_spectra(params, fri, sample_spectra(params, rng))
    }

    pub fn from_spectra(
        params: PcsParams,
        fri: FriConfig,
        spectra: Vec<M>,
    ) -> Result<Self, CastorError> {
        let code = Code::new(params, spectra.clone())?;
        let backend = FriBackend::setup(fri, params.log_cols(), &spectra)?;
        let promoted_spectra = if (&spectra as &dyn Any).is::<Vec<F>>() {
            None
        } else {
            Some(spectra.par_iter().copied().map(F::from).collect())
        };
        Ok(Self {
            code,
            backend,
            spectra,
            promoted_spectra,
        })
    }

    pub fn parameter_root(&self) -> Digest {
        self.backend.verifier.root
    }

    pub fn verifier(&self) -> QuasarVerifier<E, M> {
        QuasarVerifier {
            column: self.code.verifier,
            fri: self.backend.verifier,
            field: PhantomData,
        }
    }

    pub fn commit<W: QaField>(
        &self,
        witness: Vec<W>,
    ) -> Result<(Digest, CommitmentData<W, E>), CastorError>
    where
        E: ExtensionField<W>,
    {
        self.code.commit(witness)
    }

    pub fn commit_full(
        &self,
        witness: Vec<E>,
    ) -> Result<(Digest, CommitmentData<E, E>), CastorError> {
        self.code.commit_full(witness)
    }

    #[expect(
        clippy::type_complexity,
        reason = "The opening carries its claim and timings."
    )]
    pub fn eval<W: QaField>(
        &self,
        data: &CommitmentData<W, E>,
        point: &[F],
        transcript: &mut Transcript,
    ) -> Result<(F, QuasarProof<E, M, W>, EvalTimings), CastorError>
    where
        F: Algebra<W>,
    {
        self.code.check_input(data, point)?;
        let params = self.code.verifier.params;
        let (row_point, col_point) = point.split_at(params.log_rows());

        let start = Instant::now();
        let row_evaluations = evaluate_rows(&data.witness, params.rows(), &eq_table(col_point));
        let y = dot(&eq_table(row_point), &row_evaluations);
        let row_time = start.elapsed();
        let rho = self
            .verifier()
            .begin::<W>(&data.root(), point, y, &row_evaluations, transcript);

        let start = Instant::now();
        let mut spectrum = fold_row(&data.witness, &rho);
        let row_combinations = start.elapsed();

        let start = Instant::now();
        forward_transform(&mut spectrum, 1);
        let encoded = self.backend.commit(&spectrum);
        let spectrum_root = encoded.root();
        let spectrum_commit = start.elapsed();

        let start = Instant::now();
        let indices = self.verifier().get_columns(&spectrum_root, transcript);
        let columns = data.open_columns(&indices);
        let tau = column_challenge(&columns, transcript);
        let main_opening = start.elapsed();

        let start = Instant::now();
        let (selectors, power) =
            selector_tables(params.log_cols(), params.inv_rate(), &indices, tau);
        let selectors_time = start.elapsed();

        let start = Instant::now();
        let spectra = self.promoted_spectra.as_deref().unwrap_or_else(|| {
            (&self.spectra as &dyn Any)
                .downcast_ref::<Vec<F>>()
                .unwrap()
        });
        let cubic = sumcheck::prove(&spectrum, selectors, spectra, col_point, power, transcript);
        let sumcheck_time = start.elapsed();

        let start = Instant::now();
        let beta = backend_challenge(&cubic.endpoints, transcript);
        let fri = self.backend.prove(
            &encoded,
            spectrum,
            &self.spectra,
            &cubic.point,
            beta,
            transcript,
        );
        let fri_opening = start.elapsed();
        let proof = QuasarProof {
            row_evaluations,
            spectrum_root,
            columns,
            sumcheck: cubic.messages,
            endpoints: cubic.endpoints,
            fri,
        };
        Ok((
            y,
            proof,
            EvalTimings {
                row_evaluations: row_time,
                row_combinations,
                spectrum_commit,
                selectors: selectors_time,
                sumcheck: sumcheck_time,
                main_opening,
                fri_opening,
            },
        ))
    }

    pub fn verify<W: QaField>(
        &self,
        root: &Digest,
        point: &[F],
        y: F,
        proof: &QuasarProof<E, M, W>,
        transcript: &mut Transcript,
    ) -> Result<VerifyTimings, CastorError>
    where
        F: Algebra<W>,
    {
        self.verifier().verify(root, point, y, proof, transcript)
    }
}

impl<E: QaField, M: QaField + EncodableLevel> QuasarVerifier<E, M>
where
    F: Algebra<E> + Algebra<M>,
{
    pub fn verify<W: QaField>(
        &self,
        root: &Digest,
        point: &[F],
        y: F,
        proof: &QuasarProof<E, M, W>,
        transcript: &mut Transcript,
    ) -> Result<VerifyTimings, CastorError>
    where
        F: Algebra<W>,
    {
        use CastorError::InvalidProof;
        let params = self.column.params;
        if point.len() != params.ell()
            || proof.row_evaluations.len() != params.rows()
            || proof.sumcheck.len() != params.log_cols()
            || proof.endpoints.len() != params.inv_rate()
        {
            return Err(InvalidProof("Quasar shape"));
        }
        let (row_point, col_point) = point.split_at(params.log_rows());
        if dot::<F>(&eq_table(row_point), &proof.row_evaluations) != y {
            return Err(InvalidProof(
                "row evaluations do not give the claimed value",
            ));
        }

        let rho = self.begin::<W>(root, point, y, &proof.row_evaluations, transcript);
        let indices = self.get_columns(&proof.spectrum_root, transcript);
        let start = Instant::now();
        self.column.verify_columns(root, &indices, &proof.columns)?;
        let main_authentication = start.elapsed();
        let tau = column_challenge(&proof.columns, transcript);

        let mut initial = F::ZERO;
        let mut power = F::ONE;
        let rows = params.rows();
        let evaluations = proof
            .columns
            .systematic
            .chunks_exact(rows)
            .map(|col| dot(&rho, col))
            .chain(
                proof
                    .columns
                    .parity
                    .chunks_exact(rows)
                    .map(|col| dot(&rho, col)),
            );
        for evaluation in evaluations {
            initial += power * evaluation;
            power *= tau;
        }
        initial += power * dot::<F>(&rho, &proof.row_evaluations);
        let (final_sum, r) = sumcheck::reduce(&proof.sumcheck, initial, transcript);

        let start = Instant::now();
        let (selectors, power) = selector_values(&r, params.inv_rate(), &indices, tau);
        let selector_lookup = start.elapsed();
        let effective = selectors[0]
            + dot::<F>(&selectors[1..], &proof.endpoints[1..])
            + power * inverse_kernel(col_point, &r);
        if final_sum != proof.endpoints[0] * effective {
            return Err(InvalidProof("cubic sumcheck endpoint"));
        }
        let start = Instant::now();
        let beta = backend_challenge(&proof.endpoints, transcript);
        self.fri.verify(
            &proof.spectrum_root,
            &r,
            &proof.endpoints,
            beta,
            &proof.fri,
            transcript,
        )?;
        let fri_verification = start.elapsed();
        Ok(VerifyTimings {
            main_authentication,
            selector_lookup,
            fri_verification,
            ..VerifyTimings::default()
        })
    }

    fn begin<W: QaField>(
        &self,
        root: &Digest,
        point: &[F],
        y: F,
        rows: &[F],
        transcript: &mut Transcript,
    ) -> Vec<F> {
        self.column
            .begin::<W, E>(transcript, b"quasar/v1", root, point, y);
        self.fri.bind_parameters(transcript);
        transcript.absorb_fields(b"row-evaluations", rows);
        transcript.fields(b"rho", rows.len())
    }

    fn get_columns(&self, spectrum_root: &Digest, transcript: &mut Transcript) -> Vec<usize> {
        transcript.absorb(b"spectrum-root", spectrum_root);
        transcript.sorted_columns(
            b"main-columns",
            self.column.params.num_queries(),
            self.column.params.encoded_cols(),
        )
    }
}

fn column_challenge<W: QaField, E: QaField>(
    columns: &MainOpening<W, E>,
    transcript: &mut Transcript,
) -> F {
    transcript.absorb_main_columns(columns);
    transcript.field(b"selector-tau")
}

fn backend_challenge(endpoints: &[F], transcript: &mut Transcript) -> F {
    transcript.absorb_fields(b"cubic-endpoint-values", endpoints);
    transcript.field(b"fri-beta")
}

#[cfg(test)]
mod tests {
    use p3_binary_field::{BinaryField16, BinaryField32, BinaryField64, TowerLevel};
    use rand::{SeedableRng, rngs::StdRng};

    use super::*;
    use crate::poly::eq_index;

    fn assert_opening<E, M, W>(pcs: &QuasarQaPcs<E, M>, witness: Vec<W>)
    where
        E: QaField + ExtensionField<M> + ExtensionField<W>,
        M: QaField + EncodableLevel,
        W: QaField,
        F: Algebra<E> + Algebra<M> + Algebra<W>,
    {
        let params = pcs.code.verifier.params;
        let z: Vec<_> = (0..params.ell() as u128)
            .map(|i| F::from_repr((89 * i + 45).pow(3) << 57))
            .collect();
        let expected: F = witness
            .iter()
            .copied()
            .enumerate()
            .map(|(i, w)| eq_index(&z, i) * w)
            .sum();
        let (root, data) = pcs.commit(witness.clone()).unwrap();
        let transcript = || Transcript::new(b"quasar/integration");
        let (y, proof, _) = pcs.eval(&data, &z, &mut transcript()).unwrap();
        assert_eq!(y, expected);
        assert_eq!(
            proof.size().opened_main_columns,
            params.num_queries() * params.rows() * E::NUM_BYTES
                - proof.columns.systematic.len() * (E::NUM_BYTES - W::NUM_BYTES)
        );
        assert_eq!(proof.size().total_bytes(), proof.to_bytes().len());
        pcs.verify(&root, &z, y, &proof, &mut transcript()).unwrap();
        let (promoted_root, promoted) = pcs
            .commit_full(witness.iter().copied().map(E::from).collect())
            .unwrap();
        let (promoted_y, promoted_proof, _) = pcs.eval(&promoted, &z, &mut transcript()).unwrap();
        assert_eq!(promoted_y, y);
        if W::NUM_BYTES == E::NUM_BYTES {
            assert_eq!(promoted_root, root);
            assert_eq!(promoted_proof.to_bytes(), proof.to_bytes());
        } else {
            assert_ne!(promoted_root, root);
        }
        pcs.verify(&promoted_root, &z, y, &promoted_proof, &mut transcript())
            .unwrap();

        let rejected = |changed: &QuasarProof<E, M, W>| {
            pcs.verify(&root, &z, y, changed, &mut transcript())
                .is_err()
        };
        assert!(
            pcs.verify(&root, &z, y + F::ONE, &proof, &mut transcript())
                .is_err()
        );
        let mut changed = proof.clone();
        changed.row_evaluations[0] += F::ONE;
        assert!(rejected(&changed));
        let mut changed = proof.clone();
        changed.spectrum_root[0] ^= 1;
        assert!(rejected(&changed));
        let mut changed = proof.clone();
        changed.columns.parity[0] += E::ONE;
        assert!(rejected(&changed));
        let mut changed = proof.clone();
        changed.columns.cap[0][0] ^= 1;
        assert!(rejected(&changed));
        for (i, &round) in proof.sumcheck.iter().enumerate() {
            let mut changed = proof.clone();
            for tampered in [
                CubicMessage {
                    at_zero: round.at_zero + F::ONE,
                    ..round
                },
                CubicMessage {
                    at_omega: round.at_omega + F::ONE,
                    ..round
                },
                CubicMessage {
                    leading: round.leading + F::ONE,
                    ..round
                },
            ] {
                changed.sumcheck[i] = tampered;
                assert!(
                    rejected(&changed),
                    "accepted a changed cubic message at round {i}"
                );
            }
        }
        for i in 0..proof.endpoints.len() {
            let mut changed = proof.clone();
            changed.endpoints[i] += F::ONE;
            assert!(rejected(&changed));
        }
        let mut changed = proof.clone();
        changed.fri.witness_packets.values[0] += F::ONE;
        assert!(rejected(&changed));
        assert!(!proof.fri.layers.is_empty());
        for (layer, opening) in proof.fri.layers.iter().enumerate() {
            let mut changed = proof.clone();
            for (i, &value) in opening.packets.values.iter().enumerate() {
                changed.fri.layers[layer].packets.values[i] += F::ONE;
                assert!(
                    rejected(&changed),
                    "accepted a changed FRI symbol at layer {layer}, index {i}"
                );
                changed.fri.layers[layer].packets.values[i] = value;
            }
        }
        let mut changed = proof.clone();
        changed.fri.sumcheck[0][1] += F::ONE;
        assert!(rejected(&changed));
        let mut changed = proof;
        changed.fri.parameter_packets.values[0] += M::ONE;
        assert!(rejected(&changed));
    }

    #[test]
    fn full_multiplier_complete_and_mixed_witness_matches() {
        let params = PcsParams::new(10, 3, 1, 34, 2).unwrap();
        let mut rng = StdRng::seed_from_u64(239);
        let pcs = QuasarQaPcs::<F>::setup(
            params,
            FriConfig {
                queries: 19,
                ..FriConfig::default()
            },
            &mut rng,
        )
        .unwrap();
        assert_opening(
            &pcs,
            (0..1024)
                .map(|i| BinaryField32::from_repr(57 * i + 3))
                .collect(),
        );
    }

    #[test]
    fn narrow_multiplier_complete_and_uses_narrow_main_commitment() {
        let params = PcsParams::new(10, 3, 2, 29, 2).unwrap();
        let mut rng = StdRng::seed_from_u64(207);
        let pcs = QuasarQaPcs::<BinaryField64>::setup(
            params,
            FriConfig {
                queries: 19,
                ..FriConfig::default()
            },
            &mut rng,
        )
        .unwrap();
        assert_opening(
            &pcs,
            (0..1024)
                .map(|i| BinaryField16::from_repr(45 * i + 3))
                .collect(),
        );
    }
}
