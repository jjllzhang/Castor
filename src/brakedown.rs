//! Direct Brakedown-style opening with random and evaluation row combinations.

use std::time::Instant;

use p3_binary_qa::{QaEncoder, QaField};
use p3_field::{Algebra, ExtensionField};
use rand::Rng;

use crate::commit::{Code, sample_spectra};
use crate::merkle::stream_fields;
use crate::poly::{dot, eq_table, evaluate, fold_rows};
use crate::{
    CastorError, ChallengeField as F, CommitmentData, Digest, EvalTimings, MultiOpening, PcsParams,
    ProofSizeBreakdown, Transcript, VerifyTimings,
};

pub struct BrakedownQaPcs<E = F, M = E> {
    code: Code<E, M>,
    challenge_rows: QaEncoder<F, M>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BrakedownProof<E = F> {
    pub proximity_row: Vec<F>,
    pub evaluation_row: Vec<F>,
    pub columns: MultiOpening<E>,
}

impl<E: QaField> BrakedownProof<E> {
    pub fn size(&self) -> ProofSizeBreakdown {
        ProofSizeBreakdown {
            row_combinations: 16 * (self.proximity_row.len() + self.evaluation_row.len()),
            opened_main_columns: self.columns.values_bytes(),
            main_multiproof: self.columns.authentication_bytes(),
            ..ProofSizeBreakdown::default()
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.size().total_bytes());
        stream_fields(&self.proximity_row, |bytes| out.extend_from_slice(bytes));
        stream_fields(&self.evaluation_row, |bytes| out.extend_from_slice(bytes));
        self.columns.write_bytes(&mut out);
        out
    }
}

impl<E, M> BrakedownQaPcs<E, M>
where
    E: QaField + ExtensionField<M>,
    M: QaField,
    F: Algebra<E> + Algebra<M>,
{
    pub fn setup(params: PcsParams, rng: &mut impl Rng) -> Result<Self, CastorError> {
        Self::from_spectra(params, sample_spectra(params, rng))
    }

    pub fn from_spectra(params: PcsParams, spectra: Vec<M>) -> Result<Self, CastorError> {
        let code = Code::new(params, spectra.clone())?;
        let challenge_rows = QaEncoder::new(params.log_cols(), params.log_inv_rate(), spectra);
        Ok(Self {
            code,
            challenge_rows,
        })
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

    /// Compute a claimed value before opening it.
    pub fn evaluate<W: QaField>(
        &self,
        data: &CommitmentData<W, E>,
        point: &[F],
    ) -> Result<F, CastorError>
    where
        F: Algebra<W>,
    {
        self.code.check_input(data, point)?;
        let (row_point, col_point) = point.split_at(self.code.verifier.params.log_rows());
        Ok(evaluate(
            &data.witness,
            &eq_table(row_point),
            &eq_table(col_point),
        ))
    }

    /// Open a known claim, accumulating both row combinations in one witness pass.
    pub fn eval<W: QaField>(
        &self,
        data: &CommitmentData<W, E>,
        point: &[F],
        y: F,
        transcript: &mut Transcript,
    ) -> Result<(BrakedownProof<E>, EvalTimings), CastorError>
    where
        F: Algebra<W>,
    {
        self.code.check_input(data, point)?;
        self.code
            .verifier
            .begin::<E>(transcript, b"brakedown/v1", &data.root(), point, y);
        let start = Instant::now();
        let rho = transcript.fields(b"rho", self.code.verifier.params.rows());
        let (row_point, col_point) = point.split_at(self.code.verifier.params.log_rows());
        let chi = eq_table(row_point);
        let [v, p] = fold_rows(&data.witness, &rho, &chi);
        if dot(&eq_table(col_point), &p) != y {
            return Err(CastorError::InvalidInput("incorrect claimed value"));
        }
        let row_combinations = start.elapsed();
        let (proof, main_opening) = self.finish_open(data, transcript, v, p);
        Ok((
            proof,
            EvalTimings {
                row_combinations,
                main_opening,
                ..EvalTimings::default()
            },
        ))
    }

    fn finish_open<W>(
        &self,
        data: &CommitmentData<W, E>,
        transcript: &mut Transcript,
        v: Vec<F>,
        p: Vec<F>,
    ) -> (BrakedownProof<E>, std::time::Duration) {
        transcript.absorb_fields(b"proximity-row", &v);
        transcript.absorb_fields(b"evaluation-row", &p);
        let start = Instant::now();
        let indices = transcript.sorted_columns(
            b"main-columns",
            self.code.verifier.params.num_queries(),
            self.code.verifier.params.encoded_cols(),
        );
        let columns = data.open_columns(&indices);
        let proof = BrakedownProof {
            proximity_row: v,
            evaluation_row: p,
            columns,
        };
        (proof, start.elapsed())
    }

    pub fn verify(
        &self,
        root: &Digest,
        point: &[F],
        y: F,
        proof: &BrakedownProof<E>,
        transcript: &mut Transcript,
    ) -> Result<VerifyTimings, CastorError> {
        let params = self.code.verifier.params;
        if point.len() != params.ell()
            || proof.proximity_row.len() != params.cols()
            || proof.evaluation_row.len() != params.cols()
        {
            return Err(CastorError::InvalidProof("Brakedown shape"));
        }
        let (row_point, col_point) = point.split_at(params.log_rows());
        if dot(&eq_table(col_point), &proof.evaluation_row) != y {
            return Err(CastorError::InvalidProof("evaluation row value"));
        }

        self.code
            .verifier
            .begin::<E>(transcript, b"brakedown/v1", root, point, y);
        let rho = transcript.fields(b"rho", params.rows());
        transcript.absorb_fields(b"proximity-row", &proof.proximity_row);
        transcript.absorb_fields(b"evaluation-row", &proof.evaluation_row);
        let indices =
            transcript.sorted_columns(b"main-columns", params.num_queries(), params.encoded_cols());
        let start = Instant::now();
        self.code
            .verifier
            .verify_columns(root, &indices, &proof.columns)?;
        let main_authentication = start.elapsed();

        let start = Instant::now();
        let u = self.challenge_rows.encode(&proof.proximity_row);
        let w = self.challenge_rows.encode(&proof.evaluation_row);
        let reencode = start.elapsed();
        let chi = eq_table(row_point);
        for (&j, column) in indices
            .iter()
            .zip(proof.columns.values.chunks_exact(params.rows()))
        {
            if u[j] != dot(&rho, column) || w[j] != dot(&chi, column) {
                return Err(CastorError::InvalidProof("row encoding consistency"));
            }
        }
        Ok(VerifyTimings {
            main_authentication,
            reencode,
            ..VerifyTimings::default()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p3_binary_field::{BinaryField16, BinaryField64, TowerLevel};
    use p3_field::PrimeCharacteristicRing;

    #[test]
    fn full_and_mixed_witnesses_have_identical_proofs_at_same_storage_width() {
        let params = PcsParams::new(8, 4, 2, 30, 2).unwrap();
        let spectra: Vec<_> = (0..48)
            .map(|i| BinaryField64::from_repr((i * 157 + 11) << 32))
            .collect();
        let pcs =
            BrakedownQaPcs::<BinaryField64, BinaryField64>::from_spectra(params, spectra.clone())
                .unwrap();
        let wide = BrakedownQaPcs::<F, BinaryField64>::from_spectra(params, spectra).unwrap();
        let witness: Vec<_> = (0..256)
            .map(|i| BinaryField16::from_repr(i * 131 + 7))
            .collect();
        let (root, mixed) = pcs.commit(witness.clone()).unwrap();
        let (promoted_root, promoted) = pcs
            .commit_full(witness.iter().copied().map(BinaryField64::from).collect())
            .unwrap();
        let (wide_root, wide_data) = wide
            .commit_full(witness.iter().copied().map(F::from).collect())
            .unwrap();
        let z: Vec<_> = (0..8).map(|i| F::from_repr(97 * i + 18)).collect();
        let transcript = || Transcript::new(b"brakedown/test");
        let y = pcs.evaluate(&mixed, &z).unwrap();
        let (proof, _) = pcs.eval(&mixed, &z, y, &mut transcript()).unwrap();
        let promoted_y = pcs.evaluate(&promoted, &z).unwrap();
        let (promoted_proof, _) = pcs
            .eval(&promoted, &z, promoted_y, &mut transcript())
            .unwrap();
        assert_eq!(
            (root, y, &proof),
            (promoted_root, promoted_y, &promoted_proof)
        );
        assert_eq!(proof.to_bytes().len(), proof.size().total_bytes());
        assert_eq!(proof.size().opened_main_columns, 30 * 16 * 8);
        pcs.verify(&root, &z, y, &proof, &mut transcript()).unwrap();
        let wide_y = wide.evaluate(&wide_data, &z).unwrap();
        let (wide_proof, _) = wide
            .eval(&wide_data, &z, wide_y, &mut transcript())
            .unwrap();
        assert_eq!(wide_y, y);
        wide.verify(&wide_root, &z, y, &wide_proof, &mut transcript())
            .unwrap();

        assert!(
            pcs.verify(&root, &z, y + F::ONE, &proof, &mut transcript())
                .is_err()
        );
        let mut changed = proof.clone();
        changed.evaluation_row[0] += F::ONE;
        assert!(
            pcs.verify(&root, &z, y, &changed, &mut transcript())
                .is_err()
        );
        let mut changed = proof.clone();
        changed.proximity_row[0] += F::ONE;
        assert!(
            pcs.verify(&root, &z, y, &changed, &mut transcript())
                .is_err()
        );
        let mut changed = proof.clone();
        changed.columns.values[0] += BinaryField64::ONE;
        assert!(
            pcs.verify(&root, &z, y, &changed, &mut transcript())
                .is_err()
        );
        let mut changed = proof;
        changed.columns.cap[0][0] ^= 1;
        assert!(
            pcs.verify(&root, &z, y, &changed, &mut transcript())
                .is_err()
        );
    }
}
