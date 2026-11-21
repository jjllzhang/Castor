//! Shared column-interleaved QA commitment.

use core::any::Any;

use p3_binary_qa::{QaEncoder, QaField};
use p3_field::ExtensionField;
use p3_maybe_rayon::prelude::*;
use rand::Rng;

use crate::{
    CastorError, ChallengeField as F, Digest, PcsParams, Transcript,
    merkle::{MainOpening, MerkleTree, update_fields, verify_split},
};

#[derive(Clone, Copy)]
pub(crate) struct ColumnVerifier {
    pub(crate) params: PcsParams,
    pub(crate) code_id: Digest,
}

pub(crate) struct Code<E, M> {
    pub(crate) verifier: ColumnVerifier,
    rows: QaEncoder<E, M>,
}

impl<E, M> Code<E, M>
where
    E: QaField + ExtensionField<M>,
    M: QaField,
{
    pub(crate) fn new(params: PcsParams, spectra: Vec<M>) -> Result<Self, CastorError> {
        if spectra.len() != params.encoded_cols() - params.cols() {
            return Err(CastorError::InvalidParameters(
                "spectra have the wrong length",
            ));
        }
        let mut hasher = blake3::Hasher::new();
        hasher.update(b"castor/spectra/v1");
        hasher.update(&(M::NUM_BYTES as u64).to_le_bytes());
        update_fields(&mut hasher, &spectra);
        let id = *hasher.finalize().as_bytes();
        let rows = QaEncoder::new(params.log_cols(), params.log_inv_rate(), spectra);
        Ok(Self {
            verifier: ColumnVerifier {
                params,
                code_id: id,
            },
            rows,
        })
    }

    pub(crate) fn commit<W: QaField>(
        &self,
        witness: Vec<W>,
    ) -> Result<(Digest, CommitmentData<W, E>), CastorError>
    where
        E: ExtensionField<W>,
    {
        self.check_witness(&witness)?;
        let log_rows = self.verifier.params.log_rows();
        let codeword = if let Some(full) = (&witness as &dyn Any).downcast_ref::<Vec<E>>() {
            self.rows.interleaved_encode_batch(full, log_rows)
        } else {
            self.rows.subfield_encode_batch(&witness, log_rows)
        };
        Ok(self.finish_commit(witness, codeword))
    }

    pub(crate) fn commit_full(
        &self,
        witness: Vec<E>,
    ) -> Result<(Digest, CommitmentData<E, E>), CastorError> {
        self.commit::<E>(witness)
    }

    fn check_witness<W: QaField>(&self, witness: &[W]) -> Result<(), CastorError> {
        if witness.len() != self.verifier.params.num_values() {
            return Err(CastorError::InvalidInput("witness has the wrong length"));
        }
        if E::NUM_BYTES != W::NUM_BYTES.max(M::NUM_BYTES) {
            return Err(CastorError::InvalidInput(
                "codeword field must match the wider of witness and multiplier",
            ));
        }
        Ok(())
    }

    fn finish_commit<W: QaField>(
        &self,
        witness: Vec<W>,
        mut codeword: Vec<E>,
    ) -> (Digest, CommitmentData<W, E>) {
        let rows = self.verifier.params.rows();
        let tree = if W::NUM_BYTES == E::NUM_BYTES {
            // The full-width encoder already has the exact Merkle leaf layout.
            MerkleTree::new(&codeword, rows)
        } else {
            // Each contiguous parity block is disjoint from the block being
            // overwritten. Finish shifting one block before reusing its space.
            let mut blocks = codeword.chunks_exact_mut(witness.len());
            let mut dest = blocks.next().unwrap();
            for source in blocks {
                dest.par_chunks_mut(4096)
                    .zip(source.par_chunks(4096))
                    .with_min_len(16)
                    .for_each(|(dst, src)| dst.copy_from_slice(src));
                dest = source;
            }
            codeword.truncate(codeword.len() - witness.len());
            codeword.shrink_to_fit();
            MerkleTree::new_split(&witness, &codeword, rows)
        };
        let root = tree.root();
        (
            root,
            CommitmentData {
                witness,
                encoded: codeword,
                tree,
                params: self.verifier.params,
                code_id: self.verifier.code_id,
            },
        )
    }

    pub(crate) fn check_input<W>(
        &self,
        data: &CommitmentData<W, E>,
        point: &[F],
    ) -> Result<(), CastorError> {
        if data.params != self.verifier.params
            || data.code_id != self.verifier.code_id
            || point.len() != self.verifier.params.ell()
        {
            return Err(CastorError::InvalidInput(
                "commitment data or point does not match setup",
            ));
        }
        Ok(())
    }
}

impl ColumnVerifier {
    pub(crate) fn begin<W: QaField, E: QaField>(
        &self,
        transcript: &mut Transcript,
        scheme: &[u8],
        root: &Digest,
        point: &[F],
        y: F,
    ) {
        transcript.absorb(b"scheme", scheme);
        let p = self.params;
        let shape = [
            p.ell(),
            p.log_rows(),
            p.log_inv_rate(),
            p.num_queries(),
            p.cap_log(),
            E::NUM_BYTES,
        ];
        transcript.absorb(
            b"shape",
            &shape.map(|v| v as u64).map(u64::to_le_bytes).concat(),
        );
        if W::NUM_BYTES != E::NUM_BYTES {
            transcript.absorb(b"systematic-width", &(W::NUM_BYTES as u64).to_le_bytes());
        }
        transcript.absorb(b"spectra", &self.code_id);
        transcript.absorb(b"main-root", root);
        transcript.absorb_fields(b"point", point);
        transcript.absorb_fields(b"value", &[y]);
    }

    pub(crate) fn verify_columns<W: QaField, E: QaField>(
        &self,
        root: &Digest,
        indices: &[usize],
        opening: &MainOpening<W, E>,
    ) -> Result<(), CastorError> {
        if verify_split(
            root,
            indices,
            opening,
            self.params.rows(),
            self.params.cols(),
            self.params.log_encoded_cols(),
            self.params.cap_log(),
        ) {
            Ok(())
        } else {
            Err(CastorError::InvalidProof("main column authentication"))
        }
    }
}

pub(crate) fn sample_spectra<M: QaField>(params: PcsParams, rng: &mut impl Rng) -> Vec<M> {
    let mut bytes = [0; 16];
    (params.cols()..params.encoded_cols())
        .map(|_| {
            rng.fill_bytes(&mut bytes[..M::NUM_BYTES]);
            M::from_le_byte_iter(bytes.into_iter())
        })
        .collect()
}

/// Retained witness, encoded columns, and Merkle data required to open a commitment.
pub struct CommitmentData<W, E> {
    pub(crate) witness: Vec<W>,
    // The mixed-width path removes the widened systematic prefix.
    encoded: Vec<E>,
    tree: MerkleTree,
    params: PcsParams,
    code_id: Digest,
}

impl<W: QaField, E: QaField> CommitmentData<W, E> {
    pub fn root(&self) -> Digest {
        self.tree.root()
    }

    pub(crate) fn open_columns(&self, indices: &[usize]) -> MainOpening<W, E> {
        let parity_start = if W::NUM_BYTES == E::NUM_BYTES {
            self.witness.len()
        } else {
            0
        };
        self.tree.open_split(
            &self.witness,
            &self.encoded[parity_start..],
            self.params.rows(),
            indices,
            self.params.cap_log(),
        )
    }
}
