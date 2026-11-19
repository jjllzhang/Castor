//! Multilinear polynomial commitments built from binary QA codes.

pub use p3_binary_field::BinaryField128 as ChallengeField;

mod brakedown;
mod commit;
mod merkle;
mod metrics;
mod params;
mod poly;
mod transcript;

pub use brakedown::{BrakedownProof, BrakedownQaPcs};
pub use commit::CommitmentData;
pub use merkle::{Digest, MultiOpening};
pub use metrics::{EvalTimings, ProofSizeBreakdown, VerifyTimings};
pub use params::{
    CastorError, PcsParams, ProximityBoundMode, ProximityBudget, qa_distance_failure_bits,
};
pub use transcript::Transcript;
