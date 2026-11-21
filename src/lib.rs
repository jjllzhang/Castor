//! Multilinear polynomial commitments built from binary QA codes.

pub use p3_binary_field::BinaryField128 as ChallengeField;

mod brakedown;
mod commit;
mod fri_binius;
mod merkle;
mod metrics;
mod params;
mod poly;
mod quasar;
mod sumcheck;
mod tensor_kernel;
mod transcript;

pub use brakedown::{BrakedownProof, BrakedownQaPcs};
pub use commit::CommitmentData;
pub use fri_binius::{FriBudget, FriConfig, FriProof, FriRound};
pub use merkle::{Digest, MainOpening, MultiOpening};
pub use metrics::{EvalTimings, ProofSizeBreakdown, VerifyTimings};
pub use params::{
    CastorError, PcsParams, ProximityBoundMode, ProximityBudget, qa_distance_failure_bits,
};
pub use quasar::{QuasarProof, QuasarQaPcs, QuasarVerifier};
pub use sumcheck::CubicMessage;
pub use tensor_kernel::{forward_kernel, inverse_kernel};
pub use transcript::Transcript;
