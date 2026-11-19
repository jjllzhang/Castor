use std::time::Duration;

/// Component sizes in the fixed-geometry proof encoding.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ProofSizeBreakdown {
    pub row_combinations: usize,
    pub row_evaluations: usize,
    pub opened_main_columns: usize,
    pub main_multiproof: usize,
    pub spectrum_root: usize,
    pub cubic_sumcheck: usize,
    pub endpoint_claims: usize,
    pub fri_opening: usize,
}

impl ProofSizeBreakdown {
    pub fn total_bytes(self) -> usize {
        self.row_combinations
            + self.row_evaluations
            + self.opened_main_columns
            + self.main_multiproof
            + self.spectrum_root
            + self.cubic_sumcheck
            + self.endpoint_claims
            + self.fri_opening
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct EvalTimings {
    pub row_evaluations: Duration,
    pub row_combinations: Duration,
    pub spectrum_commit: Duration,
    pub selectors: Duration,
    pub sumcheck: Duration,
    pub main_opening: Duration,
    pub fri_opening: Duration,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct VerifyTimings {
    pub main_authentication: Duration,
    pub reencode: Duration,
    pub selector_lookup: Duration,
    pub fri_verification: Duration,
}
