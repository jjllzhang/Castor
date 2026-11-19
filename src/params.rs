//! Matrix geometry and explicit test parameters for a binary QA PCS.

use core::fmt;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CastorError {
    InvalidParameters(&'static str),
    InvalidInput(&'static str),
    InvalidProof(&'static str),
}

impl fmt::Display for CastorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (category, detail) = match self {
            Self::InvalidParameters(detail) => ("invalid parameters", detail),
            Self::InvalidInput(detail) => ("invalid input", detail),
            Self::InvalidProof(detail) => ("invalid proof", detail),
        };
        write!(f, "{category}: {detail}")
    }
}

impl std::error::Error for CastorError {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PcsParams {
    ell: usize,
    log_rows: usize,
    log_inv_rate: usize,
    num_queries: usize,
    cap_log: usize,
}

impl PcsParams {
    pub fn new(
        ell: usize,
        log_rows: usize,
        log_inv_rate: usize,
        num_queries: usize,
        cap_log: usize,
    ) -> Result<Self, CastorError> {
        use CastorError::InvalidParameters;

        if log_rows >= ell || !matches!(log_inv_rate, 1 | 2) {
            return Err(InvalidParameters(
                "a column bit and inverse rate 2 or 4 are required",
            ));
        }
        if ell
            .checked_add(log_inv_rate)
            .is_none_or(|bits| bits >= usize::BITS as usize - 4)
        {
            return Err(InvalidParameters("encoded matrix does not fit in a vector"));
        }
        if num_queries == 0 || num_queries > 1 << (ell - log_rows + log_inv_rate) {
            return Err(InvalidParameters(
                "query count must be in 1..=encoded columns",
            ));
        }
        if cap_log > ell - log_rows + log_inv_rate {
            return Err(InvalidParameters("Merkle cap is wider than the tree"));
        }
        Ok(Self {
            ell,
            log_rows,
            log_inv_rate,
            num_queries,
            cap_log,
        })
    }

    pub const fn ell(self) -> usize {
        self.ell
    }
    pub const fn log_rows(self) -> usize {
        self.log_rows
    }
    pub const fn log_cols(self) -> usize {
        self.ell - self.log_rows
    }
    pub const fn log_inv_rate(self) -> usize {
        self.log_inv_rate
    }
    pub const fn num_queries(self) -> usize {
        self.num_queries
    }
    pub const fn cap_log(self) -> usize {
        self.cap_log
    }
    pub const fn num_values(self) -> usize {
        1 << self.ell
    }
    pub const fn rows(self) -> usize {
        1 << self.log_rows
    }
    pub const fn cols(self) -> usize {
        1 << self.log_cols()
    }
    pub const fn inv_rate(self) -> usize {
        1 << self.log_inv_rate
    }
    pub const fn log_encoded_cols(self) -> usize {
        self.log_cols() + self.log_inv_rate
    }
    pub const fn encoded_cols(self) -> usize {
        1 << self.log_encoded_cols()
    }
}

/// Projection exception bound used to choose the column sampling radius.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProximityBoundMode {
    UniqueDecodingCA,
    Johnson15MCA,
    ClassicalThird,
}

/// Numerical planning bounds. The distance bound assumes independent uniform spectra.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProximityBudget {
    pub distance: f64,
    pub radius: usize,
    pub queries: usize,
    pub projection_bits: f64,
    pub distance_failure_bits: f64,
}

impl ProximityBudget {
    pub fn derive(
        log_cols: usize,
        inv_rate: usize,
        mult_bits: usize,
        mode: ProximityBoundMode,
        statistical_bits: usize,
        explicit_radius: Option<usize>,
    ) -> Result<Self, CastorError> {
        use CastorError::InvalidParameters;
        let distance: f64 = match (inv_rate, mult_bits) {
            (2, 128) => 0.40,
            (4, 128) => 0.70,
            (4, 64) => 0.50,
            _ => {
                return Err(InvalidParameters(
                    "no distance preset for this rate and multiplier",
                ));
            }
        };
        if log_cols == 0 || log_cols >= 48 || !(1..128).contains(&statistical_bits) {
            return Err(InvalidParameters("unsupported size or statistical budget"));
        }
        let k = 1 << log_cols;
        let bits = statistical_bits as f64;
        let distance_failure_bits = qa_distance_failure_bits(k, inv_rate, mult_bits, distance)?;
        // Equality after rounding does not establish that the bound meets the target.
        if distance_failure_bits <= bits {
            return Err(InvalidParameters(
                "QA distance bound misses the requested budget",
            ));
        }
        let length = k * inv_rate;
        let d0 = (distance * length as f64).ceil() as usize;
        let johnson = -((1.0 - distance).ln() / 3.0).exp_m1();
        let third = (d0 - 1) / 3;
        let projected = |e| projection_bits(length, distance, e, mode);
        let radius = if let Some(e) = explicit_radius {
            e
        } else {
            match mode {
                ProximityBoundMode::ClassicalThird => third,
                ProximityBoundMode::UniqueDecodingCA => {
                    ((johnson - 0.001) * length as f64) as usize
                }
                ProximityBoundMode::Johnson15MCA => {
                    let upper = ((johnson * length as f64) as usize).min((d0 - 1) / 2);
                    largest_mca_radius(length, distance, 1, upper, (128.0 - bits).exp2()).ok_or(
                        InvalidParameters("no MCA radius meets the projection budget"),
                    )?
                }
            }
        };
        if radius > (d0 - 1) / 2 || (mode == ProximityBoundMode::ClassicalThird && radius > third) {
            return Err(InvalidParameters("radius is outside the projection bound"));
        }
        let projection_bits = projected(radius);
        if projection_bits < bits {
            return Err(InvalidParameters("projection misses the requested budget"));
        }
        let sigma = radius as f64 / length as f64;
        let per_query = -(-sigma.min(distance - sigma)).ln_1p() / std::f64::consts::LN_2;
        // Opening every column has zero sampling error when the distance is positive.
        let queries = ((bits / per_query).ceil() as usize).min(length);
        Ok(Self {
            distance,
            radius,
            queries,
            projection_bits,
            distance_failure_bits,
        })
    }
}

fn largest_mca_radius(
    length: usize,
    distance: f64,
    lo: usize,
    hi: usize,
    threshold: f64,
) -> Option<usize> {
    if lo > hi {
        return None;
    }
    let [loss_low, _, slack_low] = mca_terms(length, distance, lo);
    let [loss_high, correction_high, slack_high] = mca_terms(length, distance, hi);
    // In the permitted radius range the loss grows and both the reciprocal
    // correction and slack decrease. Bound the projection over this interval.
    if loss_low + correction_high > threshold * slack_low {
        return None;
    }
    if loss_high + correction_high <= threshold * slack_high {
        return Some(hi);
    }
    let mid = lo + (hi - lo) / 2;
    if lo == hi {
        return None;
    }
    largest_mca_radius(length, distance, mid + 1, hi - 1, threshold)
        .or_else(|| largest_mca_radius(length, distance, lo, mid, threshold))
}

fn mca_terms(length: usize, distance: f64, radius: usize) -> [f64; 3] {
    let sigma = radius as f64 / length as f64;
    let a = 1.0 - sigma;
    let loss = -(1.5 * (-sigma).ln_1p()).exp_m1();
    let difference = a * sigma / (1.0 + a.sqrt());
    [
        length as f64 * loss,
        6.0_f64.max(2.0 / difference),
        a * a * a - (1.0 - distance),
    ]
}

fn projection_bits(length: usize, distance: f64, radius: usize, mode: ProximityBoundMode) -> f64 {
    if mode == ProximityBoundMode::ClassicalThird {
        return 128.0 - (radius as f64 + 1.0).log2();
    }
    let sigma = radius as f64 / length as f64;
    let a = 1.0 - sigma;
    let slack = a * a * a - (1.0 - distance);
    if radius == 0 || slack <= 0.0 {
        return 0.0;
    }
    let ratio = match mode {
        ProximityBoundMode::Johnson15MCA => {
            let [loss, correction, slack] = mca_terms(length, distance, radius);
            (loss + correction) / slack
        }
        ProximityBoundMode::UniqueDecodingCA => {
            let half = (distance * length as f64 / 2.0).floor();
            let d0 = (distance * length as f64).ceil();
            if radius as f64 >= half || half + radius as f64 >= d0 {
                return 0.0;
            }
            half.max(10.0 * (half + 1.0) / ((half + 1.0 - radius as f64) * slack))
        }
        ProximityBoundMode::ClassicalThird => unreachable!(),
    };
    (128.0 - ratio.log2()).max(0.0)
}

/// Negative base-two log of the finite-code distance failure bound, capped
/// below by zero. Assumes independent uniform spectra over the multiplier field.
/// Returns zero when the distance bound's beta > 0 hypothesis fails.
pub fn qa_distance_failure_bits(
    k: usize,
    rate: usize,
    field_bits: usize,
    delta: f64,
) -> Result<f64, CastorError> {
    if k == 0
        || rate < 2
        || field_bits == 0
        || !(0.0 < delta && delta < 1.0 - (-(field_bits as f64)).exp2())
    {
        return Err(CastorError::InvalidParameters(
            "invalid distance bound parameters",
        ));
    }
    let q_bits = field_bits as f64;
    let entropy_bits = -delta * delta.log2() - (1.0 - delta) * (1.0 - delta).log2();
    let ln_q_minus_one = q_bits * std::f64::consts::LN_2 + (-(-q_bits).exp2()).ln_1p();
    let beta_ln_q = (rate - 1) as f64 * q_bits * std::f64::consts::LN_2
        - rate as f64 * (delta * ln_q_minus_one + entropy_bits * std::f64::consts::LN_2);
    if beta_ln_q <= 0.0 {
        return Ok(0.0);
    }
    // Work in log space even when q^-beta underflows.
    let log_numerator = if beta_ln_q > 700.0 {
        (k as f64).ln() - beta_ln_q
    } else {
        (k as f64 * (-beta_ln_q).exp().ln_1p()).exp_m1().ln()
    };
    Ok(((ln_q_minus_one - log_numerator) / std::f64::consts::LN_2).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_budgets_at_m24() {
        use ProximityBoundMode::*;
        for (rate, third, mca, ca) in [
            (
                2,
                499,
                (160768, 160792, 429, 103.011),
                (163124, 423, 110.322),
            ),
            (
                4,
                269,
                (652288, 652670, 192, 103.015),
                (691152, 179, 108.515),
            ),
        ] {
            let derive = |mode, radius| ProximityBudget::derive(19, rate, 128, mode, 103, radius);
            assert_eq!(derive(ClassicalThird, None).unwrap().queries, third);
            let exact = derive(Johnson15MCA, Some(mca.0)).unwrap();
            assert_eq!(exact.queries, mca.2);
            assert!((exact.projection_bits - mca.3).abs() < 0.0005);
            let auto = derive(Johnson15MCA, None).unwrap();
            assert_eq!((auto.radius, auto.queries), (mca.1, mca.2));
            assert!(derive(Johnson15MCA, Some(auto.radius + 1)).is_err());
            let ca_auto = derive(UniqueDecodingCA, None).unwrap();
            assert_eq!((ca_auto.radius, ca_auto.queries), (ca.0, ca.1));
            assert!((ca_auto.projection_bits - ca.2).abs() < 0.0005);
        }
    }

    #[test]
    fn rejects_when_distance_bound_misses_statistical_budget() {
        let derive = |bits| {
            ProximityBudget::derive(11, 4, 64, ProximityBoundMode::ClassicalThird, bits, None)
        };
        assert_eq!(derive(103).unwrap().queries, 392);
        assert!(derive(112).is_ok());
        for bits in [113, 117] {
            assert_eq!(
                derive(bits).unwrap_err(),
                CastorError::InvalidParameters("QA distance bound misses the requested budget")
            );
        }
    }

    #[test]
    fn distance_and_mca_bounds_cover_supported_geometries() {
        for (rate, width, distance) in [(2, 128, 0.40), (4, 128, 0.70), (4, 64, 0.50)] {
            let bits = qa_distance_failure_bits(1 << 19, rate, width, distance).unwrap();
            assert!(bits > 103.0);
            for log in 2..=9 {
                let length = rate << log;
                let upper = ((distance * length as f64).ceil() as usize - 1) / 2;
                for bits in [103.0_f64, 118.0, 123.0, 127.0] {
                    let exhaustive = (1..=upper).rev().find(|&radius| {
                        projection_bits(length, distance, radius, ProximityBoundMode::Johnson15MCA)
                            >= bits
                    });
                    assert_eq!(
                        largest_mca_radius(length, distance, 1, upper, (128.0 - bits).exp2()),
                        exhaustive
                    );
                }
            }
        }
    }

    #[test]
    fn checks_matrix_and_query_shapes() {
        let p = PcsParams::new(5, 2, 2, 9, 1).unwrap();
        assert_eq!(
            (p.rows(), p.cols(), p.inv_rate(), p.encoded_cols()),
            (4, 8, 4, 32)
        );
        assert_eq!(p.num_values(), 32);
        assert!(PcsParams::new(5, 5, 1, 1, 1).is_err());
        assert!(PcsParams::new(5, 2, 2, 33, 1).is_err());
    }
}
