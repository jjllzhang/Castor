//! Reproducible per-phase benchmark of the two Castor PCS openings.

use std::error::Error;
use std::time::{Duration, Instant};

use castor::{
    BrakedownQaPcs, ChallengeField as F, EvalTimings, FriConfig, PcsParams, ProofSizeBreakdown,
    ProximityBoundMode, ProximityBudget, QuasarQaPcs, Transcript, VerifyTimings,
};
use clap::{Parser, ValueEnum};
use p3_binary_field::{BinaryField16, BinaryField32, BinaryField64};
use p3_binary_qa::QaField;
use p3_field::{Algebra, ExtensionField};
use rand::{Rng, SeedableRng, rngs::StdRng};

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Scheme {
    Brakedown,
    Quasar,
    Both,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Bound {
    Ca,
    Mca,
    Third,
}

impl From<Bound> for ProximityBoundMode {
    fn from(bound: Bound) -> Self {
        match bound {
            Bound::Ca => Self::UniqueDecodingCA,
            Bound::Mca => Self::Johnson15MCA,
            Bound::Third => Self::ClassicalThird,
        }
    }
}

#[derive(Parser)]
#[command(about = "Benchmark QA multilinear commitments and authenticated evaluations")]
struct Args {
    #[arg(long, value_enum, default_value = "both")]
    scheme: Scheme,
    #[arg(long, default_value_t = 20)]
    m: usize,
    /// Log row count for Quasar. Brakedown uses floor(m/2).
    #[arg(long, default_value_t = 5)]
    r: usize,
    #[arg(long, default_value_t = 2)]
    c: usize,
    #[arg(long, default_value_t = 32)]
    witness_bits: usize,
    #[arg(long, default_value_t = 128)]
    mult_bits: usize,
    #[arg(long, value_enum, default_value = "ca")]
    bound_mode: Bound,
    #[arg(long, default_value_t = 103)]
    statistical_bits: usize,
    #[arg(long)]
    queries: Option<usize>,
    #[arg(long, default_value_t = 0)]
    main_cap_log: usize,
    #[arg(long, default_value_t = 1)]
    fri_log_blowup: usize,
    #[arg(long, default_value_t = 2)]
    fri_fold_step: usize,
    #[arg(long, default_value_t = 320)]
    fri_queries: usize,
    #[arg(long, default_value_t = 0)]
    fri_cap_log: usize,
    #[arg(long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value_t = 1)]
    warmup: usize,
    #[arg(long, default_value_t = 3)]
    iters: usize,
    #[arg(long, default_value_t = 20261119)]
    seed: u64,
}

fn samples<W: QaField>(rng: &mut StdRng, count: usize) -> Vec<W> {
    let mut bytes = [0; 16];
    (0..count)
        .map(|_| {
            rng.fill_bytes(&mut bytes[..W::NUM_BYTES]);
            W::from_le_byte_iter(bytes.into_iter())
        })
        .collect()
}

fn geometry(args: &Args, log_rows: usize) -> Result<(PcsParams, ProximityBudget), Box<dyn Error>> {
    if log_rows >= args.m {
        return Err("row variables must be fewer than total variables".into());
    }
    let budget = ProximityBudget::derive(
        args.m - log_rows,
        args.c,
        args.mult_bits,
        args.bound_mode.into(),
        args.statistical_bits,
        None,
    )?;
    let params = PcsParams::new(
        args.m,
        log_rows,
        args.c.ilog2() as usize,
        args.queries.unwrap_or(budget.queries),
        args.main_cap_log,
    )?;
    Ok((params, budget))
}

fn print_config(
    args: &Args,
    scheme: &str,
    params: PcsParams,
    budget: ProximityBudget,
    storage_bits: usize,
) {
    println!(
        "\n{scheme} | m={} N={} h={} K={} c={} L={} | threads={} seed={}",
        params.ell(),
        params.num_values(),
        params.rows(),
        params.cols(),
        params.inv_rate(),
        params.encoded_cols(),
        p3_maybe_rayon::prelude::current_num_threads(),
        args.seed
    );
    println!(
        "  field bits  witness={} multiplier={} commitment={} challenge=128",
        args.witness_bits, args.mult_bits, storage_bits
    );
    println!(
        "  bounds  mode={:?} delta={:.2} e={} t={} projection_bits={:.3} distance_failure_bits={:.3}",
        args.bound_mode,
        budget.distance,
        budget.radius,
        params.num_queries(),
        budget.projection_bits,
        budget.distance_failure_bits
    );
    println!(
        "  main Merkle cap log={}  benchmark warmup={} iters={}",
        args.main_cap_log, args.warmup, args.iters
    );
    if args.queries.is_some() {
        println!(
            "  explicit query override (calculated t={})",
            budget.queries
        );
    }
}

fn ms(t: Duration) -> f64 {
    t.as_secs_f64() * 1000.0
}

fn print_size(size: ProofSizeBreakdown) {
    for (label, bytes) in [
        ("Row combinations", size.row_combinations),
        ("Row evaluations", size.row_evaluations),
        ("Spectrum commitment", size.spectrum_root),
        ("Opened main columns", size.opened_main_columns),
        ("Main Merkle multiproof", size.main_multiproof),
        ("Cubic sumcheck", size.cubic_sumcheck),
        ("Endpoint claims", size.endpoint_claims),
        ("FRI-Binius opening", size.fri_opening),
        ("Total", size.total_bytes()),
    ]
    .into_iter()
    .filter(|&(_, bytes)| bytes > 0)
    {
        println!(
            "  {label:<26} {bytes:>9} B  {:>9.3} KiB",
            bytes as f64 / 1024.0
        );
    }
}

fn report(
    iteration: usize,
    total: [Duration; 3],
    eval: EvalTimings,
    verify: VerifyTimings,
    size: ProofSizeBreakdown,
    input: [usize; 2],
) {
    let [commit, open, check] = total;
    let [values, bytes] = input;
    println!(
        "sample={iteration} commit_ms={:.3} open_ms={:.3} verify_ms={:.3} input_MElem_s={:.3} input_MB_s={:.3} proof_bytes={}",
        ms(commit),
        ms(open),
        ms(check),
        (values as f64 / 1e6) / commit.as_secs_f64(),
        (bytes as f64 / 1e6) / commit.as_secs_f64(),
        size.total_bytes()
    );
    println!(
        "  Eval ms  row_eval={:.3} row_combinations={:.3} spectrum_commit={:.3} selectors={:.3} cubic={:.3} main_open={:.3} FRI={:.3}",
        ms(eval.row_evaluations),
        ms(eval.row_combinations),
        ms(eval.spectrum_commit),
        ms(eval.selectors),
        ms(eval.sumcheck),
        ms(eval.main_opening),
        ms(eval.fri_opening)
    );
    println!(
        "  Verify ms  main_auth={:.3} reencode={:.3} selector_lookup={:.3} FRI={:.3}",
        ms(verify.main_authentication),
        ms(verify.reencode),
        ms(verify.selector_lookup),
        ms(verify.fri_verification)
    );
    if iteration == 1 {
        print_size(size);
    }
}

fn run<W, E, M>(args: &Args) -> Result<(), Box<dyn Error>>
where
    W: QaField,
    E: QaField + ExtensionField<W> + ExtensionField<M>,
    M: QaField + p3_binary_dft::EncodableLevel,
    F: Algebra<W> + Algebra<E> + Algebra<M>,
{
    let mut rng = StdRng::seed_from_u64(args.seed);
    let witness = samples::<W>(&mut rng, 1 << args.m);
    let point = samples::<F>(&mut rng, args.m);
    let transcript = || Transcript::new(b"castor/benchmark");
    let input = [witness.len(), witness.len() * W::NUM_BYTES];
    if args.scheme != Scheme::Quasar {
        let (params, budget) = geometry(args, args.m / 2)?;
        print_config(args, "Brakedown", params, budget, E::NUM_BYTES * 8);
        let start = Instant::now();
        let pcs = BrakedownQaPcs::<E, M>::setup(params, &mut StdRng::seed_from_u64(args.seed ^ 1))?;
        println!("  setup_ms={:.3}", ms(start.elapsed()));
        for i in 0..args.warmup + args.iters {
            let table = witness.clone();
            let start = Instant::now();
            let (root, data) = pcs.commit::<W>(table)?;
            let commit = start.elapsed();
            let start = Instant::now();
            let value = pcs.evaluate::<W>(&data, &point)?;
            let claim_time = start.elapsed();
            let start = Instant::now();
            let (proof, eval) = pcs.eval::<W>(&data, &point, value, &mut transcript())?;
            let open = start.elapsed();
            let start = Instant::now();
            let verify = pcs.verify(&root, &point, value, &proof, &mut transcript())?;
            let check = start.elapsed();
            if i >= args.warmup {
                println!("  claim_evaluation_ms={:.3}", ms(claim_time));
                report(
                    i + 1 - args.warmup,
                    [commit, open, check],
                    eval,
                    verify,
                    proof.size(),
                    input,
                );
            }
        }
    }
    if args.scheme != Scheme::Brakedown {
        let (params, budget) = geometry(args, args.r)?;
        let fri = FriConfig {
            log_blowup: args.fri_log_blowup,
            fold_step: args.fri_fold_step,
            queries: args.fri_queries,
            cap_log: args.fri_cap_log,
        };
        print_config(args, "Quasar", params, budget, E::NUM_BYTES * 8);
        let fri_budget = fri.budget(params.log_cols())?;
        println!(
            "  FRI  B={} fold_step={} paths={} cap_log={} backend_bits={:.3}",
            1 << fri.log_blowup,
            fri.fold_step,
            fri.initial_queries(params.log_cols()),
            fri.cap_log,
            fri_budget.total_bits
        );
        let start = Instant::now();
        let pcs =
            QuasarQaPcs::<E, M>::setup(params, fri, &mut StdRng::seed_from_u64(args.seed ^ 1))?;
        println!("  setup_ms={:.3}", ms(start.elapsed()));
        let verifier = pcs.verifier();
        for i in 0..args.warmup + args.iters {
            let table = witness.clone();
            let start = Instant::now();
            let (root, data) = pcs.commit::<W>(table)?;
            let commit = start.elapsed();
            let start = Instant::now();
            let (value, proof, eval) = pcs.eval::<W>(&data, &point, &mut transcript())?;
            let open = start.elapsed();
            let start = Instant::now();
            let verify = verifier.verify(&root, &point, value, &proof, &mut transcript())?;
            let check = start.elapsed();
            if i >= args.warmup {
                report(
                    i + 1 - args.warmup,
                    [commit, open, check],
                    eval,
                    verify,
                    proof.size(),
                    input,
                );
            }
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if !(2..=24).contains(&args.m)
        || !matches!(args.c, 2 | 4)
        || args.threads == 0
        || args.iters == 0
        || args.warmup.checked_add(args.iters).is_none()
    {
        return Err("expected m=2..24, c=2 or 4, and positive threads and iterations".into());
    }
    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;
    #[cfg(not(feature = "parallel"))]
    if args.threads != 1 {
        return Err("set --threads 1 with --no-default-features".into());
    }
    type B64 = BinaryField64;
    match (args.witness_bits, args.mult_bits) {
        (16, 128) => run::<BinaryField16, F, F>(&args),
        (32, 128) => run::<BinaryField32, F, F>(&args),
        (64, 128) => run::<B64, F, F>(&args),
        (128, 128) => run::<F, F, F>(&args),
        (16, 64) if args.c == 4 => run::<BinaryField16, B64, B64>(&args),
        (32, 64) if args.c == 4 => run::<BinaryField32, B64, B64>(&args),
        (64, 64) if args.c == 4 => run::<B64, B64, B64>(&args),
        (128, 64) if args.c == 4 => run::<F, F, B64>(&args),
        _ => {
            Err("expected witness bits 16/32/64/128 and multiplier bits 128, or 64 with c=4".into())
        }
    }
}
