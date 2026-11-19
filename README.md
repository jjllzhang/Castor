# Castor

Castor implements two multilinear polynomial commitment schemes based on systematic quadratic quotient ring codes over binary tower fields. Both commit to a column-interleaved QA encoding with BLAKE3 Merkle trees and use `BinaryField128` for challenges.

| Scheme | Matrix layout | Evaluation proof |
|---|---|---|
| Brakedown | Near-square | Two row combinations and authenticated columns |
| Quasar | Flat | Row evaluations, one spectrum commitment, a cubic sumcheck, authenticated columns, and one batched FRI-Binius opening |

Both schemes spend $\Theta(N\log K)$ on the main commitment at fixed QA rate, where $K$ is the length of a message row. Quasar uses $h=\Theta(\log N)$ rows and works with a single folded spectrum. Its prover's opening work on the full tables is $O(N+K\log K)=O(N)$. With fixed rates, fold step, query budgets, and cap sizes, query authentication adds at most $O(\log^2 N)$ work, and the proof size and verifier work are also $O(\log^2 N)$. Brakedown verifies by re-encoding the two supplied rows.

## Usage

The input is a Boolean evaluation table. The first `r` variables occupy the low index bits. The `h = 2^r` entries in each column are contiguous. Choose the commitment field to match the wider of the witness and multiplier fields. A witness and multiplier at most 64 bits wide use BF64 main columns and pack each opened entry into eight bytes.

```bash
cargo test --locked --release
cargo run --locked --release --bin bench_pcs -- \
  --scheme both --m 20 --r 5 --c 2 \
  --witness-bits 32 --mult-bits 128 --bound-mode ca \
  --fri-log-blowup 1 --fri-fold-step 2 --threads 8
cargo bench --locked --bench pcs_comparison
```

`--m` sets the total number of variables. `--r` applies to Quasar. Brakedown uses `floor(m/2)` row variables. `--c` is the inverse QA code rate. The BF64 multiplier preset uses `--c 4 --mult-bits 64`. Serial builds use `--no-default-features` and `--threads 1` for the CLI. The repository enables native CPU features by default. Cargo uses the revision in `Cargo.lock` for the QA encoder and the same pinned Plonky3 revision throughout the dependency graph.

## Protocols and parameters

- `BrakedownQaPcs` opens a known claim with `eval`. It binds the statement, draws the random row weights, and accumulates the two row combinations in one pass through the witness. `evaluate` computes a claim from the committed witness when needed.
- `QuasarQaPcs` binds the row evaluations before drawing the random row weights. It commits to the resulting spectrum before sampling columns and authenticates those columns before batching the claims. The spectral sumcheck uses a pruned transpose butterfly and a factored inverse kernel. At the endpoint, the verifier evaluates the selectors with four-bit lookup tables. `verifier()` produces a compact verification key.
- The backend encodes Boolean table values as novel-basis coefficients over one common Cantor domain. It packs the published spectra into one static tree in their multiplier field. A second, quadratic sumcheck draws fresh challenges for the grouped additive RS folds. Each sampled packet path authenticates the homomorphic combination through the terminal constant.
- `ProximityBudget` checks the requested per-term budget against both the finite-code distance failure bound and the row-projection bound. The distance estimate assumes independently uniform spectra and uses the QA bound for $\beta>0$. The planner rejects a distance estimate that rounds to the requested threshold. The CLI reports the main column and FRI path budgets separately. These numerical bounds pertain to the interactive model. Fiat-Shamir knowledge security requires a separate analysis. The benchmark seeds exist for repeatability.

Both schemes use sorted, distinct column queries and shared Merkle paths. Optional caps are authenticated under the corresponding root. Proof byte totals are the lengths of the fixed-geometry encoding returned by `to_bytes()`.

## Benchmark

We ran both schemes on a shared AMD EPYC 9754 host with eight threads pinned to logical CPUs 0 through 7, using Rust 1.95 and native CPU features. Each configuration uses BF32 witnesses, the CA radius, and shared Merkle authentication paths with a cap log of zero. Quasar uses 32 rows, FRI blowup two, fold step two, and 320 starting paths.

The CLI runs one verified warmup and three verified measurements for each row. The table gives median execution times in milliseconds and the one-time setup measurement. Input generation and cloning are outside Commit. Brakedown Eval opens a known claim, and the CLI reports the cost of computing that claim separately. Quasar Eval derives the claim from its row evaluations. Proof sizes omit the original commitment, public point, and claimed value.

| Variables | Inverse rate | Multiplier bits | Scheme | Main queries | Setup ms | Commit ms | Eval ms | Verify ms | Proof KiB |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 16 | 2 | 128 | Brakedown | 427 | 0.019 | 1.748 | 3.012 | 13.498 | 1718.44 |
| 16 | 2 | 128 | Quasar | 423 | 0.451 | 2.099 | 6.096 | 3.341 | 331.25 |
| 16 | 4 | 128 | Brakedown | 179 | 0.056 | 2.662 | 1.919 | 5.749 | 734.72 |
| 16 | 4 | 128 | Quasar | 179 | 0.854 | 3.767 | 7.189 | 2.294 | 241.72 |
| 16 | 4 | 64 | Brakedown | 312 | 0.021 | 1.723 | 2.268 | 7.925 | 644.34 |
| 16 | 4 | 64 | Quasar | 311 | 0.855 | 2.451 | 6.681 | 2.466 | 210.25 |
| 20 | 2 | 128 | Brakedown | 424 | 0.071 | 13.774 | 17.685 | 52.526 | 6838.72 |
| 20 | 2 | 128 | Quasar | 423 | 4.291 | 17.784 | 23.568 | 4.913 | 531.69 |
| 20 | 4 | 128 | Brakedown | 179 | 0.199 | 30.119 | 13.697 | 22.108 | 2917.25 |
| 20 | 4 | 128 | Quasar | 179 | 10.444 | 34.824 | 26.076 | 3.754 | 415.34 |
| 20 | 4 | 64 | Brakedown | 312 | 0.080 | 14.183 | 13.807 | 30.590 | 2557.03 |
| 20 | 4 | 64 | Quasar | 311 | 5.866 | 18.977 | 26.086 | 4.084 | 406.31 |

At 20 variables and rate 1/2, Quasar reduces the proof from 6,838.72 to 531.69 KiB and verification from 52.526 to 4.913 ms. Its Eval takes 23.568 ms. The FRI opening and the authenticated main columns account for most of its proof bytes. The BF64 multiplier preset at rate 1/4 also changes the target distance and therefore the main query count.
