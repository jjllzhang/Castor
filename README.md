# Castor

Castor implements two multilinear polynomial commitment schemes based on systematic quadratic quotient ring codes over binary tower fields. Both commit to a column-interleaved QA encoding with BLAKE3 Merkle trees and use `BinaryField128` for challenges.

| Scheme | Matrix layout | Evaluation proof |
|---|---|---|
| Brakedown | Near-square | Two row combinations and authenticated columns |
| Quasar | Flat | Row evaluations, one spectrum commitment, a cubic sumcheck, authenticated columns, and one batched FRI-Binius opening |

Both schemes spend $\Theta(N\log K)$ on the main commitment at fixed QA rate, where $K$ is the length of a message row. Quasar uses $h=\Theta(\log N)$ rows and works with a single folded spectrum. Its prover's opening work on the full tables is $O(N+K\log K)=O(N)$. With fixed rates, fold step, query budgets, and cap sizes, query authentication adds at most $O(\log^2 N)$ work, and the proof size and verifier work are also $O(\log^2 N)$. Brakedown verifies by re-encoding the two supplied rows.

## Usage

The input is a Boolean evaluation table. The first `r` variables occupy the low index bits. The `h = 2^r` entries in each column are contiguous. Choose the commitment field to match the wider of the witness and multiplier fields. Systematic columns use the original witness field for hashing and serialization. Parity columns use the commitment field. For example, a BF32 witness with BF128 multipliers uses four bytes for each opened systematic entry and sixteen bytes for each parity entry.

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

Both schemes use sorted, distinct column queries and shared Merkle paths. Optional caps are authenticated under the corresponding root. Proof byte totals are the lengths of the encoding returned by `to_bytes()`. Sorted queries determine how many systematic and parity columns appear in the proof. Mixed-width commitments bind the witness width in the transcript. The full-width format retains its original roots and proof bytes.

## Benchmark

The tables compare Castor e501780 with the optimization in 9661fc8. We ran native release builds on a shared AMD EPYC 9754 host with eight workers pinned to logical CPUs 0–7. Each process ran one warmup and three timed iterations, all verified. Entries give medians from the baseline and optimized builds. At m=26 the baseline executable has the same original implementation with its CLI size limit lifted.

Input generation and cloning are outside Commit. Brakedown Eval opens a known claim and the CLI logs claim evaluation separately. Quasar Eval includes row evaluation. Proof bytes omit the original commitment, public point, and claimed value. All cases use the CA radius and seed 20261119. The FRI backend uses inverse rate two, fold step two, cap log zero, and 320 paths except where stated.

### BF32 witnesses

| m | Scheme | QA inverse rate | Multiplier bits | Commit ms | Eval ms | Verify ms | Opened columns KiB | Proof KiB |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 16 | Brakedown | 2 | 128 | 1.889 → 1.489 | 2.877 → 2.331 | 13.474 → 12.113 | 1,708.00 → 1,069.00 | 1,718.44 → 1,079.50 |
| 16 | Brakedown | 4 | 64 | 1.933 → 2.021 | 2.446 → 1.864 | 7.942 → 8.200 | 624.00 → 545.00 | 644.34 → 565.12 |
| 16 | Brakedown | 4 | 128 | 2.761 → 2.824 | 2.163 → 1.756 | 5.658 → 5.347 | 716.00 → 584.00 | 734.72 → 603.19 |
| 16 | Quasar | 2 | 128 | 1.637 → 1.619 | 4.882 → 4.475 | 3.282 → 3.000 | 211.50 → 131.25 | 331.25 → 251.06 |
| 16 | Quasar | 4 | 64 | 1.919 → 2.314 | 4.691 → 5.090 | 2.470 → 2.464 | 77.75 → 68.12 | 210.25 → 201.62 |
| 16 | Quasar | 4 | 128 | 2.664 → 3.061 | 4.319 → 4.577 | 2.297 → 2.240 | 89.50 → 73.00 | 241.72 → 225.06 |
| 20 | Brakedown | 2 | 128 | 13.660 → 16.333 | 16.584 → 21.876 | 52.417 → 46.991 | 6,784.00 → 4,180.00 | 6,838.72 → 4,234.72 |
| 20 | Brakedown | 4 | 64 | 15.188 → 15.254 | 14.578 → 13.742 | 30.673 → 31.425 | 2,496.00 → 2,172.00 | 2,557.03 → 2,233.38 |
| 20 | Brakedown | 4 | 128 | 29.835 → 30.808 | 14.108 → 13.661 | 22.053 → 20.956 | 2,864.00 → 2,264.00 | 2,917.25 → 2,317.31 |
| 20 | Quasar | 2 | 128 | 16.147 → 15.551 | 23.565 → 23.274 | 4.923 → 4.592 | 211.50 → 131.62 | 531.69 → 446.56 |
| 20 | Quasar | 4 | 64 | 19.019 → 20.364 | 26.245 → 25.901 | 4.161 → 4.005 | 77.75 → 66.88 | 406.31 → 391.62 |
| 20 | Quasar | 4 | 128 | 35.053 → 36.103 | 29.170 → 25.852 | 3.729 → 3.633 | 89.50 → 70.38 | 415.34 → 391.00 |
| 24 | Brakedown | 2 | 128 | 226.927 → 213.972 | 220.793 → 209.548 | 207.204 → 186.519 | 27,072.00 → 17,856.00 | 27,247.56 → 18,030.50 |
| 24 | Brakedown | 4 | 64 | 235.482 → 256.046 | 214.482 → 212.581 | 121.641 → 124.236 | 9,952.00 → 8,624.00 | 10,128.53 → 8,799.56 |
| 24 | Brakedown | 4 | 128 | 488.542 → 504.025 | 210.081 → 207.681 | 88.108 → 82.750 | 11,456.00 → 9,248.00 | 11,616.50 → 9,407.28 |
| 24 | Quasar | 2 | 128 | 270.712 → 258.408 | 324.763 → 330.880 | 7.233 → 6.866 | 211.50 → 132.00 | 818.84 → 740.66 |
| 24 | Quasar | 4 | 64 | 328.847 → 351.063 | 390.169 → 377.278 | 6.136 → 6.086 | 77.75 → 69.38 | 671.78 → 667.00 |
| 24 | Quasar | 4 | 128 | 592.511 → 604.776 | 378.343 → 367.083 | 5.695 → 5.671 | 89.50 → 74.50 | 667.03 → 651.38 |

### Quasar witness width and scaling

| m | Witness bits | QA inverse rate | Commit ms | Eval ms | Verify ms | Opened columns KiB | Proof KiB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 16 | 2 | 15.595 → 14.406 | 25.108 → 25.099 | 5.160 → 4.605 | 211.50 → 111.31 | 529.59 → 431.00 |
| 20 | 32 | 2 | 16.147 → 15.551 | 23.565 → 23.274 | 4.923 → 4.592 | 211.50 → 131.62 | 531.69 → 446.56 |
| 20 | 64 | 2 | 19.140 → 16.942 | 19.090 → 18.710 | 4.850 → 4.597 | 211.50 → 159.50 | 530.50 → 479.38 |
| 20 | 128 | 2 | 20.334 → 17.994 | 21.166 → 20.964 | 4.877 → 4.900 | 211.50 → 211.50 | 531.94 → 531.94 |
| 24 | 128 | 2 | 273.637 → 275.018 | 302.708 → 298.269 | 7.361 → 6.935 | 211.50 → 211.50 | 814.72 → 814.72 |
| 24 | 128 | 4 | 559.287 → 559.654 | 357.827 → 324.067 | 5.787 → 5.649 | 89.50 → 89.50 | 670.50 → 670.50 |
| 26 | 32 | 2 | 1,099.177 → 1,064.133 | 1,435.435 → 1,359.608 | 8.440 → 8.694 | 211.50 → 132.38 | 992.09 → 918.91 |
| 26 | 32 | 4 | 2,344.124 → 2,408.694 | 1,709.780 → 1,588.139 | 7.526 → 6.980 | 89.50 → 70.75 | 833.72 → 812.69 |
| 26 | 128 | 2 | 1,128.807 → 1,121.503 | 1,274.403 → 1,195.909 | 8.368 → 8.417 | 211.50 → 211.50 | 993.19 → 993.19 |
| 26 | 128 | 4 | 2,270.605 → 2,275.378 | 1,522.367 → 1,425.936 | 6.932 → 6.798 | 89.50 → 89.50 | 824.84 → 824.84 |

At fixed query count, the opened-column saving depends on how many queried indices fall in the systematic block. The mixed-width encoder still produces a temporary full-width codeword. Commit compacts the parity blocks in parallel before building the main tree and retains the witness and parity blocks for opening. Full-width commitments retain the witness and the complete codeword, which Commit hashes directly. Peak RSS is measured separately and includes the temporary widened encoder output. Quasar reuses the original spectra at full width and promotes narrow spectra during Setup. The cubic sumcheck borrows the dynamic spectrum in its first round, then the FRI backend takes ownership of that buffer.

### Quasar FRI parameters

| m | Witness bits | Fold step | Cap log | Commit ms | Eval ms | Verify ms | Opened columns KiB | Proof KiB |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 32 | 2 | 0 | 16.147 → 15.551 | 23.565 → 23.274 | 4.923 → 4.592 | 211.50 → 131.62 | 531.69 → 446.56 |
| 20 | 32 | 3 | 0 | 16.101 → 15.386 | 22.710 → 22.414 | 4.688 → 4.384 | 211.50 → 132.75 | 518.16 → 440.88 |
| 20 | 32 | 4 | 0 | 16.239 → 15.849 | 22.247 → 22.002 | 4.972 → 4.626 | 211.50 → 130.88 | 567.62 → 484.69 |
| 20 | 32 | 4 | 8 | 16.254 → 16.115 | 22.392 → 22.043 | 5.006 → 4.706 | 211.50 → 138.38 | 587.50 → 515.38 |
| 24 | 128 | 2 | 0 | 273.637 → 275.018 | 302.708 → 298.269 | 7.361 → 6.935 | 211.50 → 211.50 | 814.72 → 814.72 |
| 24 | 128 | 3 | 0 | 276.248 → 276.772 | 277.408 → 277.866 | 6.581 → 6.495 | 211.50 → 211.50 | 762.66 → 762.66 |
| 24 | 128 | 4 | 0 | 276.011 → 276.780 | 273.371 → 270.080 | 6.859 → 6.868 | 211.50 → 211.50 | 806.84 → 806.84 |
| 24 | 128 | 4 | 8 | 276.274 → 275.985 | 274.744 → 272.459 | 6.835 → 6.769 | 211.50 → 211.50 | 825.53 → 825.53 |

Larger fold steps trade fewer layers for wider opened packets. The cap is serialized in full and paths are already shared. The FRI configuration is bound before sampling the main columns, so changing it also changes the outer queried indices at a fixed seed. The backend budget uses its own sampling bound. With 249 paths this term alone provides 80.16 bits, so the comparison uses 320 throughout.
