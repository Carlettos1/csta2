# Physics and numerical test plan

Executable tests now cover the current workspace and every item in
[TODO.md](TODO.md). The matrices retain the original average and edge-case
requirements; the implementation map below identifies where they are tested.
“Average” means representative valid use, not merely an arithmetic mean.

## Test design and execution

- Use deterministic transition tests with scripted RNG draws for acceptance
  boundaries; compare against independently computed probabilities or state sums.
- Use tiny enumerable models as physics oracles. Fix Hamiltonians, boundary
  conditions, degeneracies, ensemble, units, and total/per-site normalization.
- Stochastic regression tests use fixed seeds, fixed finite budgets, and declared
  tolerances. Choose tolerances from variance and effective sample size, accounting
  for multiple comparisons; never rerun until a test passes. Keep expensive
  multi-seed convergence studies separate from fast deterministic checks.
- Compare floating-point values with absolute and relative tolerances appropriate
  to scale and precision. Check NaN/infinity explicitly. Energy offsets should
  preserve acceptance and centered fluctuations within representable ranges.
- Test errors for rollback and unchanged counters, not just `is_err()`. Test
  cancellation at proposal boundaries without relying on wall-clock sleeps.
- Add compile/run examples for simple public usage. Advanced capabilities must
  not require edits to an existing minimal `State` implementation.
- Run targeted crate tests first, then workspace tests for shared changes. Use
  the commands in [AGENTS.md](AGENTS.md), including all-features checks for Serde
  and explicit release/ignored runs for longer studies. Results are recorded below.

## Coverage locations

| Requirements | Executable coverage |
| --- | --- |
| B01, S01 | `csta_core/src/tests.rs`: all six vector types, conversions/operators, exceptional values, Serde. |
| B02, S07 | `csta_montecarlo/src/lib.rs`: initialization, tuple arities, empty iteration, uniform/Gaussian moments and isotropy. |
| B03, S07 | `csta_derive/src/lib.rs` diagnostics; `csta/tests/derive.rs` seeded runtime/compile-pass assertions; 9 compile-fail doctests in `csta_derive/README.md`. |
| B04–B05, S02–S03, C01 | `csta_metropolis/src/tests.rs`: scripted acceptance, negative beta, rollback, schedules, counters, symmetric/asymmetric adapters. |
| S04–S05 | `csta_metropolis/src/tests.rs`: batch/merged moments, covariance, IID/AR(1), blocks, nonlinear thermal errors, signed magnetization. |
| B06–B08, B13 | Existing `csta_wl/src/tests.rs` transition/grid/schedule/cancellation and seeded reference tests remain active. |
| B09–B10 | `csta_wl/src/wl_data.rs` and `analysis_tests.rs`: exact thermodynamics, continuous bin mass, linear/quadratic microcanonical derivatives. |
| B11 | `csta_wl/src/models.rs` reference enumeration, cached reversibility, oscillator support. |
| B12–B13 | `csta_wl/src/par_wl.rs`: window validation, exchange/merging, frozen mixing, worker failure and thread-count reproducibility. |
| B14, S08 | `csta/tests/wang_landau.rs`, root README doctests, all-target builds and runnable examples. |
| S06, C03, C07, N04 | `csta_metropolis/src/tests.rs`: independent square-lattice bond sums, cluster equilibrium/bond limits, local energies, interacting grand partition sums and extreme occupancy. |
| C02 | `csta_metropolis/src/tests.rs`: swap rejection, joint/marginal canonical distributions, seed replay, counts, single replica and invalid ladders. |
| C04–C05 | `csta_wl/src/analysis_tests.rs`: conditional/state sums, production observations, reweighting limits, exact uncertainty propagation, IDs/incomplete runs/support mismatches. |
| C06 | `csta/tests/checkpoints.rs`: canonical, serial WL, parallel WL (including masked support), temperature exchange, corruption, zero remaining budget and cancellation/resume; relative-path/atomic-save checks in core tests. |
| C07 | `csta_core/src/tests.rs`: minimum images versus brute force, translations, ties, lattice topology and invalid deserialized geometry. |
| N01 | `csta_metropolis/src/finite_size.rs`: known crossings/peaks, no/multiple crossings, trimmed ranges, bad grids, block-bootstrap reference; thermal/correlation tests cover dependent observations. |
| N02 | `csta_wl/src/analysis_tests.rs`: exact field sums/marginals, joint sampling, partial results, sparse support, offsets, numerical limits. |
| N03 | `csta_core/src/dynamics.rs`: analytic trajectory/order/drift, independent force differences, momentum, invalid controls, singular and unstable states. |

Long seeded WL convergence remains an explicit ignored release test. Snapshot
save is synchronous: tests save cancelled sessions, reject bad saves without
replacing a valid file, and verify exact resume; the API does not interrupt file
writes on a sampler cancellation flag. Histogram weight ESS flags concentration,
but neither tests nor API can infer unobserved support or independence from data
alone. Finite-size range tests check interpolation/peak ambiguity, not a hidden
critical-point fitting algorithm.

## Current implementation: retain or add

| Check | Average cases and oracle | Edge cases and required outcome |
| --- | --- | --- |
| [x] B01 — Vectors, all dimensions and precisions | Distinct components; constructors/accessors, tuple/array round trips, casts, owned/borrowed operators, dot/norm/distance against scalar calculations; Serde round trip | Zero, signed zero, subnormals, extreme finite components, NaN/infinity; `/` and `/=` agree under the documented contract. Casts preserve each coordinate. Distinct-coordinate tests guard the corrected cast defect. |
| [x] B02 — `MonteCarlo` and `Randomizable` | Fixed-seed reproducibility, uniform scalar/component support, tuple arities 2–8, finite `take(n)` length; component moments and cross-covariance | `take(0)` consumes no draws; custom deterministic sampler; draws near support endpoints; document iterator is unbounded and default thread RNG is not reproducible. |
| [x] B03 — Derive macro | Compile and sample named/tuple/unit structs, generic types, enums, weighted enums, ranges, lengths, defaults, arithmetic modifiers and `after`; compare generated sample values with a manual implementation using identical draws | Empty enum, union, malformed/unknown attributes, missing/mixed weights, negative/all-zero/nonfinite weights, empty/reversed ranges, zero lengths, invalid field dependencies. Invalid literal inputs produce compile-time diagnostics; dynamic invalid ranges are checked when sampling. |
| [x] B04 — Metropolis transition | Script downhill/uphill/equal-energy acceptance at positive beta; two-state stationary probabilities and detailed balance; accepted/rejected apply/revert counts | Beta zero, negative beta policy, extreme energy differences, invalid energies, RNG threshold neighbors; state and parameter-cache rollback. Checked execution restores state/cache/counters on invalid energies. |
| [x] B05 — Observers and rates | Compare single, 2/3/4, and dynamic observers using a deterministic trajectory; independent schedules; direct steps plus repeated runs | Zero budget, zero stride, burn-in equal to/exceeding budget, first/last measurement, no observers, no accepted moves, all accepted moves; rates use actual attempts and stay in [0,1] when defined (S03). |
| [x] B06 — Energy grids and support | Exact discrete lookup; continuous bin interiors and edges; valid support masks and initial DOS | Empty/duplicate/unsorted levels, nonfinite bounds, zero bins, adjacent representable boundary values, inclusive final edge, outside support, single level, all-masked support, impossible allocation; explicit errors without clamping. |
| [x] B07 — Serial WL transitions | Script accepted, rejected and self-loop moves; exactly one retained-bin visit per completed proposal; compare log acceptance thresholds | Invalid energy/DOS update or counter overflow rolls back state, caches and all data; off-support rejection; extreme finite log ratios. |
| [x] B08 — WL schedules and stopping | Verify flatness and visit minima, histogram reset, phase transition, SAMC proposal clock and exact production budget | Zero production, zero warm-up, equality at flatness threshold, zero histogram, exhausted budget, cancellation before/during sampling, bad controls, bounded initial-state search; partial results never imply completion. |
| [x] B09 — DOS thermodynamics | One/two-level analytic state sums for Z, E, E², variance, C, F and entropy; normalization/anchor and DOS-offset identities | Beta zero and finite negative beta; huge offsets, absent masses, all-absent data, invalid k_b/T, invalid anchor, nearly equal large energies, unrepresentable outputs; continuous bins use mass without another width factor. |
| [x] B10 — Microcanonical quantities | Linear/quadratic log-DOS on uniform/nonuniform energies versus analytic derivatives | Fewer than three levels, missing neighbors, endpoints, zero/negative slope; check documented errors, `None`, infinity and negative temperatures. |
| [x] B11 — Reference models | Enumerate tiny periodic 1D Ising systems; compare degeneracies and thermodynamic state sums; oscillator combinatorial DOS; cached energy equals recomputation after moves | Invalid spins, too few spins, no oscillators, zero cutoff, maximum occupancy, oscillator decrement at zero self-loops, reversible changes and matching finite truncation. |
| [x] B12 — Parallel WL | Serial/single-window equivalence, multiple walkers, known DOS offsets and slope-based joins, exchange ratio, cache movement, frozen-DOS mixing and count conservation | Missing coverage/overlap, unsupported exchange states, rejected swaps, unvisited overlap, incompatible grids, zero budgets, worker error/panic, cancellation, large exchange intervals; no merged partial success. |
| [x] B13 — Reproducibility and convergence | Repeat successful serial/parallel seeded runs; compare parallel output across thread counts; seeded Ising/oscillator estimates versus exact references | Slow-mixing/insufficient-budget runs stay explicitly partial; disconnected support is not “fixed” by silently masking bins; do not demand identical cancellation timing or cross-version RNG trajectories. |
| [x] B14 — Facade and examples | Compile current root/prelude exports, derive usage and examples; run bounded small canonical/WL demos; verify dependency and feature compatibility | Minimal `State` without `Randomizable` still runs through WL; default and Serde feature builds; debug-mode Ising runs and independent bond sums guard the former boundary defects. |

## Should-have acceptance tests

| TODO | Average cases and oracle | Edge cases and required outcome |
| --- | --- | --- |
| [x] S01 | B01 across Vec2/3/4 f32/f64; e.g. cast (1,2,3) without coordinate replication; normalization produces a unit vector in the same direction | Tiny/huge finite vectors normalize when mathematically representable; zero/nonfinite policy explicit; division forms agree, precision loss bounded by target precision. |
| [x] S02 | Enumerated two-state transition matrix obeys detailed balance for accepted beta domain; cache-mutating toy model restores snapshot on rejection | Negative-beta downhill counterexample; beta=0; NaN/infinite energies; exponent overflow/underflow; equal energies; numerical failure causes no partial mutation. If negative beta is unsupported, reject before moving. |
| [x] S03 | Exact expected attempted/accepted/measurement counts for repeated runs and mixed `step`/run use; cancellation returns actual work | Zero attempts gives documented unavailable/zero rate, not accidental NaN; stride=0 rejected; burn-in endpoints unambiguous; counter overflow checked; all observer forms agree. |
| [x] S04 | Online and merged summaries match independent batch calculations; IID Gaussian error scales as n^(-1/2); stationary AR(1) fixture with known correlation tests error inflation and ESS | Empty/single/constant series, large offset with small variance, incomplete blocks, nonfinite samples, short/strongly correlated chains; insufficient data reported. Test negative correlation under the chosen ESS convention. |
| [x] S05 | Enumerate small Ising state sums at chosen beta/h; verify energy and M moments, C=k_b beta² Var(E), total chi=beta Var(M) for field term -hM, and U4=1-<M⁴>/(3<M²>²) | Zero M² gives undefined Binder result; total versus per-site scaling; |M| is not substituted into signed-M susceptibility; beta=0 and low-temperature variance; invalid site count/units; errors include sampling correlation. |
| [x] S06 | Enumerate a tiny 2D lattice with an independent unique-bond list; all-aligned energy; nonzero J/h; open and periodic cases; local updates equal full energy differences | Corners/row seams, odd sizes, minimal supported periodic dimensions, invalid dimensions/spins, size-product overflow; no unsigned underflow or duplicate bonds outside the declared small-lattice convention. |
| [x] S07 | B03 compile-pass cases; known 1:3 enum weights; isotropic vectors have zero component means and <u_i u_j>=delta_ij/d; Gaussian mean/variance | Zero-weight variants never drawn; all-zero/negative/nonfinite weights and invalid syntax produce useful diagnostics; unit directions have finite norm 1; zero/invalid Gaussian scale follows documented policy; existing uniform-vector sampling unchanged. |
| [x] S08 | Run introductory examples from a fresh workspace setup; compare printed small-model values to exact references; resolve local doc links and package-qualified commands | Partial/cancelled examples do not print success; no missing sibling-repository dependency; no implicit absolute DOS normalization, missing units, or unseeded “reproducible” example. |

## Could-have acceptance tests

| TODO | Average cases and oracle | Edge cases and required outcome |
| --- | --- | --- |
| [x] C01 | Enumerated asymmetric three-state proposal matrix reaches prescribed Boltzmann probabilities; WL correction agrees with target inverse-DOS weights; symmetric adapter reproduces old trajectories | Zero reverse probability rejects; impossible forward move, invalid log ratio and constrained self-loop semantics checked; rejected moves restore caches; existing `State` compiles unchanged. |
| [x] C02 | Enumerated two-replica joint distribution and scripted swap acceptance; fixed-seed replay; temperature-specific observables agree with exact canonical values | One replica, duplicate/invalid ladder entries under documented policy, identical energies, extreme beta gaps, rejected swap, cancellation; states/caches move together and swap counts do not inflate local counts. |
| [x] C03 | Small-lattice equilibrium distribution and moments match enumeration; cluster construction on scripted bonds; measure processed spins alongside cluster moves | Beta=0, low-temperature full cluster, periodic seams, repeated neighbor discovery; negative J, nonzero field and invalid beta rejected as appropriate; no infinite growth or double insertion. |
| [x] C04 | Enumerated joint (E,M) data recovers conditional M moments and reweighted canonical values; source-beta identity for histogram reweighting | Missing conditional bins, no overlap, concentrated weights, empty data and extreme target beta yield explicit unavailable/low-reliability diagnostics; no fabricated M from energy-only DOS. |
| [x] C05 | Independent seeded runs aligned to one anchor produce observable summaries and uncertainty; synthetic known DOS perturbations test propagation | One run cannot estimate between-run variance; missing anchors, differing support/grids and incomplete runs handled explicitly; shared replicas not counted as independent; additive DOS offsets do not inflate uncertainty. |
| [x] C06 | Continuous n+m steps equal checkpoint-at-n then resume-m, including RNG, caches, DOS, histograms and clocks; parallel checkpoint at a defined safe barrier | Corrupt/truncated data, incompatible version/model/grid/RNG, invalid counters, cancellation during save, zero remaining budget; reject without corrupting a valid snapshot; DOS-only restart does not claim exact continuation. |
| [x] C07 | Wrapping invariant under integer box translations; minimum-image displacement against brute-force nearby images; local delta E matches full recomputation | Negative positions, many-box displacement, half-box ties, zero/negative/nonfinite lengths, tiny lattices, corner moves, cache drift over long move/revert sequences; tie and neighbor multiplicity conventions explicit. |

## Neat-feature acceptance tests

| TODO | Average cases and oracle | Edge cases and required outcome |
| --- | --- | --- |
| [x] N01 | Synthetic curves with known crossing/peak; independent bootstrap/block fixtures for propagated uncertainties; tiny-lattice observables match enumeration before scaling analysis | No crossing, multiple crossings, too few sizes, mismatched temperature grids, strongly correlated data, fit-range changes; report ambiguity rather than a spurious critical point. |
| [x] N02 | Enumerate g(E0,M), marginalize to energy DOS at h=0, and compare sum g exp[-beta(E0-hM)] at nonzero field with direct state sums | Sparse/inaccessible cells, missing support, degenerate energies, field extremes, allocation limits, normalization offsets; no double field term and no changes to 1D APIs. |
| [x] N03 | Harmonic oscillator trajectory versus analytic solution; velocity-Verlet global error decreases approximately fourfold when dt halves at fixed duration; bounded NVE energy error; pair forces versus independent energy finite differences | Zero steps, invalid dt/mass, coincident particles for singular potentials, periodic crossings, large unstable dt; momentum conservation for isolated pair forces; fail/report nonfinite state without calling unstable motion physical. |
| [x] N04 | Enumerate a few-site noninteracting lattice gas; occupation p=1/(1+exp(-beta mu)) for zero site energy; interacting tiny lattice versus direct grand partition sum | Empty/full occupancy, extreme chemical potential, beta=0, rejected insertion/deletion, invalid controls; proposal corrections include move-selection factors; particle counts, caches and rollback remain consistent. |

## Completion criteria

For each roadmap item, implement its matching row and relevant baseline rows;
record the exact commands/results in the change description. Preserve the
minimal usage example and explain any intentional contract change. Keep long
statistical studies reproducible and separate from fast checks. A passing
statistical regression demonstrates agreement for its tested model and budget,
not convergence for every user-supplied system.

## Verified results

Validated with Rust/Cargo 1.98.0 and the checked-in dependency lockfile:

- `cargo fmt --all -- --check`: passed.
- `cargo check --workspace --all-targets --all-features --offline`: passed.
- `cargo test --workspace --offline`: 87 tests/doctests passed, one long study ignored.
- `cargo test --workspace --all-features --offline`: 96 tests/doctests passed, one long study ignored.
- `cargo clippy --workspace --all-targets --all-features --offline -- -D warnings`: passed.
- `cargo test -p csta_wl --release --offline -- --ignored`: the long convergence study passed.
- `cargo test -p csta --all-features --offline --example state`: reversible-proposal regression passed.
- Ran `ising` in debug mode, the `csta` simulation/derive examples in release
  mode, and the checkpoint example with its feature enabled; ran the WL Ising and oscillator binaries. All exited
  successfully. The state example now uses fixed seeds and exact reversible moves.

These counts include deterministic, seeded stochastic, integration and compile-fail
doctests; they do not imply convergence for arbitrary models. Default and optional
feature builds keep the same four required `State` methods. `Params: Clone` is the
explicit cache rollback migration described in the README.
