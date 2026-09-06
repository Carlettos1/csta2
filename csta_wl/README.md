# CSTA Wang–Landau

Wang–Landau preliminary refinement, SAMC production sampling, and parallel
replica-exchange energy windows. The implementation lives in this directory;
`../cstav2/csta_wl` exposes it to the workspace as `csta::wl`.

## References and algorithm choices

- [Moreno, Peralta and Davis, arXiv:2103.15028v2](https://arxiv.org/abs/2103.15028v2),
  section 2 and Eq. (4): overlapping windows, frozen-DOS walkers after local
  completion, replica exchanges, and joining at the smallest entropy-slope
  mismatch. The supplied [PDF](2103.15028v2.pdf) and `rwl/window.py`, `rwl/wl.py`,
  `rwl/src/simulator.cpp`, and `rwl/src/qho.cpp` were used as references.
- [Shakirov, arXiv:2402.05653v2](https://arxiv.org/abs/2402.05653v2), Eqs. (1), (3),
  and (5): bin-DOS acceptance and the SAMC update `gamma(t) = t0 / (t1 + t)`.
  The supplied [PDF](2402.05653v2.pdf) defines `t` as trial moves. The old Rust
  implementation incorrectly used visits divided by the number of bins here.

This implements the combination needed by this project, not every extension in
those papers. In particular, it uses constant DOS per continuous bin (as in the
SAMC paper), not the optional within-bin interpolation in `rwl`. It does not
implement the second paper's importance-sampling parameter reweighting.

Preliminary stages start with `initial_ln_f` (default 1), require both histogram
flatness and minimum visits, and halve the update between stages. Each completed
stage clears its histogram, preserving DOS, state, and lifetime counts. Production
starts from that state and DOS with its own trial counter at zero. By default
`t0` is the number of accessible bins and `t1 = 10*t0`; both are configurable.
Flatness is a warm-up criterion, not an accuracy certificate. `TargetReached`
means the requested production budget was completed, not proven convergence.

## Serial API

```rust
use csta_wl::{Config, RawWangLandauData, models::Ising, run};
use rand::{SeedableRng, rngs::StdRng};

let config = Config {
    preliminary_stages: 3,
    min_stage_steps: 100,
    visits_per_bin: 20,
    sampling_steps: 30_000,
    max_steps: 100_000,
    ..Config::default()
};
let result = run(
    Ising::<4>::new([1; 4])?,
    (),
    StdRng::seed_from_u64(42),
    RawWangLandauData::on_grid(Ising::<4>::grid()?)?,
    config,
    None, // or Some(&AtomicBool) for cancellation
)?;
if !result.is_complete() {
    return Err("sampling stopped before its target".into());
}
let mut data = result.data.process_data()?;
data.normalize_log_count(16_f64.ln())?;
let (mean_energy, second_moment) = data.energy_moments(0.5)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`run` accepts any `csta::State`; it does not require `Randomizable`. The supplied
state must be inside the accessible energy support. Alternatively,
`sample_in_support` tries at most the specified number of random initial states
and returns `InitialStateOutside` if none is valid. It does not perform an
unbounded search for a rare energy window.

Contracts:

- Proposals must be **symmetric**. Asymmetric proposals/Hastings corrections are
  not supported. `apply_change` and `revert_change` must be exact inverses.
- The physical Hamiltonian parameters stay fixed during sampling. `energy` may
  refresh caches in `Params`; `Params: Clone` provides an independent snapshot
  restored on rejection. A clone must not share mutable caches or model state.
  A cache in the state must be restored by `revert_change` itself.
- `sampling_steps` is the number of production trial moves; preliminary moves
  do not count toward it. `max_steps` bounds preliminary plus production moves.
  Zero production steps means no moves, including no preliminary stages. Zero
  preliminary stages explicitly skips warm-up. Other minima/budgets are positive.
- Each completed sampling proposal, including rejected and self-loop proposals,
  visits exactly one retained bin. Nonfinite energy or failed numerical updates
  revert the proposal and parameter snapshot without recording a visit.
- `bins()` and `visits()` describe the current stage. `lifetime_bins()` and
  `total_visits()` include preliminary and production visits. Diagnostics report
  phase, stages, acceptance, last update, production moves and stopping reason.
- Cancellation is checked between proposals. Budget exhaustion/cancellation return
  partial `RunResult` data with an explicit reason. Invalid model/data/numerical
  operations return `Err`. User model methods must terminate: Rust cannot
  preempt an infinite loop inside `energy`, `propose_change`, or a change method.

## Energy support and numerical data

`EnergyGrid::discrete` takes exact, finite, increasing energy levels. Matching is
exact: use a consistent energy unit (integer-valued energies are useful), or
choose an interval grid when the energies are genuinely continuous. Off-grid
proposals are rejected. One discrete level is valid.

`EnergyGrid::continuous(min, max, bins)` uses half-open intervals; the final bin
includes `max`. Endpoints and adjacent floating-point values are mapped using
stored edges. Thermodynamic quadrature uses the bin centers. Finite values
outside the interval are rejected, never clamped into endpoint bins. NaN and
infinities are errors.

`RawWangLandauData::with_support(grid, mask, initial_log_dos)` allows explicitly
inaccessible bins: masked bins must have `-infinity` log-DOS, all accessible bins
must have finite log-DOS. At least one bin must be accessible. Unvisited bins are
never automatically declared inaccessible. A zero histogram is not flat; with
fraction `p`, every accessible bin must have at least `p * mean` visits,
including equality. The default is 0.8.

Data fields are private and constructors validate their shape. Read-only getters
replace writable public vectors. `dos()` always means **ln(g)**. For intervals,
`g` is bin mass, including bin width, not a density that needs an extra width
factor during the partition sum. Reuse a previous estimate through a constructor
with fresh counters, rather than passing a finished run back into `run`.

The absolute log-DOS offset is undetermined by sampling. Use
`normalize_log_count(ln(total_states))` or `anchor(bin, ln(degeneracy))` before
using absolute free energy/entropy. Normalization preserves relative weights.
Unrepresentable positive updates return a numerical error; anchor an excessively
large initial log-DOS offset before restarting.

Thermodynamic convention: **beta = 1/(k_b T)**.

- `C = k_b * beta² * Var(E)` and `F = -ln(Z)/beta`.
- Finite negative beta/temperature are supported for the finite supplied spectrum.
- Beta zero is valid for distributions/moments/heat capacity, but not free energy.
  Entropy's temperature argument must be finite and nonzero; `k_b` must be finite
  and positive. The caller must use matching beta and temperature.
- Probabilities use shifted logarithms; variance uses centered differences instead
  of subtracting nearly equal moments. Unrepresentable results return errors.
- Empty/all-absent log-sum-exp is `-infinity`; `+infinity` dominates finite values;
  NaN is an error. Thermodynamic data requires at least one finite mass.
- Microcanonical temperature requires three levels. Endpoints and stencils
  containing absent levels return `None`; zero entropy slope gives `+infinity`.
  Negative slopes give negative temperature. Nonuniform spacing is supported.

The oscillator example uses `E = sum(n_i + 1/2)`, in units `hbar*omega = 1`.
A proposed decrement at zero becomes a self-loop instead of a forced increment,
which keeps nontrivial transitions symmetric. Its energy grid imposes a finite
**total-quanta cutoff**; exact reference values use the same truncated spectrum.
The Ising reference is a periodic 1D chain with `J = 1` and at least three spins.
Both models cache energy, with exact apply/revert updates.

## Parallel API

See [examples/parallel_ising.rs](examples/parallel_ising.rs).

`run_parallel` (also exported as `par_wl`) takes a shared global grid, a vector of
initial states, common Hamiltonian parameters, `ParallelConfig`, and an optional
cancellation flag. Windows are ordered half-open **bin-index ranges**, beginning
at zero and collectively covering the grid. Adjacent windows must extend coverage
and overlap in at least two bins. There is no automatically guessed overlap or
binning. Use a discrete grid containing only accessible levels for sparse spectra.

Each window has `walkers_per_window` independent replicas, initialized in
window-major order. Each replica has its own state, parameter cache, RNG, DOS,
histogram and adaptive schedule. Replica log-DOS estimates are anchored to a
common bin and averaged within each window after completion.

Rayon executes bounded chunks in parallel. Alternating disjoint adjacent-window
pairs attempt exchanges, pairing equal replica indices. Both states must be
inside the destination supports. The log acceptance ratio is

```text
ln g_i(E_i) - ln g_i(E_j) + ln g_j(E_j) - ln g_j(E_i)
```

State and parameter caches move together; DOS, RNG and adaptation stay in their
windows. Finished replicas continue moving with frozen DOS while others finish.
Those moves appear only in `mixing_steps`. Exchange trials and successes have
separate counters and do not add histogram visits or change the SAMC clock.
This explicit counting convention differs from the extra histogram updates at
exchange in the supplied `rwl` code.

DOS fragments are joined at the overlap point with the closest entropy
slope, aligned there, then pasted. Endpoint slopes use secants; interior slopes
use three-point derivatives. Every accessible level must have been visited in
every contributing replica. Missing overlap data or incompatible grids returns
an error rather than a fabricated merged estimate. Merged histograms sum actual
current-stage visits from all replicas/windows, independent of the paste point.

`ParallelResult.merged` is present only after all replicas reach a nonzero target
and joining succeeds. Cancellation or any exhausted replica budget aborts the
whole run without a merged DOS. Worker errors/panics propagate as errors and
signal peer workers to stop at their next proposal boundary. No model call is
forcibly interrupted. Completion/join failures never masquerade as a successful
partial merge.

With the same seed, initial states, parameters, configuration, platform and
`rand` version, successful results are independent of Rayon thread count. Seeds
are assigned before scheduling, and exchanges happen in a fixed order. This does
not promise identical trajectories to the old implementation or reproducible
partial results under asynchronous external cancellation.

## API migration and workspace integration

This is a checked API revision, not source-compatible with the old prototype.
Package versions have not been released or bumped to v3.

| Old use | Replacement |
| --- | --- |
| `wang_landau(...) -> RawWangLandauData` | Same-shaped wrapper returns `Result<RunResult<S>>`; inspect `is_complete()` before processing `result.data`. |
| `wang_landau2(...)` | Same-shaped checked wrapper; prefer `run` + `Config` for explicit support and budgets. |
| `par_wl(target_time, ...)` | `run_parallel(grid, states, params, ParallelConfig, cancel)`. |
| `data.dos`, `data.bins`, `min/max/bin_width` | Read-only `dos()`, `bins()`, and `grid().energies()`. Use validated grid constructors. |
| `WLData::from(raw)` | `raw.process_data()?`. |
| Numerical methods returning bare floats | Return `Result`; microcanonical temperatures are `Vec<Option<f64>>`. |
| Full-energy scan in example models | `models::Ising` and `models::Oscillators`, with cached energy and exact references. |

The serial wrappers interpret `target_time` as production visits per bin, rounded
up to a whole proposal count. Their default total proposal budget is 10 million;
use `run` to choose a larger explicit bound. `wang_landau` uses two preliminary
stages, corresponding to the old halvings from 1 through 0.5 to 0.25.

The dependency graph is acyclic:

```text
csta -> csta_wl -> csta_metropolis / csta_montecarlo
```

Path dependencies ensure CSTA and WL share the same `State`/`Randomizable` traits.
`csta::wl` and `csta::prelude::wl` expose the module. The workspace shim avoids a
second copy of the algorithm; this is the integration point for the planned v3.
The sibling `wang-landau` copy and the validated social-model use cases are not
modified by this implementation.

## Verification

From this directory:

```sh
cargo test --offline
cargo clippy --offline --all-targets -- -D warnings
cargo test --offline --release -- --ignored
cargo run --offline --release
cargo run --offline --release -- qho
cargo run --offline --release --example parallel_ising
cargo run --offline --release --example profile_energy
cargo test --offline --manifest-path ../cstav2/Cargo.toml --workspace
```

Tests cover deterministic transition/grid/counting failures, exact DOS
thermodynamics, serial and parallel seeded convergence, replica stitching,
thread-count reproducibility, cancellation and worker failures. The larger Ising
study is explicitly ignored in normal runs. Monte Carlo tolerance tests use
fixed seeds and fixed budgets; they are empirical regression checks, not proofs
of convergence for arbitrary models.

The energy microbenchmark performs 20,000 queries on a 4,096-spin chain. A local
release run measured approximately 55.7 ms for full recomputation versus 14.6 µs
for cached queries. This measures energy queries only, not whole-simulation
speedup; results depend on hardware/compiler. Full histogram scans were also
removed from production updates; per-proposal grid lookup is logarithmic and
recording is constant-time.
