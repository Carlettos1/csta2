# Statistical mechanics and computational physics roadmap

This roadmap is now implemented in the workspace. The original requirements are
retained below with implementation pointers. Priorities are **should** (correctness and reliable everyday science),
**could** (useful extensions), and **neat** (optional exploratory capabilities).
Every item has corresponding average and edge cases in [TEST_PLAN.md](TEST_PLAN.md).
The optional checkpoint feature is the only additional feature needed for persistence;
ordinary sampling keeps empty default features. See the README for migration details.

## Keep usage simple

- Preserve the short path: construct a model, choose a sampler and budget, run,
  inspect observables. No new mandatory traits or configuration for existing models.
- Prefer methods, small result structs, closures, and defaults over additional
  frameworks. Advanced samplers and persistence should be opt-in.
- State units, ensemble, normalization, and statistical assumptions in docs.
  Defaults may simplify configuration but must not guess physical constraints.
- Reuse existing RNG injection, `State`, observers, grids, and checked WL results.
  Avoid duplicating features: WL already provides DOS thermodynamics, cancellation,
  explicit budgets, reference models, and parallel energy-window exchange.
- Each addition needs a short runnable example and a test showing an existing
  simple workflow still compiles. Correctness fixes may require a documented
  migration, but must not demand extra setup for ordinary valid inputs.

## Should have

- [x] **S01 — Reliable vector arithmetic (`csta_core`).** Fix precision conversions:
  the original `impl_from` macros copied component zero into every component.
  Reconcile `/` and `/=` behavior for exceptional divisors, and make norms and
  normalization robust for very small/large finite vectors. Document zero and
  nonfinite behavior. Keep existing operators; add a checked normalization helper
  only if necessary. These are foundations for displacements, forces, and spins.

- [x] **S02 — Checked canonical Metropolis transitions (`csta_metropolis`).**
  Validate numerical inputs and evaluate acceptance with a stable log ratio.
  The original unconditional `delta_energy < 0` branch was incorrect for negative
  beta; explicitly support finite negative beta for bounded spectra or reject it
  clearly. Roll back parameter caches as well as state on rejection/error, as WL
  already does. Keep the usual constructor/run workflow; avoid imposing a new
  required model trait. Resolve the cache snapshot requirement explicitly.

- [x] **S03 — Trustworthy run counts and measurement timing (`csta_metropolis`).**
  Count actual attempted moves rather than divide cumulative acceptances by the
  configured `steps`. Define zero-attempt rates, repeated-run behavior, burn-in,
  pre/post-move observations, and measurement stride. Original observers used
  `i > after()` and `i % every()`, with no zero-stride validation. Add bounded,
  cancellable runs with explicit completion status, following the WL convention.
  Offer one default measurement schedule and retain convenient observer usage.

- [x] **S04 — Streaming observables and uncertainty (`csta_metropolis`, `csta`).**
  Add reusable running mean, variance, covariance, and block statistics so users
  need not retain every sample. Include autocorrelation-aware standard errors,
  an explicit autocorrelation-time convention, and effective sample size.
  Report insufficient data instead of fabricated precision. Expose a small
  summary object and closure-based measurement adapter with sensible defaults.

- [x] **S05 — Standard thermodynamic measurements (`csta`, `csta_metropolis`).**
  Supply energy/magnetization moment summaries, heat capacity, susceptibility,
  and Binder cumulants. WL already computes energy thermodynamics; reuse its
  conventions. Distinguish total from per-site quantities and signed M from |M|.
  A single summary call should handle the common zero-field spin-model workflow,
  with explicit site count and units rather than mandatory custom observers.

- [x] **S06 — Validated lattice examples (`csta/examples`, reusable model module).**
  Repair `examples/ising.rs`: unsigned neighbor subtraction, row wrapping, and
  bond counting need explicit boundary conventions. Promote a tested 2D Ising
  model with coupling J, field h, and explicit open/periodic boundaries; keep the
  current WL 1D Ising and truncated oscillators as reference systems. Provide
  simple constructors and cached local energy updates with recomputation checks.

- [x] **S07 — Reliable random initialization and derive errors (`csta_derive`,
  `csta_montecarlo`).** Weighted enum sampling already exists; validate weights
  and reject malformed inputs instead of macro panics or silent fallback.
  Cover ranges, lengths, defaults, generics, and field dependencies. Provide
  optional named isotropic-direction and Gaussian sampling helpers for continuous
  spins and particle initialization. Do not change existing uniform component
  sampling or suggest that uniform cube samples are isotropic directions.

- [x] **S08 — Reproducible physics examples and documentation (`csta`, workspace).**
  Populate the empty root README and correct stale paths/version claims in the WL
  README. Show seeded canonical and DOS runs with units, burn-in, uncertainties,
  finite-size limitations, and explicit partial-result handling. Include a small
  exact-enumeration comparison. Keep each introductory example runnable with one
  command and use current facade exports.

## Could have

- [x] **C01 — Asymmetric proposals (Metropolis and WL).** Add an optional proposal
  adapter carrying log reverse/forward proposal probability for Hastings
  corrections. Existing symmetric `State` models must need no changes. Define
  zero reverse probability and constrained-boundary moves explicitly.

- [x] **C02 — Temperature replica exchange (`csta_metropolis`).** Complement WL's
  existing energy-window exchanges with a canonical beta ladder, independent
  replicas, and swap/round-trip diagnostics. Accept an explicit ladder and seed;
  provide a simple run entry point. Keep caches with states and temperatures with
  replicas, with deterministic scheduling under a fixed execution configuration.

- [x] **C03 — Cluster moves for supported spin models.** Offer an opt-in Wolff
  sampler for ferromagnetic, zero-field Ising systems to study slow mixing near
  transitions. Reuse the validated lattice model. Reject unsupported couplings
  or fields rather than silently applying the wrong algorithm. Report cluster
  size and distinguish cluster updates from attempted single-spin moves.

- [x] **C04 — Observable reweighting and overlap diagnostics (`csta_wl`, analysis).**
  Extend existing energy-only DOS thermodynamics with conditional observable
  moments by energy bin; optionally add canonical single-histogram reweighting.
  Expose temperature evaluation as a method on collected data. Flag missing
  support and poor overlap; energy DOS alone cannot reconstruct magnetization.

- [x] **C05 — Independent-run uncertainty for DOS (`csta_wl`).** Summarize multiple
  independently seeded DOS estimates after consistent normalization; propagate
  uncertainty to observables and expose disagreement in overlaps. A collection
  of results should produce one summary. Do not treat interacting replicas or
  histogram flatness as independent evidence of convergence.

- [x] **C06 — Checkpoint and resume (samplers).** Add opt-in versioned snapshots
  containing state, parameter caches, RNG state, adaptation phase, and all clocks.
  Saving/loading should take one call for supported serializable models and RNGs.
  Distinguish exact continuation from restarting with a previous DOS estimate;
  reject incompatible snapshots and document cancellation/resume semantics.

- [x] **C07 — Geometry and local-energy helpers (`csta_core`, model utilities).**
  Add orthorhombic periodic wrapping, minimum-image displacement, and reusable
  lattice neighbor access. Offer optional local delta-energy evaluation with the
  existing full-energy calculation as fallback. Keep box conventions explicit
  and check cached results against full recomputation. Avoid requiring callers
  to manage caches or neighbor bookkeeping themselves.

## Would be neat to have

- [x] **N01 — Finite-size analysis recipes.** Build on S04–S06 to compare Binder
  curves, susceptibility peaks, and correlation estimates across lattice sizes.
  Start with a runnable example and tabular output, not a mandatory plotting or
  fitting dependency. Expose uncertainty and fit-range sensitivity; do not label
  a finite-system peak as an exact critical temperature.

- [x] **N02 — Joint density of states.** Offer an opt-in (E, M) grid for magnetic
  field reweighting, with one result method accepting beta and field. Keep 1D DOS
  unchanged; validate memory limits and sparse support. Define E as the field-free
  energy so the field term is not counted twice.

- [x] **N03 — Minimal molecular dynamics example/module.** Provide velocity-Verlet
  integration for a harmonic system, then a small pair-potential example using
  C07. Start with fixed-step NVE trajectories and energy-drift diagnostics.
  Use a small separate position/velocity/force API; stochastic `State` users
  should not acquire force or velocity requirements.

- [x] **N04 — Grand-canonical lattice gas.** Add an optional occupancy model with
  insertion/deletion moves and chemical potential, reusing C01 where proposals
  are asymmetric. Expose a simple `(beta, mu, budget)` workflow, document the
  lattice-gas ensemble, and retain fixed-particle models unchanged.

Implementation order: correctness foundations first, then measurements and models,
followed by optional samplers, analysis, checkpoints, and dynamics.

## Implementation map

| Items | Implementation and example |
| --- | --- |
| S01 | `csta_core/src/vec{2,3,4}.rs`: coordinate casts, IEEE operators, scaled normalization. |
| S02, S03, C01 | `csta_metropolis/src/lib.rs`: checked transitions, cumulative counts, schedules, fallible closures, `Hastings`; corrections also used by WL. `csta/examples/ising.rs` and `hastings.rs`. |
| S04, S05 | `csta_metropolis/src/statistics.rs`: moments/covariance, block and correlation errors, thermal summaries and block jackknife. `csta/examples/ising.rs`. |
| S06, C03, N04 | `csta_metropolis/src/models.rs`: validated Ising2D, Wolff, lattice gas with symmetric site toggles. `csta/examples/physics.rs`. |
| S07 | `csta_derive/src/lib.rs`, `csta_montecarlo/src/lib.rs`: checked attributes, dependency order, Gaussian and isotropic sampling. `csta/examples/rand_derive.rs` and `physics.rs`. |
| S08 | Root README doctests, corrected WL README and bounded seeded examples. |
| C02 | `csta_metropolis/src/replica.rs`: temperature ladder, exchanges, round trips, resumable chunks. `csta/examples/physics.rs`. |
| C04, C05 | `csta_wl/src/analysis.rs`, `session.rs`: conditional moments, histogram reweighting, independent DOS uncertainty, observed production. `csta/examples/reweight.rs`. |
| C06 | `csta_core/src/checkpoint.rs`, sampler save/load methods, `csta_wl/src/session.rs` and `par_wl.rs`: versioned atomic snapshots, PCG RNG, exact continuation. `csta/examples/checkpoint.rs`. |
| C07 | `csta_core/src/geometry.rs`, optional `State::delta_energy`, validated cached model updates. |
| N01 | `csta_metropolis/src/finite_size.rs`: explicit crossing/peak ambiguity and block bootstrap; `csta/examples/finite_size.rs` reports thermal and correlation estimates across sizes. |
| N02 | `csta_wl/src/joint.rs`: sparse declared support, joint WL/SAMC and field reweighting. `csta/examples/physics.rs`. |
| N03 | `csta_core/src/dynamics.rs`: velocity-Verlet, harmonic and tiny pair examples, energy drift. `csta/examples/physics.rs`. |

Design decisions: `Params: Clone` is the explicit rollback requirement, shared
with WL; existing model trait methods remain sufficient. Negative beta is allowed
with caller-selected normalizable spectra. Error estimates report insufficient
data. Independence of separate experiments remains a caller assertion; IDs catch
accidental reuse. Reweighting reports weight concentration and cannot certify
unobserved support. Finite-size tools do not fit or claim an exact critical point.
Checkpoints are synchronous atomic writes and can save a cancelled sampler;
resumption preserves the documented proposal/exchange clocks.
