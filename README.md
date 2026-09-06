# CSTA

Rust tools for statistical mechanics and small computational-physics experiments.
The workspace uses Rust 1.98 or newer, edition 2024. Start with the `csta` facade;
advanced functionality is available through its modules. Package version: 3.0.0
(workspace metadata; this does not imply a crates.io release).

The [library design document](docs/design.pdf) covers the full workspace, with
code examples, physical conventions, supported ranges, and limitations.
Its [LaTeX source and build instructions](docs/README.md) include a checker for
the embedded Rust examples.

## Canonical sampling

```rust
use csta::{Metropolis, Schedule, State, models::Ising2D, statistics::Blocking};
use rand::{SeedableRng, rngs::StdRng};

let model = Ising2D::aligned(4)?; // periodic 4x4, J=1, h=0
let mut run = Metropolis::with_all(model, (), 0.4, 10_000, StdRng::seed_from_u64(42));
let mut energy = Blocking::new(128)?;
let report = run.try_run_observed(
    Schedule { burn_in: 1000, stride: 1 },
    None, // optional cancellation flag
    |s, _| energy.push(s.recompute_energy() / 16.0),
)?;
assert_eq!(report.status, csta::RunStatus::Complete);
println!("energy/site: {:?}", energy.estimate());
# Ok::<(), Box<dyn std::error::Error>>(())
```

`beta = 1/(k_b T)`; examples use `k_b=1`. Finite negative beta is supported only
when the caller's physical spectrum admits a normalizable negative-temperature
ensemble. `State` retains its four required methods. Optional proposal-ratio,
local-energy, and checkpoint-validation methods have defaults. Rejected and
invalid moves restore the state and an independent clone of `Params`.

Measurements occur **after** completed moves. Burn-in skips that many moves;
then the first measurement occurs after `stride` more moves. Schedules restart
for each run, while acceptance counters accumulate actual attempts. Zero-attempt
acceptance and rejection rates are both zero. `run_with`, `run_with_2/3/4`, and
dynamic observers use this same convention. Use `try_step`/`run` to handle errors;
legacy `step`/`run_empty` convenience methods panic on invalid input.

`statistics::Thermodynamics` accepts total energy and signed total magnetization;
its summary reports energy, magnetization, heat capacity, and susceptibility
**per site**, plus the dimensionless Binder cumulant. Supply the site count and
`k_b` explicitly. Undefined Binder ratios return `None`. `Blocking` uses constant
memory; `correlated_estimate` optionally analyzes retained scalar samples.
`thermal_errors` uses a delete-one-block jackknife for nonlinear thermodynamics.
Blocks must exceed correlation times. Short, constant, or unresolved series do
not establish mixing or justify precise error bars. Statistical errors do not
include bias from incomplete equilibration or finite system size.

## DOS sampling

The checked WL/SAMC API, reference 1D Ising/oscillator models, support contracts,
and thermodynamics are documented in [csta_wl/README.md](csta_wl/README.md).

```rust
use csta::wl::{self, Config, RawWangLandauData, models::Ising};
use rand::{SeedableRng, rngs::StdRng};
let result = wl::run(
    Ising::<4>::new([1; 4])?, (), StdRng::seed_from_u64(42),
    RawWangLandauData::on_grid(Ising::<4>::grid()?)?,
    Config { preliminary_stages: 0, sampling_steps: 1000, ..Config::default() },
    None,
)?;
if !result.is_complete() { return Err("partial DOS run".into()); }
let mut dos = result.data.process_data()?;
dos.normalize_log_count(16_f64.ln())?;
println!("mean energy: {}", dos.energy_moments(0.5)?.0);
# Ok::<(), Box<dyn std::error::Error>>(())
```

`dos()` is **ln(g)**; continuous bins contain probability mass, not density per
unit energy. Normalize/anchor before absolute free energies or entropies. A
completed budget or flat histogram is not proof of convergence.

## Optional tools

| Module/API | Purpose and explicit limits |
| --- | --- |
| `models::Ising2D` | Open/periodic square lattices, J and h, reversible cached flips. Periodic dimensions must be at least 3; open dimensions may be 1. |
| `Ising2D::wolff_step` | One cluster update for J>=0, beta>=0, h=0; returns cluster size. |
| `models::LatticeGas` | Uniform site toggles with effective energy H-mu*N; symmetric grand-canonical proposals. `physical_energy` excludes chemical potential. |
| `Hastings` / `State::log_proposal_ratio` | Opt-in log reverse/forward correction, accepted by canonical and WL sampling. Zero reverse probability rejects. |
| `replica::ReplicaExchange` | Explicit increasing beta ladder, bounded local chunks, swap and round-trip counts. Sequential, deterministic schedule. |
| `statistics`, `finite_size` | Stable moments/covariance, block and correlation errors, Binder curves, crossings, peaks, block bootstrap. |
| `wl::analysis` | Conditional observable moments, canonical histogram reweighting, independent DOS-run uncertainty. Weight ESS cannot detect unseen states. |
| `wl::joint` | Explicit sparse (E0,M) support, WL sampling and field reweighting; E0 excludes the magnetic field. Caller supplies a cell budget. |
| `geometry` | Orthorhombic periodic wrapping, minimum images, unique lattice bonds. Half-box ties map to the negative side. |
| `dynamics` | Separate fixed-step NVE velocity-Verlet, harmonic and tiny nearest-image Lennard-Jones examples; no thermostat, cutoff, tail correction, or bulk-fluid accuracy claim. |
| `gaussian`, `isotropic_direction` | Explicit normal sampling and uniform directions on S²; existing component-uniform sampling is unchanged. |

## Reproducibility and checkpoints

Enable `checkpoint` only for persistence. It supplies a serializable `CheckpointRng`
(PCG64, pinned rand_pcg 0.10.2); ordinary workflows keep their existing RNGs.
`Metropolis`, `wl::Session`, `wl::ParallelSession`, and `replica::ReplicaExchange`
provide `save(path, model_tag)`
and `load(path, model_tag)`. Model tags must change when the model/Hamiltonian
schema changes. Models/parameters/RNGs must implement Serde; model parameters
must clone independently. Built-in canonical lattice models and WL reference models support this feature.
Use `ReplicaExchange::<_, CheckpointRng>::seeded(...)` for temperature-exchange
checkpoints; resume a cancelled chunk with its original `local_steps` value.
Custom persisted models should override `State::valid_state` to validate their
structural invariants. Snapshot loaders validate sampler/grid/counter invariants.

A WL session's `advance` retains its original warm-up, target, total budget, RNG,
and DOS. Serial calls count proposals; parallel calls count exchange chunks.
`Session::advance_observed` records retained production states through a closure,
so conditional moments need no custom observer trait.
Snapshots taken between successful parallel calls are at an exchange barrier.
Cancellation may be resumed; exhausting the original total budget cannot be
silently extended. Finishing an unfinished session returns an explicit partial
result. A previous DOS with fresh counters is a new run, not exact continuation.

Files include a format, implementation, Rust type/RNG and model tag. Saving uses
a sibling temporary file and rename; invalid inputs leave an existing snapshot
intact. Loading is limited to 256 MiB. Checkpoints preserve trajectories for the
same implementation, dependencies and platform; user callbacks and asynchronous
cancellation still determine their own behavior. Use trusted local snapshots.

## Runnable examples

See [EXAMPLES.md](EXAMPLES.md) for twelve physics problems, including five
standalone programs for paramagnetism, adsorption, diffusion, harmonic-chain
vibrations, and quantum-oscillator DOS thermodynamics.

Run from the workspace root:

```sh
cargo run -p csta --release --example state
cargo run -p csta --release --example ising
cargo run -p csta --release --example physics
cargo run -p csta --release --example reweight
cargo run -p csta --release --example hastings
cargo run -p csta --release --example finite_size
cargo run -p csta --features checkpoint --example checkpoint
cargo run -p csta --release --example parallel_ising
cargo run -p csta --release --example profile_energy
cargo run -p csta_wl --release
cargo run -p csta_wl --release -- qho
```

## Validation and migration

```sh
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features
cargo test --workspace
cargo test --workspace --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test -p csta_wl --release -- --ignored
```

Use `--offline` when dependencies are cached. See [TEST_PLAN.md](TEST_PLAN.md) for
coverage and [TODO.md](TODO.md) for the implemented roadmap.

Migration from the prototype: Metropolis execution now requires `Params: Clone`
for rollback; deriving Clone for a simple parameter struct is sufficient. Struct
literals for `Metropolis` need `attempted_moves` (constructors fill it). Observer
timing now follows the documented post-move schedule. Vector casts preserve all
coordinates; `/=` follows scalar IEEE division, matching `/`. `normalize` returns
zero for undefined directions; `try_normalize` returns `None`. Derive attributes
now reject invalid weights/ranges/conflicts with compiler diagnostics. Named-field
initializers are evaluated in dependency order; cyclic dependencies are errors.
Use expressions referring to fields directly and avoid shadowing field names in
initializer closures. Advanced APIs are optional and default features stay empty.

Algorithm references: [Geyer's initial-sequence estimators](https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html),
[Ising cluster algorithms](https://www.physik.uni-leipzig.de/~janke/Paper/lviv-ising-lecture-janke_corrected.pdf),
and [velocity-Verlet in LAMMPS](https://docs.lammps.org/fix_nve.html).
The correlation implementation uses monotone positive pairs and conservatively
floors tau at 1; it reports unavailable when the lag budget cannot resolve a cutoff.
