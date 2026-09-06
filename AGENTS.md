# Repository guidance

## Project layout

This is a Rust workspace for statistical mechanics simulations. The root
`Cargo.toml` defines workspace members and shared package metadata, including
Rust edition 2024 and minimum Rust version 1.98.0.

- `csta`: public facade and prelude, integration tests in `csta/tests`, and
  runnable examples in `csta/examples`.
- `csta_core`: two-, three-, and four-dimensional floating-point vectors, with
  optional Serde support.
- `csta_derive`: the `Randomizable` procedural derive macro.
- `csta_montecarlo`: the `Randomizable` trait and Monte Carlo iterator.
- `csta_metropolis`: the `State` trait, Metropolis sampling, and observers.
- `csta_wl`: Wang–Landau/SAMC sampling, energy grids, thermodynamic data,
  parallel replica exchange, reference models, and a demonstration binary.

Keep implementations in their owning crates. Preserve the dependency direction
from the `csta` facade to the algorithm and core crates; do not introduce runtime
dependencies from those crates back to `csta`. When changing public exports,
check both `csta/src/lib.rs` and `csta/src/prelude.rs`.

## Development and validation

Run commands from the repository root. Use a toolchain satisfying the manifest's
minimum Rust version; do not lower that requirement or change dependencies merely
to accommodate a local toolchain.

```sh
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features
cargo test --workspace
cargo test --workspace --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
```

Choose checks appropriate to the change. Start with `cargo test -p <crate>` for
localized behavior changes, and use workspace checks for shared traits, macros,
dependencies, or public API changes. Exercise Serde with `--all-features` when
changing vectors or feature wiring. Report any existing failures or unavailable
toolchain/dependencies accurately. Documentation-only edits do not require the
simulation test suite.

Useful targeted commands:

```sh
cargo test -p csta --test wang_landau
cargo test -p csta_wl --release -- --ignored
cargo run -p csta_wl --release
cargo run -p csta_wl --release -- qho
cargo run -p csta --release --example parallel_ising
cargo run -p csta --release --example profile_energy
```

Ignored release tests and profiling examples are deliberate, potentially longer
checks; run them when the algorithm or performance change warrants it. Add
`--offline` only when the required dependencies are already cached.

## Implementation conventions

- Follow the surrounding Rust style and use rustfmt. Keep changes focused and
  avoid unrelated formatting churn.
- Preserve shared workspace metadata and the existing path-dependency setup.
  Update `Cargo.lock` when dependency changes require it; do not edit generated
  files under `target/`.
- Put regression tests beside the affected implementation, following existing
  inline test modules or `csta_wl/src/tests.rs`. Use `csta/tests` to verify facade
  integration. Check example compilation when changing traits or derive output.
- For stochastic regression tests, use fixed seeds and bounded proposal budgets.
  Prefer exact small-model reference values and explicit numerical tolerances.
  A successful finite sampling run is not proof of convergence.

## Sampling and numerical contracts

Consult `csta_wl/README.md` and the implementation before changing sampling
behavior. Resolve example paths against this workspace and use the package-qualified
commands above.

- Wang–Landau proposals default to symmetric; asymmetric models must supply
  the optional log reverse/forward proposal ratio. Applying/reverting a change
  must restore state exactly. Rejections must also restore parameter caches; cloned
  parameters must provide independent mutable snapshots.
- Keep preliminary, production, lifetime, exchange, and frozen-DOS mixing
  counters distinct. Each completed sampling proposal, including a rejection
  or self-loop, records one retained-bin visit.
- Validate energy grids and support masks. Reject finite off-support proposals;
  treat nonfinite energies as errors. Do not silently clamp energies or infer
  that an unvisited bin is inaccessible.
- DOS values represent `ln(g)`; continuous-bin values represent bin mass.
  Preserve stable logarithmic calculations and explicit numerical errors.
- Preserve bounded execution, cancellation checks, and explicit partial-result
  stop reasons. Parallel results must not expose an apparently successful
  merged DOS when sampling or joining is incomplete.
- Preserve seeded parallel reproducibility across Rayon thread counts under
  the existing documented conditions.
