# Physics problems with CSTA

These twelve problems connect physical questions to CSTA's current APIs. Problems
1–5 have new, standalone implementations in `csta/examples`; the remaining
problems are recipes, with links to related existing examples where available.

The five programs need no arguments, new dependencies, or optional features.
Run from the workspace root with Rust 1.98 or newer:

```sh
cargo run -p csta --release --example paramagnet
cargo run -p csta --release --example lattice_gas_adsorption
cargo run -p csta --release --example brownian_diffusion
cargo run -p csta --release --example harmonic_chain
cargo run -p csta --release --example oscillator_dos
```

Each writes CSV to standard output and explanatory notes to standard error.
For example, append `> adsorption.csv` to save the adsorption table. Seeds and
budgets are fixed in the source. Edit the named constants to change the system
size or resolution; parameter sweeps are explicit arrays in `main`.

All thermodynamic examples use `k_B = 1` and `beta = 1/T`. Energies supplied to
samplers are totals; output labels specify any per-spin or per-oscillator division.
Monte Carlo proposals are not physical time steps. Numerical completion does not
establish equilibration or convergence. Statistical errors exclude model,
finite-size, and cutoff bias.

| Implemented problem | Main CSTA components | Default work |
| --- | --- | --- |
| 1. Paramagnet | Custom `State`, `Metropolis`, thermodynamic jackknife | 32 spins; 6 beta values; 20,000 burn-in + 256,000 retained proposals each |
| 2. Adsorption | `LatticeGas`, `Metropolis`, `Blocking` | 36 sites; 5 chemical potentials; 20,000 burn-in + 256,000 retained proposals each |
| 3. Diffusion | `Vec3f64`, `gaussian`, `Moments` | 4,096 trajectories, 400 increments each |
| 4. Harmonic chain | `Dynamics<1>`, custom force callback | 16 masses; 500, 1,000, 2,000, and 4,000 steps over the same duration |
| 5. Oscillator DOS | WL/SAMC, `Oscillators`, `DosEnsemble` | 3 oscillators; 33 energy levels; 4 runs of 250,000 production proposals, each capped at 1,000,000 total proposals |

## 1. Magnetization and the Schottky heat-capacity peak

**Question.** How do independent two-state magnetic moments polarize and store
heat as temperature decreases?

Use `N` spins `s_i = ±1` in a field `h`, with

```text
H = -h sum_i s_i
M = sum_i s_i
```

Implement the four required `State` methods. The example caches integer total
magnetization and proposes a uniformly selected site's replacement by either
spin value with equal probability. This symmetric proposal includes self-loops,
so it remains aperiodic at `beta = 0`. Changes preserve the old spin and revert
both the spin and the cache exactly.

Measure total `(E, M)` after every retained proposal using `Thermodynamics`.
Compare per-spin values with the exact independent-spin state sum:

```text
Z = [2 cosh(beta h)]^N
<M>/N = tanh(beta h)
<E>/N = -h tanh(beta h)
C/N = (beta h)^2 / cosh(beta h)^2
```

**Implementation:** [paramagnet.rs](csta/examples/paramagnet.rs).
The CSV includes sampled and exact energy, magnetization, and heat capacity,
plus block-jackknife errors for energy and heat capacity. Blocks contain 1,024
retained proposals. `NA` means the uncertainty estimator could not resolve an
estimate; it does not mean zero uncertainty.

**Expected result.** Magnetization approaches saturation at large positive
`beta h`; heat capacity vanishes at both small and large `beta h`, with a peak
between them. There is no inter-spin interaction or cooperative phase transition.
At low temperature, rare excited configurations require longer runs. Compare
several seeds and block lengths before interpreting the error bars.

## 2. An adsorption isotherm on independent lattice sites

**Question.** How does equilibrium coverage respond to a particle reservoir's
chemical potential, and how is that response related to number fluctuations?

Use `models::LatticeGas` with `J = 0` on a periodic 6×6 lattice. Each site has
occupation `n_i = 0` or `1`. The physical interaction energy is zero; the sampler
uses `H_eff = H - mu N_p = -mu N_p`. Uniform site toggles are symmetric proposals
on this discrete configuration space.

Sweep `mu` at fixed `beta = 1`. Measure particle number after retained proposals,
including rejections, and compare with the independent-site result:

```text
rho = <N_p>/N_sites = 1 / [1 + exp(-beta mu)]
Var(N_p) = N_sites rho (1-rho)
d rho/d mu = beta Var(N_p)/N_sites
```

**Implementation:**
[lattice_gas_adsorption.rs](csta/examples/lattice_gas_adsorption.rs).
It prints density, a block standard error, particle-number variance, and the
fluctuation-derived response alongside their exact values.

**Expected result.** Coverage rises smoothly from almost empty to almost full,
with half filling at `mu = 0` and the largest response near half filling. The
finite-run fluctuation and response estimates have no error bars in this example.
Equilibrium sites are independent, but consecutive local Monte Carlo observations
are correlated. Adding `J != 0` introduces interactions and invalidates the
independent-site formula. This is lattice adsorption, not continuum particle
insertion; continuum volume and proposal factors would require another model.

## 3. Recovering a diffusion coefficient from Brownian trajectories

**Question.** Can an ensemble of random paths recover its imposed diffusion
coefficient and the linear growth of mean-square displacement?

Generate unbounded three-dimensional paths with Gaussian increments:

```text
r(t+dt) = r(t) + sqrt(2 D dt) * xi
xi_x, xi_y, xi_z ~ independent standard normal
MSD(t) = <|r(t)-r(0)|^2> = 6 D t
D_estimate(t) = MSD(t)/(6t)
```

Use `gaussian` for increments, `Vec3f64` for positions, and one `Moments`
accumulator for each observation time. Retain one squared displacement per path
at each time; the standard error at a fixed time comes from the variation across
separately seeded paths. No trajectory history is needed.

**Implementation:** [brownian_diffusion.rs](csta/examples/brownian_diffusion.rs).
It uses `D = 0.7`, `dt = 0.02`, and reports five times from 0.5 to 8.0, including
MSD and diffusion estimates with ensemble standard errors.

**Expected result.** MSD grows linearly and `D_estimate` fluctuates around 0.7.
The rows share trajectories and are correlated; treating them as independent
observations in a line-fit uncertainty calculation would be incorrect. Positions
are unwrapped, avoiding artificial saturation from periodic coordinate wrapping.
This code composes random primitives for free overdamped diffusion; it has no
inertia, interactions, confinement, or general Langevin integrator. The MSD
approach is also described in the [LAMMPS diffusion documentation](https://docs.lammps.org/Howto_diffusion.html).

## 4. Normal modes of a periodic harmonic chain

**Question.** Does a numerical lattice vibration follow its exact normal mode,
and how does its error change when the integration step is halved?

Take a periodic ring of `N = 16` equal masses with displacements `u_i` from their
equilibrium sites. Neighboring displacements are coupled by springs:

```text
U = (k/2) sum_i (u_(i+1) - u_i)^2
F_i = k (u_(i-1) + u_(i+1) - 2u_i)
q = 2 pi mode/N
omega(q) = 2 sqrt(k/mass) |sin(q/2)|
u_i(t) = A cos(q i) cos(omega t)  [initial velocity zero]
```

Supply this potential and force as a callback to `dynamics::Dynamics<1>`.
Coordinates are displacements; periodicity is in spring connectivity, so wrapping
displacements with `PeriodicBox` would change the problem. Count each spring once
and add equal-and-opposite force contributions.

**Implementation:** [harmonic_chain.rs](csta/examples/harmonic_chain.rs).
It excites mode 1 with amplitude 0.1, `k = mass = 1`, and integrates to 1.3 exact
periods using four step counts. The output gives maximum final displacement
error, maximum absolute energy drift, total momentum, and consecutive error ratios.

**Expected result.** Momentum stays near zero. Halving `dt` approaches a fourfold
reduction in displacement error, consistent with second-order velocity-Verlet.
Tests compare forces with finite-difference potential gradients, check the rigid
translation mode and invalid inputs, and verify the analytic mode's error scaling.
The model has no anharmonic scattering or thermalization, and a small energy
drift alone is not a bound on long-time phase error.

## 5. Quantum-oscillator thermodynamics from a sampled DOS

**Question.** Can one estimated density of states reconstruct energy and heat
capacity at several temperatures, while separating sampling error from a spectral
cutoff?

Use `wl::models::Oscillators<3>` with `hbar*omega = 1`:

```text
E = Q + N_oscillators/2
Q = sum_i n_i, n_i >= 0
g(Q) = choose(Q + N_oscillators - 1, N_oscillators - 1)
```

Declare a grid with **total quanta** `Q <= 32`. Run two preliminary WL stages
followed by 250,000 SAMC production proposals, with a separate total bound. Reject
partial results and require a visit to every declared level. Normalize each run's
log-DOS to the exact number of states in that same truncated spectrum.

**Implementation:** [oscillator_dos.rs](csta/examples/oscillator_dos.rs).
Four separately seeded experiments form a `DosEnsemble`. For each beta, it prints
per-oscillator energy and heat capacity, their between-run standard errors, exact
truncated values from `Oscillators::exact_dos`, and infinite-spectrum values:

```text
<E>/N_oscillators = 1/2 + 1/[exp(beta)-1]
C/N_oscillators = beta^2 exp(-beta)/[1-exp(-beta)]^2
```

**Expected result.** Sampled values should approach the exact truncated state sum
as sampling improves. At high temperature, the finite cutoff visibly changes
energy and heat capacity relative to the infinite spectrum. Increasing the
production budget cannot remove that difference; increase `MAX_QUANTA` and repeat.
The four-run standard error does not include common DOS bias or cutoff bias.
`WLData::specific_heat` is a total heat capacity, so the program divides by the
oscillator count explicitly. Infinite-spectrum formulas here require positive
beta; zero and negative beta do not normalize that infinite oscillator ensemble.

## 6. Finite-size rounding near the two-dimensional Ising transition

**Question.** How do susceptibility peaks and Binder-cumulant crossings change
with lattice size?

Use periodic `models::Ising2D` at `J = 1, h = 0` for several sizes, such as 4, 8,
and 16. Sweep beta and sample with `wolff_step` or local `Metropolis` updates.
Record total energy and signed total magnetization. Use `Thermodynamics`,
`thermal_errors`, and `finite_size::{crossings, peak}` to construct and compare
curves.

**Deliverables.** A temperature table of per-site energy, susceptibility, Binder
cumulant, and block errors for each size; all crossing locations on a shared beta
grid. Compare local and cluster algorithms using their explicit time units.

**Related implementation:** [finite_size.rs](csta/examples/finite_size.rs).
This existing small-size demonstration is a starting point, not a fit of critical
exponents. A finite-size peak is not itself a thermodynamic-limit critical point.
Wolff requires nonnegative beta and coupling with zero field; a cluster update is
not one single-spin proposal. Periodic dimensions must be at least three.

## 7. A magnetic equation of state from joint density of states

**Question.** Can one finite-system calculation reconstruct magnetization for
many fields and temperatures?

For a small Ising model, enumerate accessible pairs `(E0, M)`, where `E0` excludes
the external field. Construct `wl::joint::JointGrid` from sorted unique pairs
with an explicit cell budget, then use `run_joint`. Field reweighting evaluates

```text
Z(beta,h) = sum_(E0,M) g(E0,M) exp[-beta (E0 - h M)]
```

**Deliverables.** Magnetization-versus-field curves, odd symmetry in `h` at finite
size, and comparison with direct state enumeration. Use separately seeded DOS
experiments to assess uncertainty.

**Related material:** the joint-DOS chapter and runnable snippets in
[the design document](docs/design.pdf), and exact joint-data construction in
[physics.rs](csta/examples/physics.rs). An energy-only DOS cannot reconstruct
magnetization without conditional or joint information. Declare support from
model constraints or exact enumeration; an unvisited cell is not evidence that
it is inaccessible. Large joint supports can be expensive, and the current joint
API has no dedicated checkpoint session.

## 8. Oscillation and energy conservation of a Lennard–Jones pair

**Question.** How does a slightly displaced bound pair vibrate around its
potential minimum, and how sensitive is the trajectory to the time step?

Use two particles with zero total momentum in a large `PeriodicBox<3>` and the
`dynamics::lennard_jones` force callback. The pair potential is
`4 epsilon [(sigma/r)^12 - (sigma/r)^6]`, whose minimum is at
`r = 2^(1/6) sigma`. Start slightly away from that separation and integrate with
`Dynamics` at several time steps.

**Deliverables.** Separation and kinetic/potential energy versus time, total
momentum, and absolute energy drift. Check forces against numerical derivatives
and keep the pair away from half-box image changes.

**Related implementation:** [physics.rs](csta/examples/physics.rs).
This is a two-body NVE problem. The helper evaluates all nearest-image pairs with
no cutoff or tail correction, and cannot be used to claim bulk-fluid accuracy.
Very close separations require smaller steps or produce explicit force errors.

## 9. Equilibrium barrier crossing with temperature replica exchange

**Question.** Does a temperature ladder improve equilibrium exploration of a
double-well coordinate compared with a cold local chain?

Implement a scalar `State` with `U(x) = a (x^2-b^2)^2`, positive `a,b`, and a
symmetric bounded displacement proposal. Store the old coordinate for exact
rejection rollback. Compare a single `Metropolis` chain with
`replica::ReplicaExchange` over a finite positive, increasing beta ladder.
Measure well occupancy, sign changes, energy overlap, and cold–hot–cold round trips.

**Deliverables.** Symmetric equilibrium well populations, histograms at each
beta, and comparison with a numerical quadrature of `exp(-beta U(x))` over a
validated sufficiently large interval. Collect measurements between explicit
exchange rounds; the temperature-exchange API exposes replicas for inspection.

**Related implementation:** the temperature ladder in
[physics.rs](csta/examples/physics.rs), using an Ising model instead of a double well.
CSTA does not provide a quadrature engine or automatically tune the ladder.
Exchange trajectories describe equilibrium sampling, not physical transition
rates. Normalizing the unbounded double well requires positive beta.

## 10. Polarization of a classical freely rotating dipole

**Question.** How does a three-dimensional unit dipole align with an external
field, and how does the classical response differ from a two-state spin?

Use `isotropic_direction` to construct unit vectors and propose independent
uniform orientations on the sphere. Implement `State` with energy `H = -h s_z`
and retain the old vector in each change. Such orientation proposals are symmetric
with respect to solid angle. Sample with `Metropolis` and average `s_z`.

**Deliverables.** Compare polarization with the Langevin function
`L(x) = coth(x) - 1/x`, where `x = beta h`; use its limit `L(0) = 0` and series
`L(x) = x/3 + O(x^3)` near zero to avoid numerical cancellation. Check isotropic
second moments `<s_x^2> = <s_y^2> = <s_z^2> = 1/3` at zero field.

**Limitations.** Component-uniform `Randomizable` vectors are not uniform sphere
orientations. Normalizing a cube sample does not fix its angular distribution.
This model is classical and noninteracting; it has neither quantum level
splitting nor collective magnetic order.

## 11. Finding the useful temperature range of histogram reweighting

**Question.** How far can a canonical trajectory be reweighted before its
estimated observables cease to agree with exact finite-system values?

Sample `wl::models::Ising<4>` canonically at one positive beta, retaining `(E, E)`
for energy or `(E, M)` for magnetization. Call `wl::analysis::reweight` at a grid
of target betas, and compare energy with `Ising::<4>::exact_dos` thermodynamics.

**Deliverables.** Reweighted observables, weight effective sample size, overlap
flags, and discrepancies from exact values. Repeat with independently generated
source trajectories to expose sampling variation.

**Related implementation:** [reweight.rs](csta/examples/reweight.rs).
Weight ESS measures weight concentration, not autocorrelation or unvisited
states. `Overlap::Adequate` is a threshold result, not proof of physical overlap.
This helper performs single-histogram reweighting; it does not implement WHAM or
MBAR. Keep the Hamiltonian fixed unless all required parameter-dependent weights
are constructed explicitly.

## 12. Splitting a density-of-states calculation into energy windows

**Question.** Does overlapping-window sampling reproduce a small Ising chain's
DOS and thermal observables, and where do the windows disagree?

Use `wl::models::Ising<8>` with the exact global discrete energy grid. Partition
its five levels into bin-index windows `0..4` and `2..5`, initialize each walker
inside its support, and use `wl::run_parallel`. Compare the normalized merged DOS
with exact enumeration and with a serial calculation under stated budgets.

**Deliverables.** DOS error by energy, window overlap disagreement, exchange
attempts/accepts, frozen-DOS mixing counts, and thermal energy across beta.
Report wall time and total proposals separately; more workers do not by themselves
imply a fair speedup comparison.

**Related implementation:** [parallel_ising.rs](csta/examples/parallel_ising.rs).
Windows must overlap in at least two bins and collectively cover the grid.
Partial runs provide no successful merged DOS. Interacting walkers are not
independent experiments for error bars; repeat complete parallel experiments
with different seeds. Rayon provides shared-memory execution, not MPI or GPU
sampling.

## Validation and further reading

The five new examples use existing CSTA features without modifying the library
API. Their custom model and force law have focused tests:

```sh
cargo test -p csta --example paramagnet --example harmonic_chain
cargo check -p csta --examples --all-features
cargo clippy -p csta --examples --all-features -- -D warnings
cargo fmt --all -- --check
```

Run each executable as well as these checks. Compare the numeric reference
columns, vary seeds and budgets, and inspect uncertainty assumptions before
using a result in a physics analysis. The defaults are demonstration experiments.

The [library design document](docs/design.pdf) details API contracts, numerical
ranges, and limitations; [TEST_PLAN.md](TEST_PLAN.md) maps broader library tests.
For background on deriving observables from state sums, see David Tong's
[statistical-mechanics fundamentals](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S1.html).
