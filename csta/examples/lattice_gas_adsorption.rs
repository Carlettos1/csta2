//! Adsorption on independent sites: density and particle-number fluctuations.
//! Run `cargo run -p csta --release --example lattice_gas_adsorption`.
use csta::{
    Metropolis, Schedule,
    models::{Boundary, LatticeGas},
    statistics::Blocking,
};
use rand::{SeedableRng, rngs::StdRng};

const SIDE: usize = 6;
const SITES: usize = SIDE * SIDE;
const BETA: f64 = 1.0;
const BURN_IN: usize = 20_000;
const SAMPLES: usize = 256_000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "mu,density,density_se,density_exact,number_variance,number_variance_exact,d_density_d_mu,d_density_d_mu_exact"
    );
    for (index, mu) in [-4.0_f64, -2.0, 0.0, 2.0, 4.0].into_iter().enumerate() {
        // J=0 makes sites independent in equilibrium; successive moves remain correlated.
        let model = LatticeGas::new(SIDE, SIDE, Boundary::Periodic, 0.0, mu, vec![false; SITES])?;
        let mut sampler = Metropolis::with_all(
            model,
            (),
            BETA,
            BURN_IN + SAMPLES,
            StdRng::seed_from_u64(72 + index as u64),
        );
        let mut particles = Blocking::new(1024)?;
        sampler.try_run_observed(
            Schedule {
                burn_in: BURN_IN,
                stride: 1,
            },
            None,
            |state, _| particles.push(state.particles() as f64),
        )?;
        let density = particles.moments().mean().ok_or("missing observations")? / SITES as f64;
        let variance = particles
            .moments()
            .population_variance()
            .ok_or("missing variance")?;
        let se = particles.estimate().map_or_else(
            || "NA".into(),
            |e| format!("{:.8}", e.standard_error / SITES as f64),
        );
        let exact_density = 1.0 / (1.0 + (-BETA * mu).exp());
        let exact_variance = SITES as f64 * exact_density * (1.0 - exact_density);
        println!(
            "{mu},{density:.8},{se},{exact_density:.8},{variance:.8},{exact_variance:.8},{:.8},{:.8}",
            BETA * variance / SITES as f64,
            BETA * exact_density * (1.0 - exact_density)
        );
    }
    eprintln!(
        "k_B=1, beta={BETA}, J=0, sites={SITES}. State::energy is H-mu*N; physical H=0. Density SE uses 1024-move blocks; response estimates have no error bars here."
    );
    Ok(())
}
