//! Quantum-oscillator heat capacity: sampled DOS, exact cutoff, and infinite spectrum.
//! Run `cargo run -p csta --release --example oscillator_dos`.
use csta::wl::{
    self, Config, RawWangLandauData,
    analysis::{DosEnsemble, IndependentDos},
    models::Oscillators,
};
use rand::{SeedableRng, rngs::StdRng};

const OSCILLATORS: usize = 3;
const MAX_QUANTA: u32 = 32;
const RUNS: u64 = 4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let exact = Oscillators::<OSCILLATORS>::exact_dos(MAX_QUANTA)?;
    let log_count = wl::log_sum_exp(exact.dos())?;
    let mut experiments = Vec::new();
    for id in 0..RUNS {
        let result = wl::run(
            Oscillators::<OSCILLATORS>::new([0; OSCILLATORS])?,
            (),
            StdRng::seed_from_u64(200 + id),
            RawWangLandauData::on_grid(exact.grid().clone())?,
            Config {
                preliminary_stages: 2,
                min_stage_steps: 500,
                visits_per_bin: 30,
                sampling_steps: 250_000,
                max_steps: 1_000_000,
                ..Config::default()
            },
            None,
        )?;
        eprintln!("run={id}: {:?}", result.diagnostics);
        if !result.is_complete() || result.data.lifetime_bins().contains(&0) {
            return Err(
                format!("run {id} did not complete with full declared support visited").into(),
            );
        }
        experiments.push(IndependentDos::from_run(id, result)?);
    }
    let ensemble = DosEnsemble::new(experiments, log_count)?;
    println!(
        "beta,energy_per_oscillator,energy_se,energy_exact_cutoff,energy_infinite,heat_capacity_per_oscillator,heat_capacity_se,heat_capacity_exact_cutoff,heat_capacity_infinite"
    );
    for beta in [0.2_f64, 0.5, 1.0, 2.0, 4.0] {
        let energy = ensemble.evaluate(|d| Ok(d.energy_moments(beta)?.0 / OSCILLATORS as f64))?;
        let heat = ensemble.evaluate(|d| Ok(d.specific_heat(beta, 1.0)? / OSCILLATORS as f64))?;
        let infinite_energy = 0.5 + 1.0 / beta.exp_m1();
        let infinite_heat = beta * beta * (-beta).exp() / (-(-beta).exp_m1()).powi(2);
        println!(
            "{beta},{:.8},{:.8},{:.8},{infinite_energy:.8},{:.8},{:.8},{:.8},{infinite_heat:.8}",
            energy.mean,
            energy
                .standard_error
                .ok_or("need independent DOS replicates")?,
            exact.energy_moments(beta)?.0 / OSCILLATORS as f64,
            heat.mean,
            heat.standard_error
                .ok_or("need independent DOS replicates")?,
            exact.specific_heat(beta, 1.0)? / OSCILLATORS as f64
        );
    }
    eprintln!(
        "k_B=hbar*omega=1, oscillators={OSCILLATORS}, TOTAL quanta <= {MAX_QUANTA}. SEs describe {RUNS} independent DOS experiments; the exact-cutoff/infinite difference is truncation bias, not sampling error."
    );
    Ok(())
}
