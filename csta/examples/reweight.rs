//! DOS observables and independent-run errors, checked against a four-spin state sum.
use csta::{
    State,
    wl::{
        self,
        analysis::{ConditionalMoments, DosEnsemble, IndependentDos},
        models::Ising,
    },
};
use rand::{SeedableRng, rngs::StdRng};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut experiments = Vec::new();
    for seed in 0..4 {
        let grid = Ising::<4>::grid()?;
        let mut conditional = ConditionalMoments::new(grid.clone());
        let mut session = wl::Session::new(
            Ising::<4>::new([1; 4])?,
            (),
            StdRng::seed_from_u64(seed),
            wl::RawWangLandauData::on_grid(grid)?,
            wl::Config {
                preliminary_stages: 2,
                min_stage_steps: 100,
                visits_per_bin: 20,
                sampling_steps: 20_000,
                max_steps: 100_000,
                ..Default::default()
            },
        )?;
        session.advance_observed(100_000, None, |s, e| {
            conditional.record(e, s.spins().iter().map(|s| *s as f64).sum())
        })?;
        let result = session.finish();
        if !result.is_complete() {
            return Err("partial independent DOS run".into());
        }
        let mut data = result.data.process_data()?;
        data.normalize_log_count(16_f64.ln())?;
        println!(
            "seed={seed}: conditional M moments at beta=0.5: {:?}",
            conditional.evaluate(&data, 0.5)?
        );
        experiments.push(IndependentDos::from_data(seed, data)?);
    }
    let summary = DosEnsemble::new(experiments, 16_f64.ln())?.energy(0.5)?;
    println!("DOS estimate across independent seeds: {summary:?}");
    println!(
        "exact enumerated mean energy: {}",
        Ising::<4>::exact_dos()?.energy_moments(0.5)?.0
    );
    let mut canonical = csta::Metropolis::with_all(
        Ising::<4>::new([1; 4])?,
        (),
        0.5,
        21_000,
        StdRng::seed_from_u64(9),
    );
    let mut samples = Vec::new();
    canonical.run_observed(
        csta::Schedule {
            burn_in: 1000,
            stride: 1,
        },
        None,
        |s, _| {
            let energy = s.energy(&mut ());
            samples.push((energy, energy));
        },
    )?;
    println!(
        "canonical energy reweighted to beta=0.6: {:?}",
        wl::analysis::reweight(&samples, 0.5, 0.6)?
    );
    println!(
        "exact beta=0.6 energy: {}",
        Ising::<4>::exact_dos()?.energy_moments(0.6)?.0
    );
    Ok(())
}
