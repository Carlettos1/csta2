use csta_wl::{
    models::Ising, run_parallel, sample_in_support, Config, ParallelConfig, RawWangLandauData,
};
use rand::{rngs::StdRng, SeedableRng};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let grid = Ising::<8>::grid()?;
    let windows = vec![0..4, 2..5];
    let mut initialization_rng = StdRng::seed_from_u64(7);
    let mut states = Vec::new();
    for window in &windows {
        let mask: Vec<_> = (0..grid.len()).map(|i| window.contains(&i)).collect();
        let dos = mask
            .iter()
            .map(|a| if *a { 0.0 } else { f64::NEG_INFINITY })
            .collect();
        let support = RawWangLandauData::with_support(grid.clone(), mask, dos)?;
        states.push(sample_in_support::<Ising<8>>(
            &mut initialization_rng,
            &(),
            &support,
            1000,
        )?);
    }
    let config = ParallelConfig {
        windows,
        walkers_per_window: 1,
        exchange_every: 100,
        seed: 42,
        run: Config {
            preliminary_stages: 3,
            min_stage_steps: 100,
            visits_per_bin: 20,
            sampling_steps: 100_000,
            max_steps: 1_000_000,
            ..Config::default()
        },
    };
    let result = run_parallel(grid, states, (), config, None)?;
    println!(
        "status={:?}, exchanges={}/{}",
        result.stop_reason, result.accepted_exchanges, result.exchange_attempts
    );
    let mut data = result.merged.ok_or("no complete merged DOS")?;
    data.normalize_log_count(256_f64.ln())?;
    let reference = Ising::<8>::exact_dos()?;
    for ((e, g), exact) in data
        .grid()
        .energies()
        .iter()
        .zip(data.dos())
        .zip(reference.dos())
    {
        println!("E={e}, ln g={g:.6}, exact={exact:.6}");
    }
    Ok(())
}
