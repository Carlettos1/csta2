use csta_wl::{
    Config, RawWangLandauData,
    models::{Ising, Oscillators},
    run,
};
use rand::{SeedableRng, rngs::StdRng};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = Config {
        preliminary_stages: 3,
        min_stage_steps: 100,
        visits_per_bin: 20,
        sampling_steps: 100_000,
        max_steps: 1_000_000,
        ..Config::default()
    };
    let (mut data, reference) = if std::env::args().nth(1).as_deref() == Some("qho") {
        let r = run(
            Oscillators::<2>::new([0; 2])?,
            (),
            StdRng::seed_from_u64(42),
            RawWangLandauData::on_grid(Oscillators::<2>::grid(10)?)?,
            config,
            None,
        )?;
        println!("Oscillator diagnostics: {:?}", r.diagnostics);
        if !r.is_complete() {
            return Err("simulation did not complete".into());
        }
        let mut data = r.data.process_data()?;
        data.anchor(0, 0.0)?;
        (data, Oscillators::<2>::exact_dos(10)?)
    } else {
        let r = run(
            Ising::<8>::new([1; 8])?,
            (),
            StdRng::seed_from_u64(42),
            RawWangLandauData::on_grid(Ising::<8>::grid()?)?,
            config,
            None,
        )?;
        println!("Ising diagnostics: {:?}", r.diagnostics);
        if !r.is_complete() {
            return Err("simulation did not complete".into());
        }
        (r.data.process_data()?, Ising::<8>::exact_dos()?)
    };
    data.normalize_log_count(csta_wl::log_sum_exp(reference.dos())?)?;
    println!("energy\tln g (WL)\tln g (exact)");
    for ((e, g), exact) in data
        .grid()
        .energies()
        .iter()
        .zip(data.dos())
        .zip(reference.dos())
    {
        println!("{e}\t{g:.6}\t{exact:.6}");
    }
    for beta in [0.0, 0.5, 1.0] {
        println!(
            "beta={beta}: <E>={:.6}, exact={:.6}, C(k_b=1)={:.6}",
            data.energy_moments(beta)?.0,
            reference.energy_moments(beta)?.0,
            data.specific_heat(beta, 1.0)?
        );
    }
    Ok(())
}
