//! Finite-size recipe: tabular Binder curves and block errors, no critical-point claim.
use csta::{
    State,
    models::Ising2D,
    statistics::{Thermodynamics, correlated_estimate, thermal_errors},
};
use rand::{SeedableRng, rngs::StdRng};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("size,beta,energy_per_site,susceptibility_per_site,binder,binder_se,energy_tau");
    for size in [3, 4, 6] {
        for beta in [0.3, 0.4, 0.5] {
            let mut model = Ising2D::aligned(size)?;
            let mut rng = StdRng::seed_from_u64(42 + size as u64);
            let mut moments = Thermodynamics::default();
            let mut observations = Vec::new();
            for i in 0..9000 {
                model.wolff_step(beta, &mut rng)?;
                if i >= 1000 {
                    let pair = (model.energy(&mut ()), model.magnetization());
                    moments.push(pair.0, pair.1)?;
                    observations.push(pair);
                }
            }
            let s = moments.summary(beta, 1.0, size * size)?;
            let correlation =
                correlated_estimate(&observations.iter().map(|p| p.0).collect::<Vec<_>>())?;
            let errors = thermal_errors(&observations, beta, 1.0, size * size, 128)?;
            println!(
                "{size},{beta},{},{},{:?},{:?},{:?}",
                s.energy,
                s.susceptibility,
                s.binder,
                errors.and_then(|e| e.binder),
                correlation.map(|c| c.tau)
            );
        }
    }
    eprintln!(
        "Compare several block sizes and beta grids; peaks/crossings at these finite sizes are not exact critical points."
    );
    Ok(())
}
