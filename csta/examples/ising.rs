//! Canonical 2D Ising: J=1, h=0, k_b=1, periodic 4x4 lattice.
use csta::{
    Metropolis, Schedule,
    models::Ising2D,
    statistics::{Blocking, Thermodynamics},
};
use rand::{SeedableRng, rngs::StdRng};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut sampler = Metropolis::with_all(
        Ising2D::aligned(4)?,
        (),
        0.4,
        30_000,
        StdRng::seed_from_u64(42),
    );
    let mut energy = Blocking::new(128)?;
    let mut thermal = Thermodynamics::default();
    let report = sampler.try_run_observed(
        Schedule {
            burn_in: 2000,
            stride: 1,
        },
        None,
        |s, _| {
            let e = s.recompute_energy();
            energy.push(e / 16.0)?;
            thermal.push(e, s.magnetization())
        },
    )?;
    println!(
        "status={:?}, attempts={}, acceptance={:.3}",
        report.status,
        report.attempted,
        sampler.accepted_rate()
    );
    println!(
        "per-site thermodynamics: {:?}",
        thermal.summary(0.4, 1.0, 16)?
    );
    println!(
        "energy/site and block standard error: {:?}",
        energy.estimate()
    );
    println!(
        "Finite-size estimate; inspect mixing and increase block size before interpreting error bars."
    );
    Ok(())
}
