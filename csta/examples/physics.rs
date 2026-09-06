//! Small optional physics workflows. Run `cargo run -p csta --example physics`.
use csta::{
    Metropolis, State,
    dynamics::{Dynamics, harmonic, lennard_jones},
    geometry::PeriodicBox,
    models::{Boundary, Ising2D, LatticeGas},
    replica::ReplicaExchange,
    wl::{
        self,
        analysis::{ConditionalMoments, DosEnsemble, IndependentDos},
        joint::{JointDos, JointGrid},
    },
};
use rand::{SeedableRng, rngs::StdRng};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = StdRng::seed_from_u64(42);
    println!("isotropic spin: {:?}", csta::isotropic_direction(&mut rng));
    println!("Gaussian velocity: {}", csta::gaussian(&mut rng, 0.0, 1.0)?);
    let mut model = Ising2D::aligned(3)?;
    println!("Wolff flipped {} spins", model.wolff_step(0.3, &mut rng)?);
    let mut replicas = ReplicaExchange::new(vec![model; 3], (), vec![0.2, 0.3, 0.4], 42)?;
    replicas.run(100, 10, None)?;
    println!(
        "temperature swaps: {}/{}",
        replicas.swap_accepts, replicas.swap_attempts
    );
    let gas = LatticeGas::new(2, 2, Boundary::Open, 0.0, 0.3, vec![false; 4])?;
    let mut sampler = Metropolis::with_all(gas, (), 1.0, 1000, StdRng::seed_from_u64(42));
    sampler.run(None)?;
    println!("lattice gas particles: {}", sampler.state.particles());
    let mut oscillator = Dynamics::new(vec![[1.0]], vec![[0.0]], vec![1.0])?;
    println!(
        "NVE harmonic drift: {:?}",
        oscillator.run(1000, 0.01, |x| harmonic(x, 1.0))?
    );
    let cell = PeriodicBox::new([10.0; 3])?;
    let mut pair = Dynamics::new(
        vec![[0.0; 3], [1.2, 0.0, 0.0]],
        vec![[0.0; 3]; 2],
        vec![1.0; 2],
    )?;
    println!(
        "NVE pair drift: {:?}",
        pair.run(100, 0.001, |x| lennard_jones(x, &cell, 1.0, 1.0))?
    );
    // Exact two-spin joint DOS; E0 excludes any magnetic field.
    let joint = JointDos::new(
        JointGrid::new(vec![(0.0, -2.0), (0.0, 0.0), (0.0, 2.0)], 3)?,
        vec![0.0, 2_f64.ln(), 0.0],
    )?;
    println!("field reweighting: {:?}", joint.evaluate(1.0, 0.2)?);
    let dos = wl::models::Ising::<4>::exact_dos()?;
    let mut conditional = ConditionalMoments::new(dos.grid().clone());
    for bits in 0..16 {
        let spins = std::array::from_fn(|i| if bits & (1 << i) == 0 { -1 } else { 1 });
        let s = wl::models::Ising::<4>::new(spins)?;
        conditional.record(s.energy(&mut ()), spins.iter().map(|x| *x as f64).sum())?;
    }
    println!(
        "conditional M moments: {:?}",
        conditional.evaluate(&dos, 0.5)?
    );
    let ensemble = DosEnsemble::new(vec![IndependentDos::from_data(0, dos)?], 16_f64.ln())?;
    println!(
        "one reference estimate (uncertainty unavailable): {:?}",
        ensemble.energy(0.5)?
    );
    Ok(())
}
