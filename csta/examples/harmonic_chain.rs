//! A periodic harmonic chain: exact normal-mode motion and Verlet convergence.
//! Run `cargo run -p csta --release --example harmonic_chain`.
use csta::dynamics::{self, Dynamics, Forces};

const SITES: usize = 16;
const SPRING: f64 = 1.0;
const MASS: f64 = 1.0;
const AMPLITUDE: f64 = 0.1;
const MODE: usize = 1;

// Coordinates are displacements from equilibrium, not wrapped particle positions.
// Each spring i -> (i+1)%N is counted once: U = k/2 sum (u[i+1]-u[i])^2.
fn chain_force(u: &[[f64; 1]]) -> dynamics::Result<Forces<1>> {
    let n = u.len();
    if n < 3 || u.iter().any(|x| !x[0].is_finite()) {
        return Err("chain needs at least three finite displacements");
    }
    let mut energy = 0.0;
    let mut forces = vec![[0.0]; n];
    for i in 0..n {
        let j = (i + 1) % n;
        let extension = u[j][0] - u[i][0];
        energy += 0.5 * SPRING * extension * extension;
        forces[i][0] += SPRING * extension;
        forces[j][0] -= SPRING * extension;
    }
    if !energy.is_finite() || forces.iter().any(|f| !f[0].is_finite()) {
        return Err("chain energy or forces overflow");
    }
    Ok((energy, forces))
}

fn simulate(steps: usize) -> dynamics::Result<(f64, f64, f64, f64)> {
    if steps == 0 {
        return Err("need at least one integration step");
    }
    let q = std::f64::consts::TAU * MODE as f64 / SITES as f64;
    let omega = 2.0 * (SPRING / MASS).sqrt() * (0.5 * q).sin().abs();
    let duration = 1.3 * std::f64::consts::TAU / omega;
    let dt = duration / steps as f64;
    let initial: Vec<_> = (0..SITES)
        .map(|i| [AMPLITUDE * (q * i as f64).cos()])
        .collect();
    let mut chain = Dynamics::new(initial.clone(), vec![[0.0]; SITES], vec![MASS; SITES])?;
    let drift = chain.run(steps, dt, chain_force)?;
    let error = chain
        .positions()
        .iter()
        .zip(initial)
        .map(|(actual, start)| (actual[0] - start[0] * (omega * duration).cos()).abs())
        .fold(0.0, f64::max);
    Ok((dt, error, drift.maximum_absolute, chain.momentum()?[0]))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "steps,dt,max_displacement_error,max_absolute_energy_drift,total_momentum,error_ratio"
    );
    let mut previous: Option<f64> = None;
    for steps in [500, 1000, 2000, 4000] {
        let (dt, error, drift, momentum) = simulate(steps)?;
        let ratio = previous.map_or_else(|| "NA".into(), |e| format!("{:.6}", e / error));
        println!("{steps},{dt:.8},{error:.10},{drift:.10},{momentum:.10},{ratio}");
        previous = Some(error);
    }
    eprintln!(
        "N={SITES}, k={SPRING}, mass={MASS}, mode={MODE}. All rows end at 1.3 exact mode periods; halving dt should approach a factor-four error reduction. This chain is harmonic and does not model thermalization."
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn forces_match_energy_gradients_and_translation_symmetry() {
        let u = [[0.1], [-0.2], [0.3], [0.0]];
        let (_, force) = chain_force(&u).unwrap();
        for i in 0..u.len() {
            let mut plus = u;
            let mut minus = u;
            plus[i][0] += 1e-6;
            minus[i][0] -= 1e-6;
            let gradient = (chain_force(&plus).unwrap().0 - chain_force(&minus).unwrap().0) / 2e-6;
            assert!((gradient + force[i][0]).abs() < 1e-10);
        }
        assert!(force.iter().map(|f| f[0]).sum::<f64>().abs() < 1e-14);
        assert_eq!(chain_force(&[[2.0]; 4]).unwrap(), (0.0, vec![[0.0]; 4]));
        assert!(chain_force(&[[0.0]; 2]).is_err());
        assert!(chain_force(&[[f64::NAN]; 3]).is_err());
    }
    #[test]
    fn normal_mode_has_second_order_error() {
        let coarse = simulate(500).unwrap();
        let fine = simulate(1000).unwrap();
        assert!((coarse.1 / fine.1 - 4.0).abs() < 0.01);
        assert!(fine.1 < 1e-5);
        assert!(fine.3.abs() < 1e-12);
        assert!(simulate(0).is_err());
    }
}
