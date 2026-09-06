//! Fixed-step velocity-Verlet for NVE trajectories, independent of Monte Carlo.
//! Force callbacks must be pure functions of positions. Errors leave particles
//! unchanged. Finite but unstable trajectories require a user-chosen drift bound.
use crate::geometry::PeriodicBox;
pub type Result<T> = std::result::Result<T, &'static str>;
pub type Forces<const D: usize> = (f64, Vec<[f64; D]>);
#[derive(Clone, Debug)]
pub struct Dynamics<const D: usize> {
    positions: Vec<[f64; D]>,
    velocities: Vec<[f64; D]>,
    masses: Vec<f64>,
}
#[derive(Clone, Copy, Debug)]
pub struct EnergyDrift {
    pub initial: f64,
    pub final_energy: f64,
    pub maximum_absolute: f64,
}
impl<const D: usize> Dynamics<D> {
    pub fn new(
        positions: Vec<[f64; D]>,
        velocities: Vec<[f64; D]>,
        masses: Vec<f64>,
    ) -> Result<Self> {
        if D == 0
            || positions.is_empty()
            || positions.len() != velocities.len()
            || positions.len() != masses.len()
            || positions
                .iter()
                .chain(&velocities)
                .flatten()
                .any(|x| !x.is_finite())
            || masses.iter().any(|m| !m.is_finite() || *m <= 0.0)
        {
            return Err("invalid particles or masses");
        }
        Ok(Self {
            positions,
            velocities,
            masses,
        })
    }
    pub fn positions(&self) -> &[[f64; D]] {
        &self.positions
    }
    pub fn velocities(&self) -> &[[f64; D]] {
        &self.velocities
    }
    pub fn momentum(&self) -> Result<[f64; D]> {
        let p = std::array::from_fn(|d| {
            self.velocities
                .iter()
                .zip(&self.masses)
                .map(|(v, m)| v[d] * m)
                .sum::<f64>()
        });
        if p.iter().any(|x| !x.is_finite()) {
            return Err("momentum overflow");
        }
        Ok(p)
    }
    pub fn kinetic_energy(&self) -> Result<f64> {
        let k = self
            .velocities
            .iter()
            .zip(&self.masses)
            .map(|(v, m)| 0.5 * m * v.iter().map(|x| x * x).sum::<f64>())
            .sum::<f64>();
        if !k.is_finite() {
            return Err("kinetic energy overflow");
        }
        Ok(k)
    }
    fn validate_forces(&self, f: &Forces<D>) -> Result<()> {
        if !f.0.is_finite()
            || f.1.len() != self.positions.len()
            || f.1.iter().flatten().any(|x| !x.is_finite())
        {
            return Err("invalid potential or forces");
        }
        Ok(())
    }
    pub fn step(
        &mut self,
        dt: f64,
        force: &impl Fn(&[[f64; D]]) -> Result<Forces<D>>,
    ) -> Result<f64> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err("time step must be finite and positive");
        }
        let old = force(&self.positions)?;
        self.validate_forces(&old)?;
        let mut next = self.clone();
        for i in 0..next.positions.len() {
            for d in 0..D {
                next.velocities[i][d] += 0.5 * dt * (old.1[i][d] / self.masses[i]);
                next.positions[i][d] += dt * next.velocities[i][d];
            }
        }
        if next
            .positions
            .iter()
            .chain(&next.velocities)
            .flatten()
            .any(|x| !x.is_finite())
        {
            return Err("unstable nonfinite trajectory");
        }
        let new = force(&next.positions)?;
        next.validate_forces(&new)?;
        for i in 0..next.positions.len() {
            for d in 0..D {
                next.velocities[i][d] += 0.5 * dt * (new.1[i][d] / self.masses[i]);
            }
        }
        let energy = next.kinetic_energy()? + new.0;
        if !energy.is_finite() {
            return Err("total energy overflow");
        }
        *self = next;
        Ok(energy)
    }
    /// Absolute energy drift avoids dividing by zero for zero reference energy.
    pub fn run(
        &mut self,
        steps: usize,
        dt: f64,
        force: impl Fn(&[[f64; D]]) -> Result<Forces<D>>,
    ) -> Result<EnergyDrift> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err("invalid time step");
        }
        let f = force(&self.positions)?;
        self.validate_forces(&f)?;
        let initial = self.kinetic_energy()? + f.0;
        if !initial.is_finite() {
            return Err("total energy overflow");
        }
        let mut report = EnergyDrift {
            initial,
            final_energy: initial,
            maximum_absolute: 0.0,
        };
        for _ in 0..steps {
            report.final_energy = self.step(dt, &force)?;
            report.maximum_absolute = report
                .maximum_absolute
                .max((report.final_energy - initial).abs());
        }
        Ok(report)
    }
}
pub fn harmonic<const D: usize>(positions: &[[f64; D]], k: f64) -> Result<Forces<D>> {
    if !k.is_finite() || k < 0.0 || positions.iter().flatten().any(|x| !x.is_finite()) {
        return Err("invalid harmonic parameters");
    }
    let energy = positions
        .iter()
        .flatten()
        .map(|x| 0.5 * k * x * x)
        .sum::<f64>();
    let forces = positions
        .iter()
        .map(|r| std::array::from_fn(|d| -k * r[d]))
        .collect::<Vec<_>>();
    if !energy.is_finite() || forces.iter().flatten().any(|x| !x.is_finite()) {
        return Err("harmonic overflow");
    }
    Ok((energy, forces))
}
/// All nearest-image pairs, no cutoff or tail correction. Intended for tiny
/// examples; O(N²), and not a production bulk-fluid interaction model.
pub fn lennard_jones<const D: usize>(
    positions: &[[f64; D]],
    cell: &PeriodicBox<D>,
    epsilon: f64,
    sigma: f64,
) -> Result<Forces<D>> {
    if positions.iter().flatten().any(|x| !x.is_finite())
        || !epsilon.is_finite()
        || epsilon <= 0.0
        || !sigma.is_finite()
        || sigma <= 0.0
    {
        return Err("invalid Lennard-Jones parameters");
    }
    let mut forces = vec![[0.0; D]; positions.len()];
    let mut energy = 0.0;
    for i in 0..positions.len() {
        for j in i + 1..positions.len() {
            let r = cell.displacement(positions[i], positions[j])?;
            let r2 = r.iter().map(|x| x * x).sum::<f64>();
            if r2 == 0.0 || !r2.is_finite() {
                return Err("coincident or unrepresentable pair separation");
            }
            let x = (sigma * sigma / r2).powi(3);
            energy += 4.0 * epsilon * x * (x - 1.0);
            let factor = 24.0 * epsilon * x * (1.0 - 2.0 * x) / r2;
            for (d, rd) in r.iter().enumerate() {
                forces[i][d] += factor * rd;
                forces[j][d] -= factor * rd;
            }
        }
    }
    if !energy.is_finite() || forces.iter().flatten().any(|x| !x.is_finite()) {
        return Err("pair force overflow");
    }
    Ok((energy, forces))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn harmonic_accuracy_and_errors() {
        let mut errors = Vec::new();
        for n in [100, 200] {
            let mut d = Dynamics::new(vec![[1.0]], vec![[0.0]], vec![1.0]).unwrap();
            let report = d.run(n, 1.0 / n as f64, |x| harmonic(x, 1.0)).unwrap();
            errors.push((d.positions()[0][0] - 1.0_f64.cos()).abs());
            assert!(report.maximum_absolute < 2e-5);
        }
        assert!((errors[0] / errors[1] - 4.0).abs() < 0.01);
        let mut d = Dynamics::new(vec![[1.0]], vec![[0.0]], vec![1.0]).unwrap();
        let before = d.positions.clone();
        assert!(d.step(f64::MAX, &|x| harmonic(x, 1.0)).is_err());
        assert_eq!(d.positions, before);
        for dt in [0.0, -1.0, f64::NAN] {
            assert!(d.run(0, dt, |x| harmonic(x, 1.0)).is_err());
        }
        assert_eq!(
            d.run(0, 0.1, |x| harmonic(x, 1.0))
                .unwrap()
                .maximum_absolute,
            0.0
        );
        assert!(Dynamics::new(vec![[0.0]], vec![[0.0]], vec![0.0]).is_err());
    }
    #[test]
    fn pair_forces_and_momentum() {
        let cell = PeriodicBox::new([10.0; 3]).unwrap();
        let x = vec![[0.0; 3], [1.3, 0.2, 0.0]];
        let (e, f) = lennard_jones(&x, &cell, 1.0, 1.0).unwrap();
        assert!(e < 0.0);
        for i in 0..2 {
            for d in 0..3 {
                let mut plus = x.clone();
                let mut minus = x.clone();
                plus[i][d] += 1e-6;
                minus[i][d] -= 1e-6;
                let derivative = (lennard_jones(&plus, &cell, 1.0, 1.0).unwrap().0
                    - lennard_jones(&minus, &cell, 1.0, 1.0).unwrap().0)
                    / 2e-6;
                assert!((derivative + f[i][d]).abs() < 1e-7);
            }
        }
        let mut sim = Dynamics::new(x, vec![[0.0; 3]; 2], vec![1.0; 2]).unwrap();
        sim.run(100, 0.001, |x| lennard_jones(x, &cell, 1.0, 1.0))
            .unwrap();
        assert!(sim.momentum().unwrap().iter().all(|p| p.abs() < 1e-14));
        assert!(lennard_jones(&[[0.0; 3]; 2], &cell, 1.0, 1.0).is_err());
    }
}
