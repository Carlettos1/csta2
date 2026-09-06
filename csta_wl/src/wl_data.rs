use crate::{EnergyGrid, Error, Result};

/// DOS values are ln(bin mass), including any degeneracy. Interval quadrature
/// uses bin centers. Normalize or anchor before interpreting absolute F or S.
#[derive(Clone, Debug, PartialEq)]
pub struct WLData {
    grid: EnergyGrid,
    dos: Vec<f64>,
    bins: Vec<u64>,
}
impl WLData {
    pub fn new(grid: EnergyGrid, dos: Vec<f64>, bins: Vec<u64>) -> Result<Self> {
        if dos.len() != grid.len() || bins.len() != grid.len() {
            return Err(Error::Invalid("DOS, histogram and grid lengths differ"));
        }
        if dos.iter().any(|v| v.is_nan() || *v == f64::INFINITY)
            || !dos.iter().any(|v| v.is_finite())
        {
            return Err(Error::Invalid(
                "DOS needs finite mass; only -infinity is allowed for absent levels",
            ));
        }
        Ok(Self { grid, dos, bins })
    }
    pub fn grid(&self) -> &EnergyGrid {
        &self.grid
    }
    pub fn dos(&self) -> &[f64] {
        &self.dos
    }
    pub fn bins(&self) -> &[u64] {
        &self.bins
    }

    /// Set sum(g) to exp(log_count), using logarithmic input to avoid overflow.
    pub fn normalize_log_count(&mut self, log_count: f64) -> Result<()> {
        if !log_count.is_finite() {
            return Err(Error::Invalid("log state count must be finite"));
        }
        let max = self.dos.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let relative: Vec<_> = self.dos.iter().map(|v| *v - max).collect();
        let z = log_sum_exp(&relative)?;
        let next: Vec<_> = relative.iter().map(|v| (*v - z) + log_count).collect();
        self.replace_dos(next)
    }
    pub fn anchor(&mut self, bin: usize, log_degeneracy: f64) -> Result<()> {
        if bin >= self.dos.len() || !self.dos[bin].is_finite() || !log_degeneracy.is_finite() {
            return Err(Error::Invalid(
                "anchor requires a finite accessible bin and log degeneracy",
            ));
        }
        let base = self.dos[bin];
        self.replace_dos(
            self.dos
                .iter()
                .map(|v| (*v - base) + log_degeneracy)
                .collect(),
        )
    }
    fn replace_dos(&mut self, next: Vec<f64>) -> Result<()> {
        if next
            .iter()
            .zip(&self.dos)
            .any(|(n, old)| old.is_finite() && !n.is_finite())
        {
            return Err(Error::Numerical("normalization overflow"));
        }
        self.dos = next;
        Ok(())
    }
    /// beta = 1/(k_b T). Any finite beta is valid for a finite supported spectrum.
    fn shifted_weights(&self, beta: f64) -> Result<(Vec<f64>, f64, f64)> {
        if !beta.is_finite() {
            return Err(Error::Invalid("beta must be finite"));
        }
        let log_base = self.dos.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut indices = (0..self.dos.len()).filter(|i| self.dos[*i].is_finite());
        let index = if beta >= 0.0 {
            indices.next()
        } else {
            indices.next_back()
        }
        .unwrap();
        let energy_base = self.grid.energies()[index];
        let weights: Vec<_> = self
            .dos
            .iter()
            .zip(self.grid.energies())
            .map(|(g, e)| {
                if *g == f64::NEG_INFINITY {
                    f64::NEG_INFINITY
                } else if beta == 0.0 {
                    *g - log_base
                } else {
                    (*g - log_base) - beta * (*e - energy_base)
                }
            })
            .collect();
        if weights.iter().any(|v| v.is_nan() || *v == f64::INFINITY)
            || !weights.iter().any(|v| v.is_finite())
        {
            return Err(Error::Numerical("unrepresentable Boltzmann weights"));
        }
        Ok((weights, log_base, energy_base))
    }
    pub fn log_partition_function(&self, beta: f64) -> Result<f64> {
        let (weights, g0, e0) = self.shifted_weights(beta)?;
        finite(
            log_sum_exp(&weights)? + g0 - beta * e0,
            "log partition function overflow",
        )
    }
    pub fn energy_distribution(&self, beta: f64) -> Result<Vec<f64>> {
        let (weights, _, _) = self.shifted_weights(beta)?;
        let z = log_sum_exp(&weights)?;
        let mut probabilities: Vec<_> = weights.iter().map(|v| (*v - z).exp()).collect();
        let sum: f64 = probabilities.iter().sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(Error::Numerical("probability normalization failed"));
        }
        probabilities.iter_mut().for_each(|v| *v /= sum);
        Ok(probabilities)
    }
    pub fn energy_moments(&self, beta: f64) -> Result<(f64, f64)> {
        let p = self.energy_distribution(beta)?;
        let mut first = 0.0;
        let mut second = 0.0;
        for (p, e) in p
            .iter()
            .zip(self.grid.energies())
            .filter(|(p, _)| **p > 0.0)
        {
            first += p * e;
            second += (p * e) * e;
        }
        Ok((
            finite(first, "mean overflow")?,
            finite(second, "second moment overflow")?,
        ))
    }
    /// Centered variance avoids catastrophic cancellation in <E²> - <E>².
    pub fn energy_variance(&self, beta: f64) -> Result<f64> {
        let p = self.energy_distribution(beta)?;
        let mean: f64 = p.iter().zip(self.grid.energies()).map(|(p, e)| p * e).sum();
        let variance: f64 = p
            .iter()
            .zip(self.grid.energies())
            .filter(|(p, _)| **p > 0.0)
            .map(|(p, e)| {
                let d = e - mean;
                (p * d) * d
            })
            .sum();
        finite(variance, "variance overflow")
    }
    /// Total heat capacity C = k_b * beta² * Var(E).
    pub fn specific_heat(&self, beta: f64, k_b: f64) -> Result<f64> {
        positive(k_b, "k_b must be finite and positive")?;
        let var = self.energy_variance(beta)?;
        if var == 0.0 {
            return Ok(0.0);
        }
        finite(
            ((beta * var.sqrt()) * k_b.sqrt()).powi(2),
            "heat capacity overflow",
        )
    }
    /// F = -ln(Z)/beta. k_b is validated; it is already included in beta.
    pub fn free_energy(&self, log_z: f64, beta: f64, k_b: f64) -> Result<f64> {
        positive(k_b, "k_b must be finite and positive")?;
        if !beta.is_finite() || beta == 0.0 || !log_z.is_finite() {
            return Err(Error::Invalid(
                "free energy requires nonzero finite beta and finite log Z",
            ));
        }
        finite(-log_z / beta, "free energy overflow")
    }
    pub fn free_energy2(&self, beta: f64, k_b: f64) -> Result<f64> {
        self.free_energy(self.log_partition_function(beta)?, beta, k_b)
    }
    /// Negative finite temperature is supported for bounded spectra; zero is not.
    pub fn entropy(&self, e_avg: f64, free_energy: f64, temperature: f64) -> Result<f64> {
        if !e_avg.is_finite()
            || !free_energy.is_finite()
            || !temperature.is_finite()
            || temperature == 0.0
        {
            return Err(Error::Invalid(
                "entropy needs finite E, F and nonzero finite temperature",
            ));
        }
        finite((e_avg - free_energy) / temperature, "entropy overflow")
    }
    pub fn microcanonical_entropy(&self, k_b: f64) -> Result<Vec<f64>> {
        positive(k_b, "k_b must be finite and positive")?;
        self.dos
            .iter()
            .map(|g| {
                if *g == f64::NEG_INFINITY {
                    Ok(f64::NEG_INFINITY)
                } else {
                    finite(k_b * g, "entropy overflow")
                }
            })
            .collect()
    }
    /// Endpoints and stencils crossing absent levels return None. At least three
    /// levels are required. A zero slope gives +infinity; negative slopes are valid.
    pub fn microcanonical_temperature(&self, k_b: f64) -> Result<Vec<Option<f64>>> {
        let s = self.microcanonical_entropy(k_b)?;
        if s.len() < 3 {
            return Err(Error::Invalid(
                "temperature derivative needs at least three levels",
            ));
        }
        let mut temp = vec![None; s.len()];
        for i in 1..s.len() - 1 {
            if s[i - 1..=i + 1].iter().all(|v| v.is_finite()) {
                let slope = derivative(self.grid.energies(), &s, i)?;
                temp[i] = Some(if slope == 0.0 {
                    f64::INFINITY
                } else {
                    1.0 / slope
                });
            }
        }
        Ok(temp)
    }
}
fn positive(x: f64, message: &'static str) -> Result<()> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(Error::Invalid(message))
    }
}
fn finite(x: f64, message: &'static str) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::Numerical(message))
    }
}
/// Nonuniform three-point derivative; one-sided secants at endpoints for joining.
pub(crate) fn derivative(e: &[f64], v: &[f64], i: usize) -> Result<f64> {
    if e.len() < 2 || e.len() != v.len() || i >= e.len() {
        return Err(Error::InsufficientOverlap);
    }
    let slope = if i == 0 {
        (v[1] - v[0]) / (e[1] - e[0])
    } else if i == e.len() - 1 {
        (v[i] - v[i - 1]) / (e[i] - e[i - 1])
    } else {
        let left = e[i] - e[i - 1];
        let right = e[i + 1] - e[i];
        ((v[i] - v[i - 1]) / left) * (right / (left + right))
            + ((v[i + 1] - v[i]) / right) * (left / (left + right))
    };
    finite(slope, "nonfinite entropy slope")
}
/// Empty/all -infinity inputs return -infinity. +infinity dominates finite values.
/// NaN is always an error. Callers requiring positive finite mass validate it.
pub fn log_sum_exp(values: &[f64]) -> Result<f64> {
    if values.iter().any(|v| v.is_nan()) {
        return Err(Error::Numerical("NaN in log-sum-exp"));
    }
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return Ok(max);
    }
    finite(
        max + values.iter().map(|v| (*v - max).exp()).sum::<f64>().ln(),
        "log-sum-exp overflow",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn data(e: Vec<f64>, g: Vec<f64>) -> WLData {
        let n = e.len();
        WLData::new(EnergyGrid::discrete(e).unwrap(), g, vec![0; n]).unwrap()
    }
    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-10 * (1.0 + b.abs()), "{a} != {b}");
    }
    #[test]
    fn log_sum_exp_extremes() {
        close(log_sum_exp(&[0.0, 0.0]).unwrap(), 2_f64.ln());
        close(log_sum_exp(&[1000.0, 1000.0]).unwrap(), 1000.0 + 2_f64.ln());
        close(log_sum_exp(&[-1000.0]).unwrap(), -1000.0);
        assert_eq!(log_sum_exp(&[]).unwrap(), f64::NEG_INFINITY);
        assert_eq!(
            log_sum_exp(&[f64::NEG_INFINITY; 2]).unwrap(),
            f64::NEG_INFINITY
        );
        assert_eq!(log_sum_exp(&[f64::INFINITY, 0.0]).unwrap(), f64::INFINITY);
        assert!(log_sum_exp(&[f64::NAN, f64::INFINITY]).is_err());
        close(log_sum_exp(&[f64::NEG_INFINITY, 4.0]).unwrap(), 4.0);
    }
    #[test]
    fn one_and_two_level_exact_thermodynamics() {
        let one = data(vec![2.0], vec![3_f64.ln()]);
        assert_eq!(one.energy_moments(2.0).unwrap(), (2.0, 4.0));
        assert_eq!(one.specific_heat(1e300, 2.0).unwrap(), 0.0);
        let two = data(vec![0.0, 2.0], vec![0.0, 3_f64.ln()]);
        for beta in [-2.0, 0.0, 0.5, 3.0] {
            let w = 3.0 * (-2.0_f64 * beta).exp();
            let p = w / (1.0 + w);
            let (e, e2) = two.energy_moments(beta).unwrap();
            close(e, 2.0 * p);
            close(e2, 4.0 * p);
            close(
                two.specific_heat(beta, 2.0).unwrap(),
                2.0 * beta * beta * 4.0 * p * (1.0 - p),
            );
            if beta != 0.0 {
                close(two.free_energy2(beta, 2.0).unwrap(), -(1.0 + w).ln() / beta);
            }
        }
        assert!(two.free_energy2(0.0, 1.0).is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(two.specific_heat(1.0, bad).is_err());
        }
        assert!(two.energy_distribution(f64::NAN).is_err());
        assert!(two.energy_distribution(f64::INFINITY).is_err());
        assert!(two.entropy(1.0, 0.0, 0.0).is_err());
        close(two.entropy(1.0, 0.0, -2.0).unwrap(), -0.5);
    }
    #[test]
    fn offsets_normalization_and_extreme_beta() {
        let original = data(vec![-1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]);
        let mut shifted = data(vec![-1.0, 0.0, 1.0], vec![10000.0, 10001.0, 10000.0]);
        for b in [-1e300, -2.0, 0.0, 2.0, 1e300] {
            let p = shifted.energy_distribution(b).unwrap();
            close(p.iter().sum(), 1.0);
            for (a, b) in p.iter().zip(original.energy_distribution(b).unwrap()) {
                close(*a, b);
            }
            close(
                shifted.energy_variance(b).unwrap(),
                original.energy_variance(b).unwrap(),
            );
        }
        close(
            shifted.log_partition_function(1.0).unwrap()
                - original.log_partition_function(1.0).unwrap(),
            10000.0,
        );
        shifted.normalize_log_count(8_f64.ln()).unwrap();
        close(log_sum_exp(shifted.dos()).unwrap(), 8_f64.ln());
        shifted.anchor(0, 2_f64.ln()).unwrap();
        close(shifted.dos()[0], 2_f64.ln());
        assert!(shifted.anchor(5, 0.0).is_err());
    }
    #[test]
    fn shape_and_mass_validation() {
        let g = EnergyGrid::discrete(vec![0.0]).unwrap();
        assert!(WLData::new(g.clone(), vec![0.0], vec![]).is_err());
        for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(WLData::new(g.clone(), vec![v], vec![0]).is_err());
        }
        let d = data(vec![0.0, 1.0], vec![0.0, f64::NEG_INFINITY]);
        assert_eq!(d.energy_distribution(1.0).unwrap(), vec![1.0, 0.0]);
    }
    #[test]
    fn centered_variance_and_microcanonical_edges() {
        let d = data(vec![1e10, 1e10 + 1.0], vec![0.0; 2]);
        close(d.energy_variance(0.0).unwrap(), 0.25);
        assert!(d.microcanonical_temperature(1.0).is_err());
        assert!(
            data(vec![0.0], vec![0.0])
                .microcanonical_temperature(1.0)
                .is_err()
        );
        for slope in [-2.0, 0.0, 2.0] {
            let d = data(vec![0.0, 1.0, 3.0], vec![0.0, slope, 3.0 * slope]);
            let t = d.microcanonical_temperature(2.0).unwrap();
            assert_eq!(t[0], None);
            assert_eq!(t[2], None);
            if slope == 0.0 {
                assert_eq!(t[1], Some(f64::INFINITY));
            } else {
                close(t[1].unwrap(), 1.0 / (2.0 * slope));
            }
        }
        assert_eq!(
            data(vec![0.0, 1.0, 2.0], vec![0.0, f64::NEG_INFINITY, 0.0])
                .microcanonical_temperature(1.0)
                .unwrap(),
            vec![None; 3]
        );
    }
}
