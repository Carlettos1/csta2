//! Stable streaming moments and conservative uncertainty estimates.
use crate::{Error, Result};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Moments {
    n: u64,
    mean: f64,
    m2: f64,
}
impl Moments {
    pub fn count(&self) -> u64 {
        self.n
    }
    pub fn mean(&self) -> Option<f64> {
        (self.n > 0).then_some(self.mean)
    }
    pub fn variance(&self) -> Option<f64> {
        (self.n > 1).then(|| self.m2 / (self.n - 1) as f64)
    }
    pub fn population_variance(&self) -> Option<f64> {
        (self.n > 0).then(|| self.m2 / self.n as f64)
    }
    pub fn push(&mut self, x: f64) -> Result<()> {
        self.merge(&Self {
            n: 1,
            mean: x,
            m2: 0.0,
        })
    }
    pub fn merge(&mut self, other: &Self) -> Result<()> {
        if other.n == 0 {
            return Ok(());
        }
        if !other.mean.is_finite() || !other.m2.is_finite() {
            return Err(Error("nonfinite observation"));
        }
        if self.n == 0 {
            *self = *other;
            return Ok(());
        }
        let n = self
            .n
            .checked_add(other.n)
            .ok_or(Error("sample count overflow"))?;
        let d = other.mean - self.mean;
        let mean = self.mean + d * (other.n as f64 / n as f64);
        let m2 = self.m2 + other.m2 + d * d * (self.n as f64 * (other.n as f64 / n as f64));
        if !mean.is_finite() || !m2.is_finite() {
            return Err(Error("moment overflow"));
        }
        *self = Self { n, mean, m2 };
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct Covariance {
    x: Moments,
    y: Moments,
    c: f64,
}
impl Covariance {
    pub fn push(&mut self, x: f64, y: f64) -> Result<()> {
        let mut next = *self;
        next.x.push(x)?;
        next.y.push(y)?;
        next.c += (x - self.x.mean) * (y - next.y.mean);
        if !next.c.is_finite() {
            return Err(Error("covariance overflow"));
        }
        *self = next;
        Ok(())
    }
    pub fn covariance(&self) -> Option<f64> {
        (self.x.n > 1).then(|| self.c / (self.x.n - 1) as f64)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Estimate {
    pub mean: f64,
    pub standard_error: f64,
    pub effective_samples: f64,
    pub samples: u64,
}
/// Fixed nonoverlapping blocks, constant memory. Unfinished blocks are excluded
/// from uncertainty, and at least 8 completed blocks are required. Choose a block
/// longer than the correlation scale; this cannot diagnose equilibration.
#[derive(Clone, Debug)]
pub struct Blocking {
    size: u64,
    current: Moments,
    blocks: Moments,
    all: Moments,
    completed: Moments,
}
impl Default for Blocking {
    fn default() -> Self {
        Self::new(64).expect("valid default")
    }
}
impl Blocking {
    pub fn new(size: u64) -> Result<Self> {
        if size == 0 {
            return Err(Error("block size must be positive"));
        }
        Ok(Self {
            size,
            current: Moments::default(),
            blocks: Moments::default(),
            all: Moments::default(),
            completed: Moments::default(),
        })
    }
    pub fn push(&mut self, x: f64) -> Result<()> {
        let mut next = self.clone();
        next.current.push(x)?;
        next.all.push(x)?;
        if next.current.n == next.size {
            next.blocks.push(next.current.mean)?;
            next.completed.merge(&next.current)?;
            next.current = Moments::default();
        }
        *self = next;
        Ok(())
    }
    pub fn moments(&self) -> &Moments {
        &self.all
    }
    pub fn incomplete_samples(&self) -> u64 {
        self.current.n
    }
    pub fn estimate(&self) -> Option<Estimate> {
        if self.blocks.n < 8 {
            return None;
        }
        let variance = self.blocks.variance()?;
        // A constant trajectory gives no evidence that other states were explored.
        if variance <= 0.0 {
            return None;
        }
        let se2 = variance / self.blocks.n as f64;
        if se2 <= 0.0 || !se2.is_finite() {
            return None;
        }
        Some(Estimate {
            mean: self.blocks.mean,
            standard_error: se2.sqrt(),
            effective_samples: (self.completed.variance()? / se2)
                .min((self.blocks.n * self.size) as f64),
            samples: self.blocks.n * self.size,
        })
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Correlation {
    pub estimate: Estimate,
    pub tau: f64,
    pub max_lag: usize,
}
/// Initial-positive-pair autocorrelation estimate (reversible stationary chains).
/// tau = 1 + 2 sum rho, ESS = n/tau; conservatively floor tau at 1, so negative
/// correlations never claim ESS > n. None means too short/constant/unresolved.
/// O(n * min(n/2, 4096)) work; retained samples are explicitly opt-in.
pub fn correlated_estimate(values: &[f64]) -> Result<Option<Correlation>> {
    let mut m = Moments::default();
    for x in values {
        m.push(*x)?;
    }
    let n = values.len();
    if n < 32 || m.m2 == 0.0 {
        return Ok(None);
    }
    let gamma0 = m.m2 / n as f64;
    let cov = |lag: usize| {
        values[..n - lag]
            .iter()
            .zip(&values[lag..])
            .map(|(a, b)| (a - m.mean) * (b - m.mean))
            .sum::<f64>()
            / n as f64
    };
    let mut sum = 0.0_f64;
    let mut previous = f64::INFINITY;
    for k in (0..(n / 2).min(4096) - 1).step_by(2) {
        let pair = (cov(k) + cov(k + 1)) / gamma0;
        if !pair.is_finite() {
            return Err(Error("autocovariance overflow"));
        }
        if pair <= 0.0 {
            let tau = (2.0 * sum - 1.0).max(1.0);
            if (n as f64) < 20.0 * tau {
                return Ok(None);
            }
            return Ok(Some(Correlation {
                estimate: Estimate {
                    mean: m.mean,
                    standard_error: (gamma0 * tau / n as f64).sqrt(),
                    effective_samples: n as f64 / tau,
                    samples: n as u64,
                },
                tau,
                max_lag: k + 1,
            }));
        }
        previous = previous.min(pair);
        sum += previous;
    }
    Ok(None)
}

/// Total E and signed total M are inputs. k_b is explicit; beta=1/(k_b T).
#[derive(Clone, Debug, Default)]
pub struct Thermodynamics {
    energy: Moments,
    magnetization: Moments,
    m2: Moments,
    m4: Moments,
}
#[derive(Clone, Copy, Debug)]
pub struct ThermalSummary {
    pub energy: f64,
    pub magnetization: f64,
    pub heat_capacity: f64,
    pub susceptibility: f64,
    pub binder: Option<f64>,
    pub samples: u64,
}
impl Thermodynamics {
    pub fn push(&mut self, energy: f64, magnetization: f64) -> Result<()> {
        let mut next = self.clone();
        next.energy.push(energy)?;
        next.magnetization.push(magnetization)?;
        next.m2.push(magnetization.powi(2))?;
        next.m4.push(magnetization.powi(4))?;
        *self = next;
        Ok(())
    }
    pub fn summary(&self, beta: f64, k_b: f64, sites: usize) -> Result<ThermalSummary> {
        if !beta.is_finite() || !k_b.is_finite() || k_b <= 0.0 || sites == 0 || self.energy.n == 0 {
            return Err(Error("invalid units, site count or empty observations"));
        }
        let v = self.energy.population_variance().unwrap();
        let heat = if v == 0.0 || beta == 0.0 {
            0.0
        } else {
            k_b * beta * beta * v / sites as f64
        };
        let chi = beta * self.magnetization.population_variance().unwrap() / sites as f64;
        let binder = if self.m2.mean > 0.0 {
            Some(1.0 - (self.m4.mean / self.m2.mean) / self.m2.mean / 3.0)
        } else {
            None
        };
        if !heat.is_finite() || !chi.is_finite() || binder.is_some_and(|x| !x.is_finite()) {
            return Err(Error("thermodynamic overflow"));
        }
        Ok(ThermalSummary {
            energy: self.energy.mean / sites as f64,
            magnetization: self.magnetization.mean / sites as f64,
            heat_capacity: heat,
            susceptibility: chi,
            binder,
            samples: self.energy.n,
        })
    }
}

/// Delete-one-block jackknife standard errors for nonlinear thermal summaries.
/// Blocks must exceed the correlation scale. At least eight complete blocks;
/// trailing observations are excluded from both estimates and errors.
#[derive(Clone, Copy, Debug)]
pub struct ThermalErrors {
    pub energy: f64,
    pub magnetization: f64,
    pub heat_capacity: f64,
    pub susceptibility: f64,
    pub binder: Option<f64>,
    pub blocks: usize,
    pub used_samples: usize,
}
pub fn thermal_errors(
    samples: &[(f64, f64)],
    beta: f64,
    k_b: f64,
    sites: usize,
    block: usize,
) -> Result<Option<ThermalErrors>> {
    if block == 0 {
        return Err(Error("block size must be positive"));
    }
    let n = samples.len() / block;
    let mut validation = Thermodynamics::default();
    for (e, m) in samples {
        validation.push(*e, *m)?;
    }
    if !samples.is_empty() {
        validation.summary(beta, k_b, sites)?;
    } else if !beta.is_finite() || !k_b.is_finite() || k_b <= 0.0 || sites == 0 {
        return Err(Error("invalid units"));
    }
    if n < 8 {
        return Ok(None);
    }
    if n > 4096 {
        return Err(Error("use larger jackknife blocks (at most 4096 blocks)"));
    }
    let blocks: Vec<_> = samples[..n * block]
        .chunks_exact(block)
        .map(|chunk| {
            let mut t = Thermodynamics::default();
            for (e, m) in chunk {
                t.push(*e, *m)?;
            }
            Ok(t)
        })
        .collect::<Result<_>>()?;
    let mut stats = [Moments::default(); 5];
    let mut binder_defined = true;
    for omitted in 0..n {
        let mut t = Thermodynamics::default();
        for (i, b) in blocks.iter().enumerate() {
            if i != omitted {
                t.energy.merge(&b.energy)?;
                t.magnetization.merge(&b.magnetization)?;
                t.m2.merge(&b.m2)?;
                t.m4.merge(&b.m4)?;
            }
        }
        let s = t.summary(beta, k_b, sites)?;
        for (m, x) in stats[..4].iter_mut().zip([
            s.energy,
            s.magnetization,
            s.heat_capacity,
            s.susceptibility,
        ]) {
            m.push(x)?;
        }
        if let Some(b) = s.binder {
            stats[4].push(b)?;
        } else {
            binder_defined = false;
        }
    }
    let se = |i: usize| (stats[i].m2 * ((n - 1) as f64 / n as f64)).sqrt();
    // A constant observed energy and magnetization cannot establish mixing.
    if blocks.iter().all(|b| {
        b.energy.m2 == 0.0
            && b.magnetization.m2 == 0.0
            && b.energy.mean == blocks[0].energy.mean
            && b.magnetization.mean == blocks[0].magnetization.mean
    }) {
        return Ok(None);
    }
    Ok(Some(ThermalErrors {
        energy: se(0),
        magnetization: se(1),
        heat_capacity: se(2),
        susceptibility: se(3),
        binder: binder_defined.then(|| se(4)),
        blocks: n,
        used_samples: n * block,
    }))
}
