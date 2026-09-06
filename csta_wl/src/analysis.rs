//! Observable reweighting and uncertainty from explicitly independent runs.
use crate::{EnergyGrid, Error, Result, WLData, log_sum_exp};
use csta_metropolis::statistics::Moments;

#[derive(Clone, Debug)]
pub struct ConditionalMoments {
    grid: EnergyGrid,
    bins: Vec<Moments>,
    squares: Vec<Moments>,
}
#[derive(Clone, Copy, Debug)]
pub struct ObservableMoments {
    pub mean: f64,
    pub second_moment: f64,
}
impl ConditionalMoments {
    pub fn new(grid: EnergyGrid) -> Self {
        Self {
            bins: vec![Moments::default(); grid.len()],
            squares: vec![Moments::default(); grid.len()],
            grid,
        }
    }
    /// Record a retained state, including repeated states after rejection. Only
    /// use samples whose within-bin conditional distribution is the desired one.
    pub fn record(&mut self, energy: f64, value: f64) -> Result<()> {
        let i = self
            .grid
            .bin(energy)?
            .ok_or(Error::Invalid("observation outside grid"))?;
        let mut m = self.bins[i];
        let mut s = self.squares[i];
        m.push(value)
            .map_err(|_| Error::Numerical("invalid observable"))?;
        s.push(value * value)
            .map_err(|_| Error::Numerical("observable square overflow"))?;
        self.bins[i] = m;
        self.squares[i] = s;
        Ok(())
    }
    pub fn evaluate(&self, dos: &WLData, beta: f64) -> Result<ObservableMoments> {
        if self.grid != *dos.grid() {
            return Err(Error::Invalid("conditional data grid mismatch"));
        }
        let probabilities = dos.energy_distribution(beta)?;
        let mut mean = 0.0;
        let mut second = 0.0;
        for (i, p) in probabilities.iter().enumerate() {
            // Missing finite support is an error even if weights underflow at this beta.
            if dos.dos()[i].is_finite() {
                mean += p * self.bins[i]
                    .mean()
                    .ok_or(Error::Invalid("missing conditional support"))?;
                second += p * self.squares[i].mean().unwrap();
            }
        }
        if !mean.is_finite() || !second.is_finite() {
            return Err(Error::Numerical("observable overflow"));
        }
        Ok(ObservableMoments {
            mean,
            second_moment: second,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Overlap {
    Adequate,
    Low,
}
#[derive(Clone, Copy, Debug)]
pub struct Reweighted {
    pub mean: f64,
    pub weight_effective_samples: f64,
    pub overlap: Overlap,
}
/// Canonical single-histogram importance reweighting. Weight ESS detects weight
/// concentration, not autocorrelation or unseen states; Adequate is not proof of
/// overlap. The caller must establish stationary source data and physical support.
pub fn reweight(samples: &[(f64, f64)], source_beta: f64, target_beta: f64) -> Result<Reweighted> {
    if samples.is_empty()
        || !source_beta.is_finite()
        || !target_beta.is_finite()
        || samples
            .iter()
            .any(|(e, o)| !e.is_finite() || !o.is_finite())
    {
        return Err(Error::Invalid("invalid reweighting input"));
    }
    let reference = samples[0].0;
    let logs: Vec<_> = samples
        .iter()
        .map(|(e, _)| {
            if source_beta == target_beta {
                0.0
            } else {
                (source_beta - target_beta) * (e - reference)
            }
        })
        .collect();
    if logs.iter().any(|x| !x.is_finite()) {
        return Err(Error::Numerical("reweighting log overflow"));
    }
    let z = log_sum_exp(&logs)?;
    let mut mean = 0.0;
    let mut sum2 = 0.0;
    for ((_, value), w) in samples.iter().zip(logs) {
        let p = (w - z).exp();
        mean += p * value;
        sum2 += p * p;
    }
    if !mean.is_finite() {
        return Err(Error::Numerical("reweighted observable overflow"));
    }
    let ess = 1.0 / sum2;
    Ok(Reweighted {
        mean,
        weight_effective_samples: ess,
        overlap: if ess < 20.0 || ess < 0.1 * samples.len() as f64 {
            Overlap::Low
        } else {
            Overlap::Adequate
        },
    })
}

/// One completed independent experiment. IDs prevent accidentally including the
/// same run twice; callers must not label interacting walkers as independent runs.
pub struct IndependentDos {
    id: u64,
    data: WLData,
}
impl IndependentDos {
    pub fn from_run<S: crate::State>(id: u64, run: crate::RunResult<S>) -> Result<Self> {
        if !run.is_complete() || run.diagnostics.sampling_steps == 0 {
            return Err(Error::Invalid("incomplete independent run"));
        }
        Self::from_data(id, run.data.process_data()?)
    }
    /// For externally completed or exact reference estimates. This assertion of
    /// independence/completion belongs to the caller.
    pub fn from_data(id: u64, data: WLData) -> Result<Self> {
        Ok(Self { id, data })
    }
}
pub struct DosEnsemble {
    runs: Vec<WLData>,
}
#[derive(Clone, Copy, Debug)]
pub struct EnsembleEstimate {
    pub mean: f64,
    pub standard_error: Option<f64>,
    pub runs: usize,
}
impl DosEnsemble {
    /// Normalize every independent estimate to the same total state count.
    pub fn new(runs: Vec<IndependentDos>, log_count: f64) -> Result<Self> {
        if runs.is_empty() {
            return Err(Error::Invalid("empty DOS ensemble"));
        }
        let mut ids = std::collections::HashSet::new();
        let reference = &runs[0].data;
        for r in &runs {
            if !ids.insert(r.id)
                || r.data.grid() != reference.grid()
                || r.data
                    .dos()
                    .iter()
                    .zip(reference.dos())
                    .any(|(a, b)| a.is_finite() != b.is_finite())
            {
                return Err(Error::Invalid("duplicate run or incompatible DOS support"));
            }
        }
        let mut data = Vec::with_capacity(runs.len());
        for mut r in runs {
            r.data.normalize_log_count(log_count)?;
            data.push(r.data);
        }
        Ok(Self { runs: data })
    }
    pub fn evaluate(
        &self,
        mut observable: impl FnMut(&WLData) -> Result<f64>,
    ) -> Result<EnsembleEstimate> {
        let mut moments = Moments::default();
        for r in &self.runs {
            moments
                .push(observable(r)?)
                .map_err(|_| Error::Numerical("invalid ensemble observable"))?;
        }
        Ok(EnsembleEstimate {
            mean: moments.mean().unwrap(),
            standard_error: moments
                .variance()
                .map(|v| (v / self.runs.len() as f64).sqrt()),
            runs: self.runs.len(),
        })
    }
    pub fn energy(&self, beta: f64) -> Result<EnsembleEstimate> {
        self.evaluate(|d| Ok(d.energy_moments(beta)?.0))
    }
    /// Per-bin scatter exposes disagreements, including energy-window overlaps.
    pub fn log_dos_uncertainty(&self) -> Result<Vec<Option<EnsembleEstimate>>> {
        (0..self.runs[0].grid().len())
            .map(|i| {
                if self.runs[0].dos()[i].is_finite() {
                    self.evaluate(|d| Ok(d.dos()[i])).map(Some)
                } else {
                    Ok(None)
                }
            })
            .collect()
    }
}
