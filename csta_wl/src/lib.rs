//! Checked Wang–Landau warm-up followed by SAMC, and replica-exchange windows.
//!
//! The SAMC counter is **trial moves**, following arXiv:2402.05653v2 Eq. (5).
//! `State` proposals default to symmetric; asymmetric models supply a log proposal
//! ratio. `apply_change`/`revert_change` must be exact
//! inverses and model methods must terminate. Energy may cache into `Params`;
//! rejected moves restore a cloned parameter snapshot. Clone must be independent
//! of mutable caches (no shared interior-mutability side effects).

#![doc = include_str!("../README.md")]

pub use csta_metropolis::State;
pub use csta_montecarlo::Randomizable;
use rand::RngExt;
use std::sync::atomic::{AtomicBool, Ordering};

mod grid;
pub mod models;
mod par_wl;
mod wl_data;
pub use grid::*;
pub use par_wl::*;
pub use wl_data::*;

pub type Result<T> = std::result::Result<T, Error>;
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    Invalid(&'static str),
    InvalidEnergy,
    InitialStateOutside,
    CounterOverflow,
    Numerical(&'static str),
    InsufficientOverlap,
    WorkerFailed,
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(s) | Self::Numerical(s) => f.write_str(s),
            Self::InvalidEnergy => f.write_str("model returned nonfinite energy"),
            Self::InitialStateOutside => f.write_str("initial state is outside accessible support"),
            Self::CounterOverflow => f.write_str("visit counter or budget overflow"),
            Self::InsufficientOverlap => {
                f.write_str("DOS fragments need visited overlap and valid slopes")
            }
            Self::WorkerFailed => f.write_str("parallel worker panicked"),
        }
    }
}
impl std::error::Error for Error {}

#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug)]
pub struct RawWangLandauData {
    grid: EnergyGrid,

    /// ln(g), finite on accessible bins and -infinity elsewhere.
    #[cfg_attr(
        feature = "checkpoint",
        serde(with = "csta_core::checkpoint::float_bits")
    )]
    dos: Vec<f64>,
    bins: Vec<u64>,
    lifetime_bins: Vec<u64>,
    accessible: Vec<bool>,
    visits: u64,
    total_visits: u64,
}
impl RawWangLandauData {
    pub fn new(n_bins: usize, min: f64, max: f64) -> Result<Self> {
        Self::on_grid(EnergyGrid::continuous(min, max, n_bins)?)
    }
    pub fn on_grid(grid: EnergyGrid) -> Result<Self> {
        let n = grid.len();
        Self::with_support(grid, vec![true; n], vec![0.0; n])
    }
    pub fn with_support(grid: EnergyGrid, accessible: Vec<bool>, dos: Vec<f64>) -> Result<Self> {
        let n = grid.len();
        if accessible.len() != n || dos.len() != n || !accessible.iter().any(|v| *v) {
            return Err(Error::Invalid(
                "grid, mask and DOS must have equal nonzero lengths and accessible support",
            ));
        }
        if dos.iter().zip(&accessible).any(|(v, a)| {
            if *a {
                !v.is_finite()
            } else {
                *v != f64::NEG_INFINITY
            }
        }) {
            return Err(Error::Invalid(
                "DOS must be finite on support and -infinity outside it",
            ));
        }
        Ok(Self {
            grid,
            dos,
            accessible,
            bins: vec![0; n],
            lifetime_bins: vec![0; n],
            visits: 0,
            total_visits: 0,
        })
    }
    pub fn from_g_e_0(dos: Vec<f64>, min: f64, max: f64) -> Result<Self> {
        Self::with_support(
            EnergyGrid::continuous(min, max, dos.len())?,
            vec![true; dos.len()],
            dos,
        )
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
    pub fn lifetime_bins(&self) -> &[u64] {
        &self.lifetime_bins
    }
    pub fn accessible(&self) -> &[bool] {
        &self.accessible
    }
    pub fn visits(&self) -> u64 {
        self.visits
    }
    pub fn total_visits(&self) -> u64 {
        self.total_visits
    }
    pub fn active_bins(&self) -> usize {
        self.accessible.iter().filter(|v| **v).count()
    }
    pub fn energy_to_bin(&self, energy: f64) -> Result<Option<usize>> {
        Ok(self.grid.bin(energy)?.filter(|i| self.accessible[*i]))
    }
    pub fn get(&self, energy: f64) -> Result<Option<u64>> {
        Ok(self.energy_to_bin(energy)?.map(|i| self.bins[i]))
    }
    pub fn mean(&self) -> f64 {
        self.bins
            .iter()
            .zip(&self.accessible)
            .filter(|(_, a)| **a)
            .map(|(v, _)| *v as f64)
            .sum::<f64>()
            / self.active_bins() as f64
    }
    pub fn is_flat(&self) -> bool {
        self.flat_at(0.8)
    }
    pub fn flat_at(&self, fraction: f64) -> bool {
        let mean = self.mean();
        fraction.is_finite()
            && fraction > 0.0
            && fraction <= 1.0
            && mean > 0.0
            && self
                .bins
                .iter()
                .zip(&self.accessible)
                .all(|(v, a)| !a || *v as f64 >= fraction * mean)
    }
    pub fn clear_hist(&mut self) {
        self.bins.fill(0);
        self.visits = 0;
    }
    fn check_record(&self, bin: usize, delta: f64) -> Result<()> {
        if !delta.is_finite() || delta < 0.0 {
            return Err(Error::Invalid("DOS update must be finite and nonnegative"));
        }
        if self.total_visits == u64::MAX
            || self.visits == u64::MAX
            || self.bins[bin] == u64::MAX
            || self.lifetime_bins[bin] == u64::MAX
        {
            return Err(Error::CounterOverflow);
        }
        if !(self.dos[bin] + delta).is_finite() {
            return Err(Error::Numerical("DOS update overflow"));
        }
        if delta > 0.0 && self.dos[bin] + delta == self.dos[bin] {
            return Err(Error::Numerical(
                "DOS update below floating-point resolution; anchor initial DOS",
            ));
        }
        Ok(())
    }
    fn record(&mut self, bin: usize, delta: f64) {
        self.dos[bin] += delta;
        self.bins[bin] += 1;
        self.lifetime_bins[bin] += 1;
        self.visits += 1;
        self.total_visits += 1;
    }
    pub fn process_data(self) -> Result<WLData> {
        WLData::new(self.grid, self.dos, self.bins)
    }
}

/// Zero preliminary stages explicitly skips warm-up. Other minima must be
/// positive. `sampling_steps` counts only SAMC proposals; zero does no work,
/// including no warm-up. `max_steps` bounds warm-up + SAMC proposals together.
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub preliminary_stages: usize,
    pub min_stage_steps: u64,
    pub visits_per_bin: u64,
    pub flatness: f64,
    pub initial_ln_f: f64,
    pub sampling_steps: u64,
    pub max_steps: u64,
    /// Defaults are t0 = number of accessible bins, t1 = 10*t0.
    pub t0: Option<f64>,
    pub t1: Option<f64>,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            preliminary_stages: 5,
            min_stage_steps: 1000,
            visits_per_bin: 500,
            flatness: 0.8,
            initial_ln_f: 1.0,
            sampling_steps: 100_000,
            max_steps: 10_000_000,
            t0: None,
            t1: None,
        }
    }
}
impl Config {
    fn validate(&self, bins: usize) -> Result<(u64, f64, f64)> {
        if self.min_stage_steps == 0
            || self.visits_per_bin == 0
            || self.max_steps == 0
            || !self.flatness.is_finite()
            || !(0.0 < self.flatness && self.flatness <= 1.0)
            || !self.initial_ln_f.is_finite()
            || self.initial_ln_f <= 0.0
            || self.preliminary_stages > 1024
            || (self.preliminary_stages > 0
                && self.initial_ln_f * 0.5_f64.powi((self.preliminary_stages - 1) as i32) == 0.0)
        {
            return Err(Error::Invalid(
                "invalid run budget, stage count, flatness or update",
            ));
        }
        let visits = self
            .visits_per_bin
            .checked_mul(bins as u64)
            .ok_or(Error::CounterOverflow)?;
        let t0 = self.t0.unwrap_or(bins as f64);
        let t1 = self.t1.unwrap_or(10.0 * t0);
        if !t0.is_finite()
            || !t1.is_finite()
            || t0 <= 0.0
            || t1 <= 0.0
            || !(t0 / t1).is_finite()
            || t0 / (t1 + self.sampling_steps as f64) <= 0.0
        {
            return Err(Error::Invalid(
                "SAMC scales must give finite positive updates",
            ));
        }
        Ok((visits.max(self.min_stage_steps), t0, t1))
    }
}
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Phase {
    Preliminary,
    Sampling,
}
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopReason {
    TargetReached,
    BudgetExhausted,
    Cancelled,
}
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct Diagnostics {
    pub stop_reason: Option<StopReason>,
    pub phase: Phase,
    pub completed_stages: usize,
    pub proposals: u64,
    pub accepted: u64,
    pub sampling_steps: u64,
    pub last_update: f64,
    pub flat: bool,
    /// Parallel-only moves after this window finishes adapting; excluded from DOS/histograms.
    pub mixing_steps: u64,
}
#[derive(Debug)]
pub struct RunResult<S: State> {
    pub state: S,
    pub params: S::Params,
    pub data: RawWangLandauData,
    pub diagnostics: Diagnostics,
}
impl<S: State> RunResult<S> {
    pub fn is_complete(&self) -> bool {
        self.diagnostics.stop_reason == Some(StopReason::TargetReached)
    }
}

#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "checkpoint",
    serde(bound(
        serialize = "S: serde::Serialize, S::Params: serde::Serialize, R: serde::Serialize",
        deserialize = "S: serde::Deserialize<'de>, S::Params: serde::Deserialize<'de>, R: serde::Deserialize<'de>"
    ))
)]
pub(crate) struct Walker<S: State, R> {
    state: S,
    params: S::Params,
    rng: R,
    energy: f64,
    bin: usize,
    data: RawWangLandauData,
    config: Config,
    minimum_visits: u64,
    t0: f64,
    t1: f64,
    ln_f: f64,
    diagnostics: Diagnostics,
}
impl<S: State, R: RngExt> Walker<S, R>
where
    S::Params: Clone,
{
    fn new(
        state: S,
        mut params: S::Params,
        rng: R,
        data: RawWangLandauData,
        config: Config,
    ) -> Result<Self> {
        let (minimum_visits, t0, t1) = config.validate(data.active_bins())?;
        if data.total_visits != 0 {
            return Err(Error::Invalid(
                "start from fresh counters; supply initial DOS for refinement",
            ));
        }
        let energy = state.energy(&mut params);
        let bin = data
            .energy_to_bin(energy)?
            .ok_or(Error::InitialStateOutside)?;
        let phase = if config.preliminary_stages == 0 {
            Phase::Sampling
        } else {
            Phase::Preliminary
        };
        let stop_reason = (config.sampling_steps == 0).then_some(StopReason::TargetReached);
        Ok(Self {
            state,
            params,
            rng,
            energy,
            bin,
            data,
            minimum_visits,
            t0,
            t1,
            ln_f: config.initial_ln_f,
            config,
            diagnostics: Diagnostics {
                stop_reason,
                phase,
                completed_stages: 0,
                proposals: 0,
                accepted: 0,
                sampling_steps: 0,
                last_update: 0.0,
                flat: false,
                mixing_steps: 0,
            },
        })
    }
    /// Atomic wrt model energy errors: revert the change and restore Params.
    fn transition(&mut self, delta: Option<f64>) -> Result<bool> {
        let backup = self.params.clone();
        let change = self.state.propose_change(&mut self.rng);
        let correction = self.state.log_proposal_ratio(&change);
        if correction.is_nan() || correction == f64::INFINITY {
            return Err(Error::Invalid("invalid proposal log ratio"));
        }
        let local = self.state.delta_energy(&change, &self.params);
        self.state.apply_change(change.clone());
        let energy = local.map_or_else(|| self.state.energy(&mut self.params), |d| self.energy + d);
        let new_bin = match self.data.energy_to_bin(energy) {
            Ok(bin) => bin,
            Err(e) => {
                self.state.revert_change(change);
                self.params = backup;
                return Err(e);
            }
        };
        let accepted = new_bin.is_some_and(|i| {
            accept(
                if correction == f64::NEG_INFINITY {
                    correction
                } else {
                    self.data.dos[self.bin] - self.data.dos[i] + correction
                },
                &mut self.rng,
            )
        });
        // Validate the retained-bin update before committing any model changes.
        if let Some(d) = delta {
            let retained = if accepted { new_bin.unwrap() } else { self.bin };
            if let Err(error) = self.data.check_record(retained, d) {
                self.state.revert_change(change);
                self.params = backup;
                return Err(error);
            }
        }
        if accepted {
            self.energy = energy;
            self.bin = new_bin.unwrap();
        } else {
            self.state.revert_change(change);
            self.params = backup;
        }
        if let Some(d) = delta {
            self.data.record(self.bin, d);
        }
        Ok(accepted)
    }
    fn advance(&mut self, cancel: Option<&AtomicBool>) -> Result<()> {
        if self.diagnostics.stop_reason.is_some() {
            return Ok(());
        }
        if cancelled(cancel) {
            self.diagnostics.stop_reason = Some(StopReason::Cancelled);
            return Ok(());
        }
        if self.diagnostics.proposals >= self.config.max_steps {
            self.diagnostics.stop_reason = Some(StopReason::BudgetExhausted);
            return Ok(());
        }
        // Eq. (5): t is completed SAMC trial moves, not visits divided by bins.
        let delta = match self.diagnostics.phase {
            Phase::Preliminary => self.ln_f,
            Phase::Sampling => self.t0 / (self.t1 + self.diagnostics.sampling_steps as f64),
        };
        let accepted = self.transition(Some(delta))?;
        self.diagnostics.proposals += 1;
        self.diagnostics.accepted += u64::from(accepted);
        self.diagnostics.last_update = delta;
        match self.diagnostics.phase {
            Phase::Preliminary => {
                if self.data.visits >= self.minimum_visits
                    && self.data.flat_at(self.config.flatness)
                {
                    self.diagnostics.completed_stages += 1;
                    self.ln_f *= 0.5;
                    self.data.clear_hist();
                    self.diagnostics.flat = false;
                    if self.diagnostics.completed_stages == self.config.preliminary_stages {
                        self.diagnostics.phase = Phase::Sampling;
                    }
                }
            }
            Phase::Sampling => {
                self.diagnostics.sampling_steps += 1;
                if self.diagnostics.sampling_steps == self.config.sampling_steps {
                    self.diagnostics.stop_reason = Some(StopReason::TargetReached);
                }
            }
        }
        if self.diagnostics.stop_reason.is_none()
            && self.diagnostics.proposals == self.config.max_steps
        {
            self.diagnostics.stop_reason = Some(StopReason::BudgetExhausted);
        }
        Ok(())
    }
    fn finish(mut self) -> RunResult<S> {
        self.diagnostics.flat = self.data.flat_at(self.config.flatness);
        RunResult {
            state: self.state,
            params: self.params,
            data: self.data,
            diagnostics: self.diagnostics,
        }
    }
}
fn accept(log_ratio: f64, rng: &mut impl RngExt) -> bool {
    // log comparison avoids exp overflow; equality rejects (u is in [0,1)).
    log_ratio >= 0.0 || rng.random::<f64>().ln() < log_ratio
}
fn cancelled(flag: Option<&AtomicBool>) -> bool {
    flag.is_some_and(|f| f.load(Ordering::Relaxed))
}

/// Preferred entry point. A supplied initial state must be inside the support.
/// Cancellation and budget exhaustion return partial data with an explicit reason.
/// Model/data errors return Err; no result is presented as complete.
pub fn run<S: State, R: RngExt>(
    state: S,
    params: S::Params,
    rng: R,
    data: RawWangLandauData,
    config: Config,
    cancel: Option<&AtomicBool>,
) -> Result<RunResult<S>>
where
    S::Params: Clone,
{
    let mut walker = Walker::new(state, params, rng, data, config)?;
    while walker.diagnostics.stop_reason.is_none() {
        walker.advance(cancel)?;
    }
    Ok(walker.finish())
}

/// Bounded random initialization; use a supplied state for rare windows.
pub fn sample_in_support<S: State + Randomizable>(
    rng: &mut impl RngExt,
    params: &S::Params,
    data: &RawWangLandauData,
    attempts: usize,
) -> Result<S>
where
    S::Params: Clone,
{
    if attempts == 0 {
        return Err(Error::Invalid("initialization attempts must be positive"));
    }
    for _ in 0..attempts {
        let state = S::sample(rng);
        if data
            .energy_to_bin(state.energy(&mut params.clone()))?
            .is_some()
        {
            return Ok(state);
        }
    }
    Err(Error::InitialStateOutside)
}
fn target_steps(time: f64, bins: usize) -> Result<u64> {
    let steps = time * bins as f64;
    if bins == 0
        || !time.is_finite()
        || time < 0.0
        || !steps.is_finite()
        || steps >= u64::MAX as f64
    {
        return Err(Error::Invalid(
            "target time must produce a finite nonnegative proposal count",
        ));
    }
    Ok(steps.ceil() as u64)
}
/// Compatibility-shaped convenience wrapper. Returns a checked run result.
/// `target_time` is SAMC visits per bin; warm-up is separate and bounded.
#[allow(clippy::too_many_arguments)]
pub fn wang_landau2<S: State + Randomizable>(
    target_time: f64,
    preliminary_sweeps: usize,
    preliminary_runs: usize,
    k_var: usize,
    params: S::Params,
    min: f64,
    max: f64,
    bins: usize,
) -> Result<RunResult<S>>
where
    S::Params: Clone,
{
    let config = Config {
        preliminary_stages: preliminary_runs,
        min_stage_steps: preliminary_sweeps as u64,
        visits_per_bin: k_var as u64,
        sampling_steps: target_steps(target_time, bins)?,
        ..Config::default()
    };
    config.validate(bins)?;
    let data = RawWangLandauData::new(bins, min, max)?;
    let mut rng = rand::rng();
    let state = sample_in_support::<S>(&mut rng, &params, &data, 1000)?;
    run(state, params, rng, data, config, None)
}
/// The previous 0.3 switch threshold corresponds to two halvings from ln_f=1.
pub fn wang_landau<S: State + Randomizable>(
    target_time: f64,
    sweeps: usize,
    params: S::Params,
    min: f64,
    max: f64,
    bins: usize,
) -> Result<RunResult<S>>
where
    S::Params: Clone,
{
    wang_landau2::<S>(target_time, sweeps, 2, 500, params, min, max, bins)
}

pub mod analysis;
#[cfg(test)]
mod analysis_tests;
pub mod joint;
mod session;
#[cfg(test)]
mod tests;
pub use session::Session;
