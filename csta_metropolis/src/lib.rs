//! This module is for metropoli + montecarlo simulations

use crate::observer::*;
use rand::{RngExt, rngs::ThreadRng};

pub mod models;
pub mod observer;
pub mod replica;
pub mod statistics;

pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, Clone, PartialEq)]
pub struct Error(pub &'static str);
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for Error {}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunStatus {
    Complete,
    Cancelled,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RunReport {
    pub attempted: usize,
    pub accepted: usize,
    pub status: RunStatus,
}
/// Measurements follow completed moves: burn-in moves are skipped, then every
/// `stride`th retained move is measured. The schedule restarts for each run.
#[derive(Clone, Copy, Debug)]
pub struct Schedule {
    pub burn_in: usize,
    pub stride: usize,
}
impl Default for Schedule {
    fn default() -> Self {
        Self {
            burn_in: 0,
            stride: 1,
        }
    }
}
impl Schedule {
    pub fn validate(self) -> Result<()> {
        if self.stride == 0 {
            Err(Error("measurement stride must be positive"))
        } else {
            Ok(())
        }
    }
    pub fn measures(self, completed: usize) -> bool {
        completed > self.burn_in && (completed - self.burn_in).is_multiple_of(self.stride)
    }
}
/// Valid finite energies and beta; overflow of a log ratio has the limiting
/// acceptance meaning. Beta zero is handled before subtracting extreme energies.
pub fn boltzmann_log_ratio(beta: f64, old: f64, new: f64) -> Result<f64> {
    if !beta.is_finite() || !old.is_finite() || !new.is_finite() {
        return Err(Error("beta and energies must be finite"));
    }
    if beta == 0.0 || old == new {
        return Ok(0.0);
    }
    let delta = old - new;
    Ok(if delta.is_finite() {
        beta * delta
    } else {
        (beta * old) - (beta * new)
    })
}
pub fn accept_log(log_ratio: f64, rng: &mut impl RngExt) -> bool {
    log_ratio >= 0.0 || rng.random::<f64>().ln() < log_ratio
}

pub trait State {
    type Params;
    type Change: Clone;
    /// State changes must be reversible; model methods must terminate.
    fn energy(&self, params: &mut Self::Params) -> f64;
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change;
    fn apply_change(&mut self, change: Self::Change) /* -> ModificationError */;
    fn revert_change(&mut self, change: Self::Change);
    /// Override for persisted models with structural invariants. Energy is also
    /// checked by checkpoint loaders. No extra method is required for ordinary runs.
    fn valid_state(&self, _params: &Self::Params) -> bool {
        true
    }
    /// Optional log(q(reverse)/q(forward)), evaluated before applying change.
    /// Negative infinity means no reverse move. Positive infinity and NaN are invalid.
    fn log_proposal_ratio(&self, _change: &Self::Change) -> f64 {
        0.0
    }
    /// Optional local energy difference. None uses full recomputation.
    /// Override only when accepted changes do not require refreshing Params caches.
    fn delta_energy(&self, _change: &Self::Change, _params: &Self::Params) -> Option<f64> {
        None
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
#[derive(Clone)]
pub struct Metropolis<S: State, R: RngExt> {
    pub state: S,
    pub params: S::Params,
    pub beta: f64,
    pub steps: usize,
    pub accepted_moves: usize,
    pub attempted_moves: usize,
    pub rng: R,
}

impl<S> Metropolis<S, ThreadRng>
where
    S: State + Default,
    S::Params: Default,
{
    pub fn with_steps(beta: f64, steps: usize) -> Self {
        Self::with_rng(beta, steps, rand::rng())
    }

    pub fn with_steps_no_beta(steps: usize) -> Self {
        Self::with_steps(1.0, steps)
    }
}

impl<S> Metropolis<S, ThreadRng>
where
    S: State,
    S::Params: Default,
{
    pub fn with_state(state: S, beta: f64, steps: usize) -> Self {
        Self::with_state_params(state, S::Params::default(), beta, steps)
    }

    pub fn with_state_no_beta(state: S, steps: usize) -> Self {
        Self::with_state(state, 1.0, steps)
    }
}

impl<S: State> Metropolis<S, ThreadRng> {
    pub fn with_state_params(state: S, params: S::Params, beta: f64, steps: usize) -> Self {
        Self::with_all(state, params, beta, steps, rand::rng())
    }

    pub fn with_state_params_no_beta(state: S, params: S::Params, steps: usize) -> Self {
        Self::with_state_params(state, params, 1.0, steps)
    }
}

impl<S, R> Metropolis<S, R>
where
    S: State,
    S::Params: Default,
    R: RngExt,
{
    pub fn with_state_rng(state: S, beta: f64, steps: usize, rng: R) -> Self {
        Self::with_all(state, S::Params::default(), beta, steps, rng)
    }

    pub fn with_state_rng_no_beta(state: S, steps: usize, rng: R) -> Self {
        Self::with_state_rng(state, 1.0, steps, rng)
    }
}

impl<S, R> Metropolis<S, R>
where
    S: State + Default,
    S::Params: Default,
    R: RngExt,
{
    pub fn with_rng(beta: f64, steps: usize, rng: R) -> Self {
        Self::with_all(S::default(), S::Params::default(), beta, steps, rng)
    }

    pub fn with_rng_no_beta(steps: usize, rng: R) -> Self {
        Self::with_rng(1.0, steps, rng)
    }
}

impl<S: State, R: RngExt> Metropolis<S, R> {
    pub fn with_all(state: S, params: S::Params, beta: f64, steps: usize, rng: R) -> Self {
        Self {
            state,
            params,
            beta,
            steps,
            accepted_moves: 0,
            attempted_moves: 0,
            rng,
        }
    }

    pub fn with_all_no_beta(state: S, params: S::Params, steps: usize, rng: R) -> Self {
        Self::with_all(state, params, 1.0, steps, rng)
    }

    pub fn accepted_rate(&self) -> f64 {
        if self.attempted_moves == 0 {
            0.0
        } else {
            self.accepted_moves as f64 / self.attempted_moves as f64
        }
    }
    pub fn rejected_rate(&self) -> f64 {
        if self.attempted_moves == 0 {
            0.0
        } else {
            1.0 - self.accepted_rate()
        }
    }
}

impl<S: State, R: RngExt> Metropolis<S, R>
where
    S::Params: Clone,
{
    /// Checked move. Params::clone must not share mutable caches. Errors restore
    /// model/cache/counters, but RNG draws already consumed are not rewound.
    pub fn try_step(&mut self) -> Result<bool> {
        if !self.beta.is_finite() {
            return Err(Error("beta must be finite"));
        }
        if self.attempted_moves == usize::MAX
            || self.accepted_moves == usize::MAX
            || self.accepted_moves > self.attempted_moves
        {
            return Err(Error("invalid or overflowing move counters"));
        }
        let backup = self.params.clone();
        let old = self.state.energy(&mut self.params);
        if !old.is_finite() {
            self.params = backup;
            return Err(Error("nonfinite initial energy"));
        }
        let change = self.state.propose_change(&mut self.rng);
        let correction = self.state.log_proposal_ratio(&change);
        if correction.is_nan() || correction == f64::INFINITY {
            self.params = backup;
            return Err(Error("invalid proposal log ratio"));
        }
        let local = self.state.delta_energy(&change, &self.params);
        self.state.apply_change(change.clone());
        let new = local.map_or_else(|| self.state.energy(&mut self.params), |d| old + d);
        let log_ratio = boltzmann_log_ratio(self.beta, old, new);
        let log_ratio = match log_ratio {
            Ok(r) if correction == f64::NEG_INFINITY => {
                let _ = r;
                f64::NEG_INFINITY
            }
            Ok(r) => r + correction,
            Err(e) => {
                self.state.revert_change(change);
                self.params = backup;
                return Err(e);
            }
        };
        if log_ratio.is_nan() {
            self.state.revert_change(change);
            self.params = backup;
            return Err(Error("indeterminate acceptance ratio"));
        }
        let accepted = accept_log(log_ratio, &mut self.rng);
        if accepted {
            self.accepted_moves += 1;
        } else {
            self.state.revert_change(change);
            self.params = backup;
        }
        self.attempted_moves += 1;
        Ok(accepted)
    }
    /// Convenience wrapper; use try_step to handle invalid numerical models.
    pub fn step(&mut self) {
        self.try_step().expect("invalid Metropolis move");
    }

    pub fn run(&mut self, cancel: Option<&std::sync::atomic::AtomicBool>) -> Result<RunReport> {
        self.run_observed(Schedule::default(), cancel, |_, _| {})
    }
    pub fn run_observed(
        &mut self,
        schedule: Schedule,
        cancel: Option<&std::sync::atomic::AtomicBool>,
        mut observe: impl FnMut(&S, &S::Params),
    ) -> Result<RunReport> {
        self.try_run_observed(schedule, cancel, |s, p| {
            observe(s, p);
            Ok(())
        })
    }
    /// Fallible streaming observer. A measurement error occurs after its move
    /// committed; counters/state remain at that move, and execution stops.
    pub fn try_run_observed(
        &mut self,
        schedule: Schedule,
        cancel: Option<&std::sync::atomic::AtomicBool>,
        mut observe: impl FnMut(&S, &S::Params) -> Result<()>,
    ) -> Result<RunReport> {
        schedule.validate()?;
        if !self.beta.is_finite() {
            return Err(Error("beta must be finite"));
        }
        let start = self.accepted_moves;
        for i in 0..self.steps {
            if cancel.is_some_and(|c| c.load(std::sync::atomic::Ordering::Relaxed)) {
                return Ok(RunReport {
                    attempted: i,
                    accepted: self.accepted_moves - start,
                    status: RunStatus::Cancelled,
                });
            }
            self.try_step()?;
            if schedule.measures(i + 1) {
                observe(&self.state, &self.params)?;
            }
        }
        Ok(RunReport {
            attempted: self.steps,
            accepted: self.accepted_moves - start,
            status: RunStatus::Complete,
        })
    }
    pub fn run_empty(&mut self) {
        self.run(None).expect("invalid Metropolis run");
    }
    pub fn run_with<O: Observer<S>>(&mut self) -> Vec<O::Observation> {
        let mut out = Vec::new();
        self.run_observed(
            Schedule {
                burn_in: O::after(),
                stride: O::every(),
            },
            None,
            |s, p| out.push(O::measure(s, p)),
        )
        .expect("invalid observer run");
        out
    }
    #[allow(clippy::type_complexity)]
    pub fn run_with_2<O1, O2>(&mut self) -> (Vec<O1::Observation>, Vec<O2::Observation>)
    where
        O1: Observer<S>,
        O2: Observer<S>,
    {
        let mut o1 = Vec::new();
        let s1 = Schedule {
            burn_in: O1::after(),
            stride: O1::every(),
        };
        s1.validate().expect("invalid observer stride");
        let mut o2 = Vec::new();
        let s2 = Schedule {
            burn_in: O2::after(),
            stride: O2::every(),
        };
        s2.validate().expect("invalid observer stride");
        for i in 0..self.steps {
            self.step();
            if s1.measures(i + 1) {
                o1.push(O1::measure(&self.state, &self.params));
            }
            if s2.measures(i + 1) {
                o2.push(O2::measure(&self.state, &self.params));
            }
        }
        (o1, o2)
    }
    #[allow(clippy::type_complexity)]
    pub fn run_with_3<O1, O2, O3>(
        &mut self,
    ) -> (
        Vec<O1::Observation>,
        Vec<O2::Observation>,
        Vec<O3::Observation>,
    )
    where
        O1: Observer<S>,
        O2: Observer<S>,
        O3: Observer<S>,
    {
        let mut o1 = Vec::new();
        let s1 = Schedule {
            burn_in: O1::after(),
            stride: O1::every(),
        };
        s1.validate().expect("invalid observer stride");
        let mut o2 = Vec::new();
        let s2 = Schedule {
            burn_in: O2::after(),
            stride: O2::every(),
        };
        s2.validate().expect("invalid observer stride");
        let mut o3 = Vec::new();
        let s3 = Schedule {
            burn_in: O3::after(),
            stride: O3::every(),
        };
        s3.validate().expect("invalid observer stride");
        for i in 0..self.steps {
            self.step();
            if s1.measures(i + 1) {
                o1.push(O1::measure(&self.state, &self.params));
            }
            if s2.measures(i + 1) {
                o2.push(O2::measure(&self.state, &self.params));
            }
            if s3.measures(i + 1) {
                o3.push(O3::measure(&self.state, &self.params));
            }
        }
        (o1, o2, o3)
    }
    #[allow(clippy::type_complexity)]
    pub fn run_with_4<O1, O2, O3, O4>(
        &mut self,
    ) -> (
        Vec<O1::Observation>,
        Vec<O2::Observation>,
        Vec<O3::Observation>,
        Vec<O4::Observation>,
    )
    where
        O1: Observer<S>,
        O2: Observer<S>,
        O3: Observer<S>,
        O4: Observer<S>,
    {
        let mut o1 = Vec::new();
        let s1 = Schedule {
            burn_in: O1::after(),
            stride: O1::every(),
        };
        s1.validate().expect("invalid observer stride");
        let mut o2 = Vec::new();
        let s2 = Schedule {
            burn_in: O2::after(),
            stride: O2::every(),
        };
        s2.validate().expect("invalid observer stride");
        let mut o3 = Vec::new();
        let s3 = Schedule {
            burn_in: O3::after(),
            stride: O3::every(),
        };
        s3.validate().expect("invalid observer stride");
        let mut o4 = Vec::new();
        let s4 = Schedule {
            burn_in: O4::after(),
            stride: O4::every(),
        };
        s4.validate().expect("invalid observer stride");
        for i in 0..self.steps {
            self.step();
            if s1.measures(i + 1) {
                o1.push(O1::measure(&self.state, &self.params));
            }
            if s2.measures(i + 1) {
                o2.push(O2::measure(&self.state, &self.params));
            }
            if s3.measures(i + 1) {
                o3.push(O3::measure(&self.state, &self.params));
            }
            if s4.measures(i + 1) {
                o4.push(O4::measure(&self.state, &self.params));
            }
        }
        (o1, o2, o3, o4)
    }
    pub fn run_with_n<Obs>(
        &mut self,
        obs: Vec<Box<dyn DynObserver<S, Observation = Obs>>>,
    ) -> Vec<Vec<Obs>> {
        let schedules: Vec<_> = obs
            .iter()
            .map(|o| Schedule {
                burn_in: o.after(),
                stride: o.every(),
            })
            .collect();
        for s in &schedules {
            s.validate().expect("invalid observer stride");
        }
        let mut out: Vec<Vec<Obs>> = (0..obs.len()).map(|_| Vec::new()).collect();
        for i in 0..self.steps {
            self.step();
            for (j, o) in obs.iter().enumerate() {
                if schedules[j].measures(i + 1) {
                    out[j].push(o.measure(&self.state, &self.params));
                }
            }
        }
        out
    }
}

/// Add a proposal correction without changing an existing model implementation.
#[derive(Clone)]
pub struct Hastings<S, F> {
    pub state: S,
    pub log_ratio: F,
}
impl<S: State, F: Fn(&S, &S::Change) -> f64> State for Hastings<S, F> {
    type Params = S::Params;
    type Change = S::Change;
    fn energy(&self, p: &mut Self::Params) -> f64 {
        self.state.energy(p)
    }
    fn propose_change(&self, r: &mut impl RngExt) -> Self::Change {
        self.state.propose_change(r)
    }
    fn apply_change(&mut self, c: Self::Change) {
        self.state.apply_change(c)
    }
    fn revert_change(&mut self, c: Self::Change) {
        self.state.revert_change(c)
    }
    fn log_proposal_ratio(&self, c: &Self::Change) -> f64 {
        (self.log_ratio)(&self.state, c)
    }
    fn delta_energy(&self, c: &Self::Change, p: &Self::Params) -> Option<f64> {
        self.state.delta_energy(c, p)
    }
    fn valid_state(&self, p: &Self::Params) -> bool {
        self.state.valid_state(p)
    }
}
#[cfg(test)]
mod tests;

#[cfg(feature = "checkpoint")]
impl<S: State, R: RngExt> Metropolis<S, R>
where
    S: serde::Serialize + serde::de::DeserializeOwned,
    S::Params: Clone + serde::Serialize + serde::de::DeserializeOwned,
    R: serde::Serialize + serde::de::DeserializeOwned,
{
    pub fn save(
        &self,
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.validate_checkpoint()?;
        csta_core::checkpoint::save(path, model, self)
    }
    pub fn load(
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let s: Self = csta_core::checkpoint::load(path, model)?;
        s.validate_checkpoint()?;
        Ok(s)
    }
    fn validate_checkpoint(&self) -> Result<()> {
        if !self.beta.is_finite()
            || self.accepted_moves > self.attempted_moves
            || !self.state.valid_state(&self.params)
            || !self.state.energy(&mut self.params.clone()).is_finite()
        {
            return Err(Error("invalid Metropolis snapshot"));
        }
        Ok(())
    }
}
pub mod finite_size;
