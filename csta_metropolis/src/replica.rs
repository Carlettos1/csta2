//! Canonical temperature exchange with a fixed, sequential exchange order.
use crate::{Error, Metropolis, Result, RunStatus, State, accept_log, boltzmann_log_ratio};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use std::sync::atomic::{AtomicBool, Ordering};

#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "checkpoint",
    serde(bound(
        serialize = "S: serde::Serialize, S::Params: serde::Serialize, R: serde::Serialize",
        deserialize = "S: serde::Deserialize<'de>, S::Params: serde::Deserialize<'de>, R: serde::Deserialize<'de>"
    ))
)]
pub struct ReplicaExchange<S: State, R: RngExt = StdRng> {
    replicas: Vec<Metropolis<S, R>>,
    rng: R,
    local_cursor: usize,
    pending_steps: usize,
    parity: usize,
    labels: Vec<usize>,
    seen_hot: Vec<bool>,
    seen_cold: Vec<bool>,
    pub swap_attempts: u64,
    pub swap_accepts: u64,
    pub round_trips: u64,
}
impl<S: State, R: RngExt + SeedableRng> ReplicaExchange<S, R>
where
    S::Params: Clone,
{
    /// Ladder must be finite, strictly increasing, and nonempty. Same Hamiltonian
    /// parameters in all replicas; Params clones must have independent caches.
    pub fn seeded(states: Vec<S>, params: S::Params, betas: Vec<f64>, seed: u64) -> Result<Self> {
        if states.is_empty()
            || states.len() != betas.len()
            || betas.iter().any(|b| !b.is_finite())
            || betas.windows(2).any(|b| b[0] >= b[1])
        {
            return Err(Error("invalid beta ladder or replica count"));
        }
        let mut rng = R::seed_from_u64(seed);
        let n = states.len();
        let mut replicas = Vec::with_capacity(n);
        for (s, b) in states.into_iter().zip(betas) {
            let mut p = params.clone();
            if !s.energy(&mut p).is_finite() {
                return Err(Error("invalid replica energy"));
            }
            replicas.push(Metropolis::with_all(
                s,
                p,
                b,
                0,
                R::seed_from_u64(rng.random()),
            ));
        }
        let mut seen_cold = vec![false; n];
        seen_cold[n - 1] = true;
        Ok(Self {
            replicas,
            rng,
            local_cursor: 0,
            pending_steps: 0,
            parity: 0,
            labels: (0..n).collect(),
            seen_hot: vec![false; n],
            seen_cold,
            swap_attempts: 0,
            swap_accepts: 0,
            round_trips: 0,
        })
    }
    pub fn replicas(&self) -> &[Metropolis<S, R>] {
        &self.replicas
    }
    pub fn exchange(&mut self, i: usize) -> Result<bool> {
        if i >= self.replicas.len().saturating_sub(1) {
            return Err(Error("invalid exchange pair"));
        }
        if self.swap_attempts == u64::MAX || self.swap_accepts == u64::MAX {
            return Err(Error("swap counter overflow"));
        }
        let (left, right) = self.replicas.split_at_mut(i + 1);
        let a = &mut left[i];
        let b = &mut right[0];
        let mut pa = a.params.clone();
        let mut pb = b.params.clone();
        let ea = a.state.energy(&mut pa);
        let eb = b.state.energy(&mut pb);
        // log acceptance = (beta_i-beta_j)*(E_i-E_j).
        let ratio = boltzmann_log_ratio(a.beta - b.beta, ea, eb)?;
        let accepted = accept_log(ratio, &mut self.rng);
        self.swap_attempts += 1;
        if accepted {
            std::mem::swap(&mut a.state, &mut b.state);
            std::mem::swap(&mut a.params, &mut b.params);
            self.labels.swap(i, i + 1);
            self.swap_accepts += 1;
        }
        Ok(accepted)
    }
    pub fn run(
        &mut self,
        rounds: usize,
        local_steps: usize,
        cancel: Option<&AtomicBool>,
    ) -> Result<RunStatus> {
        if local_steps == 0 {
            return Err(Error("local exchange chunk must be positive"));
        }
        if self.pending_steps != 0 && self.pending_steps != local_steps {
            return Err(Error("resume with the original local chunk size"));
        }
        let total = local_steps
            .checked_mul(self.replicas.len())
            .ok_or(Error("exchange chunk overflow"))?;
        for _ in 0..rounds {
            self.pending_steps = local_steps;
            while self.local_cursor < total {
                if cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                    return Ok(RunStatus::Cancelled);
                }
                self.replicas[self.local_cursor / local_steps].try_step()?;
                self.local_cursor += 1;
            }
            for i in (self.parity..self.replicas.len().saturating_sub(1)).step_by(2) {
                self.exchange(i)?;
            }
            if self.replicas.len() > 1 {
                let hot = self.labels[0];
                let cold = *self.labels.last().unwrap();
                if self.seen_cold[hot] {
                    self.seen_hot[hot] = true;
                }
                if self.seen_hot[cold] {
                    self.round_trips = self
                        .round_trips
                        .checked_add(1)
                        .ok_or(Error("round-trip counter overflow"))?;
                    self.seen_hot[cold] = false;
                }
                self.seen_cold[cold] = true;
            }
            self.parity ^= 1;
            self.local_cursor = 0;
            self.pending_steps = 0;
        }
        Ok(RunStatus::Complete)
    }
}

impl<S: State> ReplicaExchange<S>
where
    S::Params: Clone,
{
    pub fn new(states: Vec<S>, params: S::Params, betas: Vec<f64>, seed: u64) -> Result<Self> {
        Self::seeded(states, params, betas, seed)
    }
}
#[cfg(feature = "checkpoint")]
impl<S: State, R: RngExt + SeedableRng> ReplicaExchange<S, R>
where
    S: serde::Serialize + serde::de::DeserializeOwned,
    S::Params: Clone + serde::Serialize + serde::de::DeserializeOwned,
    R: serde::Serialize + serde::de::DeserializeOwned,
{
    fn validate(&self) -> Result<()> {
        let n = self.replicas.len();
        let mut labels = self.labels.clone();
        labels.sort_unstable();
        if n == 0
            || self.parity > 1
            || self.labels.len() != n
            || self.seen_hot.len() != n
            || self.seen_cold.len() != n
            || labels != (0..n).collect::<Vec<_>>()
            || self.swap_accepts > self.swap_attempts
            || self.replicas.windows(2).any(|w| w[0].beta >= w[1].beta)
        {
            return Err(Error("invalid temperature-exchange checkpoint"));
        }
        if self.pending_steps == 0 {
            if self.local_cursor != 0 {
                return Err(Error("invalid exchange cursor"));
            }
        } else if self
            .pending_steps
            .checked_mul(n)
            .is_none_or(|total| self.local_cursor > total)
        {
            return Err(Error("invalid pending exchange chunk"));
        }
        for r in &self.replicas {
            r.validate_checkpoint()?;
        }
        Ok(())
    }
    pub fn save(
        &self,
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.validate()?;
        csta_core::checkpoint::save(path, model, self)
    }
    pub fn load(
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let s: Self = csta_core::checkpoint::load(path, model)?;
        s.validate()?;
        Ok(s)
    }
}
