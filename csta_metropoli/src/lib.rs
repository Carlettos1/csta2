//! This module is for metropoli + montecarlo simulations

use crate::observer::*;
use rand::{Rng, rngs::ThreadRng};

pub mod observer;

pub trait State {
    type Params;
    type Change: Clone;
    // todo: consider modification error
    // maybe is better to propose_change to be lightweight, and offload most
    // of the computation to apply or revert change. Returning an error can be useful for
    // solving the issue (thinking of panics inside apply_change and revert_change)
    // type ModificationError

    fn energy(&self, params: &mut Self::Params) -> f64;
    fn propose_change(&self, rng: &mut impl Rng) -> Self::Change;
    fn apply_change(&mut self, change: Self::Change) /* -> ModificationError */;
    fn revert_change(&mut self, change: Self::Change) /* -> ModificationError */;
}

pub struct Metropolis<S: State, R: Rng> {
    pub state: S,
    pub params: S::Params,
    pub beta: f64,
    pub steps: usize,
    pub accepted_moves: usize,
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
    R: Rng,
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
    R: Rng,
{
    pub fn with_rng(beta: f64, steps: usize, rng: R) -> Self {
        Self::with_all(S::default(), S::Params::default(), beta, steps, rng)
    }

    pub fn with_rng_no_beta(steps: usize, rng: R) -> Self {
        Self::with_rng(1.0, steps, rng)
    }
}

impl<S: State, R: Rng> Metropolis<S, R> {
    pub fn with_all(state: S, params: S::Params, beta: f64, steps: usize, rng: R) -> Self {
        Self {
            state,
            params,
            beta,
            steps,
            accepted_moves: 0,
            rng,
        }
    }

    pub fn with_all_no_beta(state: S, params: S::Params, steps: usize, rng: R) -> Self {
        Self::with_all(state, params, 1.0, steps, rng)
    }

    /// Metropolis algorithm
    pub fn step(&mut self) {
        let change = self.state.propose_change(&mut self.rng);
        let old_energy = self.state.energy(&mut self.params);
        self.state.apply_change(change.clone());
        let new_energy = self.state.energy(&mut self.params);
        let delta_energy = new_energy - old_energy;
        if delta_energy < 0.0 || self.rng.random::<f64>() < (-self.beta * delta_energy).exp() {
            self.accepted_moves += 1;
        } else {
            self.state.revert_change(change);
        }
    }

    pub fn run_empty(&mut self) {
        for _ in 0..self.steps {
            self.step();
        }
    }

    pub fn run_with<O: Observer<S>>(&mut self) -> Vec<O::Observation> {
        let mut measures: Vec<O::Observation> = Vec::new();
        for i in 0..self.steps {
            if i > O::after() && i % O::every() == 0 {
                measures.push(O::measure(&self.state, &self.params));
            }
            self.step();
        }
        measures
    }

    pub fn run_with_2<O1, O2>(&mut self) -> (Vec<O1::Observation>, Vec<O2::Observation>)
    where
        O1: Observer<S>,
        O2: Observer<S>,
    {
        let mut o1_measures: Vec<O1::Observation> = Vec::new();
        let mut o2_measures: Vec<O2::Observation> = Vec::new();

        for i in 0..self.steps {
            if i > O1::after() && i % O1::every() == 0 {
                o1_measures.push(O1::measure(&self.state, &self.params));
            }
            if i > O2::after() && i % O2::every() == 0 {
                o2_measures.push(O2::measure(&self.state, &self.params));
            }
            self.step();
        }

        (o1_measures, o2_measures)
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
        let mut o1_measures: Vec<O1::Observation> = Vec::new();
        let mut o2_measures: Vec<O2::Observation> = Vec::new();
        let mut o3_measures: Vec<O3::Observation> = Vec::new();

        for i in 0..self.steps {
            if i > O1::after() && i % O1::every() == 0 {
                o1_measures.push(O1::measure(&self.state, &self.params));
            }
            if i > O2::after() && i % O2::every() == 0 {
                o2_measures.push(O2::measure(&self.state, &self.params));
            }
            if i > O3::after() && i % O3::every() == 0 {
                o3_measures.push(O3::measure(&self.state, &self.params));
            }
            self.step();
        }

        (o1_measures, o2_measures, o3_measures)
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
        let mut o1_measures: Vec<O1::Observation> = Vec::new();
        let mut o2_measures: Vec<O2::Observation> = Vec::new();
        let mut o3_measures: Vec<O3::Observation> = Vec::new();
        let mut o4_measures: Vec<O4::Observation> = Vec::new();

        for i in 0..self.steps {
            if i > O1::after() && i % O1::every() == 0 {
                o1_measures.push(O1::measure(&self.state, &self.params));
            }
            if i > O2::after() && i % O2::every() == 0 {
                o2_measures.push(O2::measure(&self.state, &self.params));
            }
            if i > O3::after() && i % O3::every() == 0 {
                o3_measures.push(O3::measure(&self.state, &self.params));
            }
            if i > O4::after() && i % O4::every() == 0 {
                o4_measures.push(O4::measure(&self.state, &self.params));
            }
            self.step();
        }

        (o1_measures, o2_measures, o3_measures, o4_measures)
    }

    /// Makes a run for a vec of observers
    /// All observers needs to have the same Observation type
    /// And need to be dyn-capable
    /// All Observers are made dyn using DynObserver
    pub fn run_with_n<Obs>(
        &mut self,
        obs: Vec<Box<dyn DynObserver<S, Observation = Obs>>>,
    ) -> Vec<Vec<Obs>> {
        let mut measures: Vec<Vec<Obs>> = Vec::new();
        for _ in obs.iter() {
            measures.push(Vec::new());
        }

        for i in 0..self.steps {
            for (j, o) in obs.iter().enumerate() {
                if i > o.after() && i % o.every() == 0 {
                    measures[j].push(o.measure(&self.state, &self.params));
                }
            }
            self.step();
        }
        measures
    }

    pub fn accepted_rate(&self) -> f64 {
        self.accepted_moves as f64 / self.steps as f64
    }

    pub fn rejected_rate(&self) -> f64 {
        1.0 - self.accepted_rate()
    }
}
