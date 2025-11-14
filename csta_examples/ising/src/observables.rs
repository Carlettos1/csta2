use csta::observer::Observer;

use crate::{Ising, Spin};

pub struct Magnetization;

impl Observer<Ising> for Magnetization {
    type Observation = f64;

    fn after() -> usize {
        0 // will measure magnetization from the beggining
    }

    fn every() -> usize {
        10 // will measure every 10 steps
    }

    fn measure(state: &Ising, _params: &<Ising as csta::State>::Params) -> Self::Observation {
        state
            .states
            .iter()
            .map(|s| match s {
                Spin::Up => 1.0,
                Spin::Down => -1.0,
            })
            .sum::<f64>()
            / state.states.len() as f64
    }
}
