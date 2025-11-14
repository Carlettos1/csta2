use csta::{Metropolis, MonteCarlo, State, csta_derive::Randomizable};

use crate::observables::Magnetization;

mod observables;

fn main() {
    // init montecarlo
    let mc = MonteCarlo::<Ising, _>::default();

    // make 10 states
    mc.take(10).enumerate().for_each(|(i, ising)| {
        // init metropolis
        let mut metropolis = Metropolis::with_state(ising, (i as f64 + 1.0) / 5.0, 2_000);

        // running 1 observer
        let magnetizations = metropolis.run_with::<Magnetization>();

        // running 2 observers
        metropolis.run_with_2::<Magnetization, Magnetization>();

        // running n observers
        metropolis.run_with_n(vec![
            Box::new(Magnetization),
            Box::new(Magnetization),
            Box::new(Magnetization),
            Box::new(Magnetization),
            Box::new(Magnetization),
            Box::new(Magnetization),
        ]);

        // show first 1 observer results
        println!("{:?}", magnetizations);
    });
}

#[derive(Randomizable, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Debug, Clone)]
pub struct IsingParams {
    j: f64,
}

impl Default for IsingParams {
    fn default() -> Self {
        IsingParams { j: 2.0 }
    }
}

#[derive(Randomizable, Debug)]
pub struct Ising {
    #[csta(default = 10)]
    w: usize,
    #[csta(range(5..10))]
    h: usize,
    #[csta(len(w * h))]
    pub states: Vec<Spin>,
}

impl Spin {
    fn flip(&mut self) {
        match self {
            Spin::Down => *self = Spin::Up,
            Spin::Up => *self = Spin::Down,
        }
    }

    fn mul(&self, other: &Self) -> f64 {
        if self == other { 1.0 } else { -1.0 }
    }
}

impl State for Ising {
    type Change = usize;
    type Params = IsingParams;

    fn propose_change(&self, rng: &mut impl rand::Rng) -> Self::Change {
        rng.random_range(0..self.w * self.h)
    }

    fn apply_change(&mut self, change: Self::Change) {
        self.states[change].flip();
    }

    fn revert_change(&mut self, change: Self::Change) {
        self.states[change].flip();
    }

    fn energy(&self, params: &mut Self::Params) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.w * self.h {
            for j in [i + 1, i - 1, i + self.w, i - self.w] {
                if let Some(other) = self.states.get(j) {
                    energy -= params.j * self.states[i].mul(other);
                }
            }
        }
        energy
    }
}
