use csta::{csta_derive::Randomizable, prelude::*};
use rand::rngs::ThreadRng;

fn main() {
    MonteCarlo::<Something, _>::default()
        .take(10)
        .for_each(|s| {
            let mut metropoli = Metropolis::with_state(s, 1.5, 1000);
            metropoli.run_empty();
        });
}

#[derive(Debug, Randomizable)]
struct BoolState(#[csta(default)] bool);

#[derive(Debug, Randomizable)]
struct Something {
    #[csta(after(BoolState(rng.random_bool(0.5))))]
    state: BoolState,
}

impl State for Something {
    type Params = ();
    type Change = bool;

    fn energy(&self, _params: &mut Self::Params) -> f64 {
        if self.state.0 { -1.0 } else { 1.0 }
    }

    fn propose_change(&self, rng: &mut impl rand::Rng) -> Self::Change {
        rng.random_bool(0.5)
    }

    fn apply_change(&mut self, change: Self::Change) {
        self.state.0 = change;
    }

    fn revert_change(&mut self, change: Self::Change) {
        self.state.0 = !change;
    }
}

struct Algo {
    metropolis: Metropolis<Something, ThreadRng>,
}

impl State for Algo {
    type Change = f64;
    type Params = ();

    fn propose_change(&self, rng: &mut impl rand::Rng) -> Self::Change {
        rng.random::<f64>()
    }

    fn apply_change(&mut self, change: Self::Change) {
        self.metropolis.beta += change;
    }

    fn revert_change(&mut self, change: Self::Change) {
        self.metropolis.beta -= change;
    }

    fn energy(&self, params: &mut Self::Params) -> f64 {
        self.metropolis.state.energy(params)
    }
}

impl Randomizable for Algo {
    fn sample<R: rand::Rng + ?Sized>(_rng: &mut R) -> Self {
        let mut rng = rand::rng();
        Algo {
            metropolis: Metropolis::with_all(Something::sample(&mut rng), (), 1.5, 1_000, rng),
        }
    }
}

fn using_algo() {
    MonteCarlo::<Algo, _>::default().take(10).for_each(|algo| {
        let beta = algo.metropolis.beta;
        let mut metropolis = Metropolis::with_state(algo, beta, 100);
        metropolis.run_empty();
    });
}
