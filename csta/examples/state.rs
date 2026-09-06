//! A minimal reversible model: two levels E=0 and E=1, k_b=1.
use csta::{
    Metropolis, MonteCarlo, Schedule, State, csta_derive::Randomizable, statistics::Blocking,
};
use rand::{RngExt, SeedableRng, rngs::StdRng};

#[derive(Clone, Copy, Debug, PartialEq, Randomizable)]
enum Spin {
    Down,
    Up,
}
#[derive(Clone, Debug, Randomizable)]
struct TwoLevel {
    spin: Spin,
}
impl State for TwoLevel {
    type Params = ();
    // Keep the original value so self-loops and nontrivial changes both revert exactly.
    type Change = (Spin, Spin);
    fn energy(&self, _: &mut ()) -> f64 {
        f64::from(self.spin == Spin::Up)
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        (
            self.spin,
            if rng.random_bool(0.5) {
                Spin::Up
            } else {
                Spin::Down
            },
        )
    }
    fn apply_change(&mut self, (_, next): Self::Change) {
        self.spin = next;
    }
    fn revert_change(&mut self, (old, _): Self::Change) {
        self.spin = old;
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let initial = MonteCarlo::<TwoLevel, _>::new(StdRng::seed_from_u64(42));
    for (seed, state) in initial.take(3).enumerate() {
        let mut sampler =
            Metropolis::with_all(state, (), 1.0, 21_000, StdRng::seed_from_u64(seed as u64));
        let mut energy = Blocking::new(64)?;
        sampler.try_run_observed(
            Schedule {
                burn_in: 1000,
                stride: 1,
            },
            None,
            |s, _| energy.push(s.energy(&mut ())),
        )?;
        println!(
            "seed={seed}, mean energy={:?}, exact={}",
            energy.estimate(),
            1.0 / (1.0 + 1.0_f64.exp())
        );
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn every_proposal_reverts_exactly() {
        for old in [Spin::Down, Spin::Up] {
            for next in [Spin::Down, Spin::Up] {
                let mut s = TwoLevel { spin: old };
                s.apply_change((old, next));
                assert_eq!(s.spin, next);
                s.revert_change((old, next));
                assert_eq!(s.spin, old);
            }
        }
    }
}
