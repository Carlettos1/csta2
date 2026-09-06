//! Biased independence proposals corrected to canonical two-state probabilities.
use csta::{Hastings, Metropolis, Schedule, State};
use rand::{RngExt, SeedableRng, rngs::StdRng};
struct Bit(bool);
impl State for Bit {
    type Params = ();
    type Change = (bool, bool);
    fn energy(&self, _: &mut ()) -> f64 {
        f64::from(self.0)
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        (self.0, rng.random_bool(0.75))
    }
    fn apply_change(&mut self, (_, next): Self::Change) {
        self.0 = next;
    }
    fn revert_change(&mut self, (old, _): Self::Change) {
        self.0 = old;
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Hastings {
        state: Bit(false),
        log_ratio: |_: &Bit, &(old, next): &(bool, bool)| {
            let probability = |b| if b { 0.75_f64 } else { 0.25 };
            (probability(old) / probability(next)).ln()
        },
    };
    let mut sampler = Metropolis::with_all(model, (), 1.0, 21_000, StdRng::seed_from_u64(42));
    let mut occupied = 0;
    sampler.run_observed(
        Schedule {
            burn_in: 1000,
            stride: 1,
        },
        None,
        |s, _| occupied += usize::from(s.state.0),
    )?;
    println!(
        "occupied fraction={}, exact={}",
        occupied as f64 / 20_000.0,
        1.0 / (1.0 + 1.0_f64.exp())
    );
    Ok(())
}
