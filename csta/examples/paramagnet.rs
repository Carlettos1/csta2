//! Independent spins in a field: magnetization and the Schottky heat-capacity peak.
//! Run `cargo run -p csta --release --example paramagnet`. See EXAMPLES.md, problem 1.
use csta::{
    Metropolis, Schedule, State,
    statistics::{Thermodynamics, thermal_errors},
};
use rand::{RngExt, SeedableRng, rngs::StdRng};

const SITES: usize = 32;
const FIELD: f64 = 1.0;
const BURN_IN: usize = 20_000;
const SAMPLES: usize = 256_000;
const BLOCK: usize = 1024;

#[derive(Clone, Debug, PartialEq)]
struct Paramagnet {
    spins: [i8; SITES],
    magnetization: i32,
}
impl Paramagnet {
    fn aligned() -> Self {
        Self {
            spins: [1; SITES],
            magnetization: SITES as i32,
        }
    }
}
impl State for Paramagnet {
    type Params = f64; // field h; H = -h M, with k_B = 1
    type Change = (usize, i8, i8); // site, old spin, proposed spin
    fn energy(&self, field: &mut f64) -> f64 {
        -*field * f64::from(self.magnetization)
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        let site = rng.random_range(0..SITES);
        // Symmetric replacement includes self-loops, also at beta = 0.
        let next = if rng.random_bool(0.5) { 1 } else { -1 };
        (site, self.spins[site], next)
    }
    fn apply_change(&mut self, (site, old, next): Self::Change) {
        self.spins[site] = next;
        self.magnetization += i32::from(next - old);
    }
    fn revert_change(&mut self, (site, old, next): Self::Change) {
        self.spins[site] = old;
        self.magnetization += i32::from(old - next);
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "beta,energy_per_spin,energy_se,energy_exact,magnetization_per_spin,magnetization_exact,heat_capacity_per_spin,heat_capacity_se,heat_capacity_exact"
    );
    for (index, beta) in [0.0, 0.25, 0.5, 1.0, 2.0, 4.0].into_iter().enumerate() {
        let mut sampler = Metropolis::with_all(
            Paramagnet::aligned(),
            FIELD,
            beta,
            BURN_IN + SAMPLES,
            StdRng::seed_from_u64(42 + index as u64),
        );
        let mut moments = Thermodynamics::default();
        let mut pairs = Vec::with_capacity(SAMPLES);
        sampler.try_run_observed(
            Schedule {
                burn_in: BURN_IN,
                stride: 1,
            },
            None,
            |state, field| {
                let magnetization = f64::from(state.magnetization);
                let energy = -field * magnetization;
                moments.push(energy, magnetization)?;
                pairs.push((energy, magnetization));
                Ok(())
            },
        )?;
        let summary = moments.summary(beta, 1.0, SITES)?;
        let error = thermal_errors(&pairs, beta, 1.0, SITES, BLOCK)?;
        let se = |value: Option<f64>| value.map_or_else(|| "NA".into(), |v| format!("{v:.8}"));
        let x = beta * FIELD;
        let exact_m = x.tanh();
        let exact_c = x * x / x.cosh().powi(2);
        println!(
            "{beta},{:.8},{},{:.8},{:.8},{exact_m:.8},{:.8},{},{exact_c:.8}",
            summary.energy,
            se(error.map(|e| e.energy)),
            -FIELD * exact_m,
            summary.magnetization,
            summary.heat_capacity,
            se(error.map(|e| e.heat_capacity))
        );
    }
    eprintln!(
        "k_B=1, h={FIELD}, spins={SITES}, samples={SAMPLES}, block={BLOCK}; SEs are block-jackknife estimates. Vary blocks and seeds before drawing precision conclusions."
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn replacements_and_self_loops_restore_spins_and_cached_magnetization() {
        let mut state = Paramagnet::aligned();
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..500 {
            let before = state.clone();
            let change = state.propose_change(&mut rng);
            state.apply_change(change);
            let sum: i32 = state.spins.iter().map(|s| i32::from(*s)).sum();
            assert_eq!(state.magnetization, sum);
            let mut field = FIELD;
            assert_eq!(state.energy(&mut field), -FIELD * f64::from(sum));
            state.revert_change(change);
            assert_eq!(state, before);
            state.apply_change(change);
        }
    }
}
