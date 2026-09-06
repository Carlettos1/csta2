//! Small reference models. Both use symmetric proposals and cached O(1) energy.
use crate::{EnergyGrid, Error, Randomizable, Result, State, WLData};
use rand::RngExt;

/// Periodic one-dimensional Ising chain with coupling J=1, N>=3.
#[derive(Clone, Debug, PartialEq)]
pub struct Ising<const N: usize> {
    spins: [i8; N],
    energy: i64,
}
impl<const N: usize> Ising<N> {
    pub fn new(spins: [i8; N]) -> Result<Self> {
        if N < 3 || spins.iter().any(|s| *s != -1 && *s != 1) {
            return Err(Error::Invalid("Ising needs N>=3 and spins +/-1"));
        }
        let energy = -(0..N)
            .map(|i| i64::from(spins[i]) * i64::from(spins[(i + 1) % N]))
            .sum::<i64>();
        Ok(Self { spins, energy })
    }
    pub fn grid() -> Result<EnergyGrid> {
        if N < 3 {
            return Err(Error::Invalid("Ising needs N>=3"));
        }
        EnergyGrid::discrete((0..=N / 2).map(|k| -(N as f64) + 4.0 * k as f64).collect())
    }
    pub fn spins(&self) -> &[i8; N] {
        &self.spins
    }
    pub fn recompute_energy(&self) -> f64 {
        -(0..N)
            .map(|i| f64::from(self.spins[i]) * f64::from(self.spins[(i + 1) % N]))
            .sum::<f64>()
    }
    /// Enumerate small systems independently of the DOS combinatorial formula.
    pub fn exact_dos() -> Result<WLData> {
        if !(3..=20).contains(&N) {
            return Err(Error::Invalid("exact enumeration supports 3<=N<=20"));
        }
        let grid = Self::grid()?;
        let mut counts = vec![0_u64; grid.len()];
        for bits in 0..(1_u64 << N) {
            let spins = std::array::from_fn(|i| if bits & (1 << i) == 0 { -1 } else { 1 });
            let state = Self::new(spins)?;
            counts[grid.bin(state.recompute_energy())?.unwrap()] += 1;
        }
        WLData::new(
            grid,
            counts.iter().map(|c| (*c as f64).ln()).collect(),
            counts,
        )
    }
}
impl<const N: usize> Randomizable for Ising<N> {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Self::new(std::array::from_fn(|_| {
            if rng.random_bool(0.5) { 1 } else { -1 }
        }))
        .expect("Ising N>=3")
    }
}
impl<const N: usize> State for Ising<N> {
    type Params = ();
    type Change = (usize, i64);
    fn energy(&self, _: &mut ()) -> f64 {
        self.energy as f64
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        let i = rng.random_range(0..N);
        let delta = 2
            * i64::from(self.spins[i])
            * (i64::from(self.spins[(i + N - 1) % N]) + i64::from(self.spins[(i + 1) % N]));
        (i, delta)
    }
    fn apply_change(&mut self, (i, delta): Self::Change) {
        self.spins[i] *= -1;
        self.energy += delta;
    }
    fn revert_change(&mut self, (i, delta): Self::Change) {
        self.spins[i] *= -1;
        self.energy -= delta;
    }
}

/// N distinguishable quantum oscillators with hbar*omega=1:
/// E=sum(n_i+1/2), n_i>=0. A downward proposal at zero is a self-loop,
/// preserving symmetry of every nontrivial proposal (unlike forced upward moves).
#[derive(Clone, Debug, PartialEq)]
pub struct Oscillators<const N: usize> {
    occupation: [u32; N],
    quanta: u64,
}
impl<const N: usize> Oscillators<N> {
    pub fn new(occupation: [u32; N]) -> Result<Self> {
        if N == 0 {
            return Err(Error::Invalid("need at least one oscillator"));
        }
        let quanta = occupation
            .iter()
            .try_fold(0_u64, |sum, n| sum.checked_add(u64::from(*n)))
            .ok_or(Error::CounterOverflow)?;
        if quanta > (1_u64 << 52) {
            return Err(Error::Invalid("occupation energy exceeds exact f64 range"));
        }
        Ok(Self { occupation, quanta })
    }
    pub fn occupation(&self) -> &[u32; N] {
        &self.occupation
    }
    pub fn grid(max_quanta: u32) -> Result<EnergyGrid> {
        if N == 0 {
            return Err(Error::Invalid("need at least one oscillator"));
        }
        EnergyGrid::discrete(
            (0..=max_quanta)
                .map(|q| f64::from(q) + N as f64 * 0.5)
                .collect(),
        )
    }
    /// Exact bin masses for the truncated total-quanta window, C(q+N-1,N-1).
    pub fn exact_dos(max_quanta: u32) -> Result<WLData> {
        let grid = Self::grid(max_quanta)?;
        let mut log_g = 0.0;
        let dos = (0..=max_quanta)
            .map(|q| {
                if q > 0 {
                    log_g += (f64::from(q) + N as f64 - 1.0).ln() - f64::from(q).ln();
                }
                log_g
            })
            .collect();
        let n = grid.len();
        WLData::new(grid, dos, vec![0; n])
    }
}
impl<const N: usize> Randomizable for Oscillators<N> {
    fn sample<R: RngExt + ?Sized>(_rng: &mut R) -> Self {
        Self::new([0; N]).expect("oscillator N>0")
    }
}
impl<const N: usize> State for Oscillators<N> {
    type Params = ();
    type Change = (usize, u32, u32);
    fn energy(&self, _: &mut ()) -> f64 {
        self.quanta as f64 + N as f64 * 0.5
    }
    fn propose_change(&self, rng: &mut impl RngExt) -> Self::Change {
        let i = rng.random_range(0..N);
        let old = self.occupation[i];
        let next = if rng.random_bool(0.5) {
            old.saturating_add(1)
        } else {
            old.saturating_sub(1)
        };
        (i, old, next)
    }
    fn apply_change(&mut self, (i, old, next): Self::Change) {
        self.quanta = self.quanta - u64::from(old) + u64::from(next);
        self.occupation[i] = next;
    }
    fn revert_change(&mut self, (i, old, next): Self::Change) {
        self.quanta = self.quanta - u64::from(next) + u64::from(old);
        self.occupation[i] = old;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};
    #[test]
    fn ising_enumeration_and_cached_energy() {
        let exact = Ising::<4>::exact_dos().unwrap();
        assert_eq!(exact.bins(), &[2, 12, 2]);
        assert_eq!(exact.grid().energies(), &[-4.0, 0.0, 4.0]);
        assert_eq!(exact.bins().iter().sum::<u64>(), 16);
        let mut state = Ising::<7>::new([1; 7]).unwrap();
        let mut rng = StdRng::seed_from_u64(41);
        for _ in 0..1000 {
            let before = state.clone();
            let change = state.propose_change(&mut rng);
            state.apply_change(change);
            assert_eq!(state.energy(&mut ()), state.recompute_energy());
            state.revert_change(change);
            assert_eq!(state, before);
            state.apply_change(change);
        }
        assert!(Ising::<2>::new([1, -1]).is_err());
        assert!(Ising::<3>::new([0, 1, 1]).is_err());
    }
    #[test]
    fn exact_ising_thermodynamics_match_state_sums() {
        let exact = Ising::<4>::exact_dos().unwrap();
        for beta in [-1.0, 0.0, 0.5, 2.0] {
            let mut z = 0.0;
            let mut e = 0.0;
            let mut e2 = 0.0;
            for bits in 0..16 {
                let state = Ising::<4>::new(std::array::from_fn(|i| {
                    if bits & (1 << i) == 0 { -1 } else { 1 }
                }))
                .unwrap();
                let energy = state.recompute_energy();
                let w = (-beta * energy).exp();
                z += w;
                e += energy * w;
                e2 += energy * energy * w;
            }
            let (a, b) = exact.energy_moments(beta).unwrap();
            assert!((a - e / z).abs() < 1e-12);
            assert!((b - e2 / z).abs() < 1e-12);
            assert!((exact.log_partition_function(beta).unwrap() - z.ln()).abs() < 1e-12);
        }
    }
    #[test]
    fn oscillator_support_symmetry_and_reversion() {
        let mut rng = StdRng::seed_from_u64(71);
        for initial in [[0, 0], [2, 3], [u32::MAX, 0]] {
            let mut state = Oscillators::<2>::new(initial).unwrap();
            for _ in 0..100 {
                let before = state.clone();
                let c = state.propose_change(&mut rng);
                state.apply_change(c);
                assert_eq!(
                    state.energy(&mut ()),
                    state
                        .occupation
                        .iter()
                        .map(|n| f64::from(*n) + 0.5)
                        .sum::<f64>()
                );
                state.revert_change(c);
                assert_eq!(state, before);
            }
        }
        let ground = Oscillators::<1>::new([0]).unwrap();
        let mut self_loops = 0;
        for _ in 0..100 {
            if ground.propose_change(&mut rng).2 == 0 {
                self_loops += 1;
            }
        }
        assert!(self_loops > 0 && self_loops < 100);
        let exact = Oscillators::<2>::exact_dos(3).unwrap();
        for (i, g) in exact.dos().iter().enumerate() {
            assert!((g.exp() - (i + 1) as f64).abs() < 1e-12);
        }
        assert!(Oscillators::<0>::new([]).is_err());
    }
}
