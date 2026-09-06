use super::*;
use rand::{SeedableRng, rngs::StdRng};
#[derive(Clone, Debug)]
struct Two {
    value: bool,
    high: f64,
}
impl State for Two {
    type Params = f64;
    type Change = bool;
    fn energy(&self, p: &mut f64) -> f64 {
        *p = if self.value { self.high } else { 0.0 };
        *p
    }
    fn propose_change(&self, _: &mut impl RngExt) -> bool {
        self.value
    }
    fn apply_change(&mut self, _: bool) {
        self.value = !self.value;
    }
    fn revert_change(&mut self, c: bool) {
        self.value = c;
    }
}
fn sampler(beta: f64) -> Metropolis<Two, StdRng> {
    Metropolis::with_all(
        Two {
            value: false,
            high: 1.0,
        },
        0.0,
        beta,
        100,
        StdRng::seed_from_u64(42),
    )
}
#[test]
fn rollback_counts_and_negative_beta() {
    let mut m = sampler(1000.0);
    m.params = 123.0;
    assert!(!m.try_step().unwrap());
    assert_eq!(m.params, 123.0);
    assert!(!m.state.value);
    m.beta = -1000.0;
    assert!(m.try_step().unwrap());
    assert!(!m.try_step().unwrap());
    assert_eq!(m.params, 1.0);
    assert_eq!(m.attempted_moves, 3);
    assert_eq!(m.accepted_rate(), 1.0 / 3.0);
    m.state.high = f64::NAN;
    assert!(m.try_step().is_err());
    assert_eq!(m.attempted_moves, 3);
    assert_eq!(m.params, 1.0);
    m.state.high = 1.0;
    m.beta = f64::NAN;
    assert!(m.try_step().is_err());
    m.beta = 0.0;
    m.attempted_moves = usize::MAX;
    assert!(m.try_step().is_err());
    assert_eq!(sampler(1.0).accepted_rate(), 0.0);
    for beta in [-2.0, 0.0, 2.0] {
        let mut m = sampler(beta);
        let mut sum = 0.0;
        for _ in 0..50_000 {
            m.try_step().unwrap();
            sum += f64::from(m.state.value);
        }
        assert!((sum / 50_000.0 - 1.0 / (1.0 + beta.exp())).abs() < 0.015);
    }
    assert_eq!(boltzmann_log_ratio(0.0, -f64::MAX, f64::MAX).unwrap(), 0.0);
}
struct Obs;
impl Observer<Two> for Obs {
    type Observation = bool;
    fn measure(s: &Two, _: &f64) -> bool {
        s.value
    }
    fn every() -> usize {
        2
    }
    fn after() -> usize {
        1
    }
}
#[test]
fn schedules_and_cancellation() {
    let mut m = sampler(0.0);
    m.steps = 6;
    let obs = m.run_with::<Obs>();
    assert_eq!(obs, vec![true, true]);
    assert_eq!(m.run_with_2::<Obs, Obs>(), (obs.clone(), obs.clone()));
    assert_eq!(
        m.run_with_3::<Obs, Obs, Obs>(),
        (obs.clone(), obs.clone(), obs.clone())
    );
    assert_eq!(
        m.run_with_4::<Obs, Obs, Obs, Obs>(),
        (obs.clone(), obs.clone(), obs.clone(), obs.clone())
    );
    assert_eq!(m.run_with_n(vec![Box::new(Obs)]), vec![obs]);
    assert_eq!(m.attempted_moves, 30);
    assert_eq!(m.accepted_rate(), 1.0);
    assert_eq!(m.rejected_rate(), 0.0);
    assert!(
        m.run_observed(
            Schedule {
                burn_in: 0,
                stride: 0
            },
            None,
            |_, _| {}
        )
        .is_err()
    );
    let flag = std::sync::atomic::AtomicBool::new(false);
    let mut count = 0;
    let report = m
        .run_observed(Schedule::default(), Some(&flag), |_, _| {
            count += 1;
            flag.store(true, std::sync::atomic::Ordering::Relaxed)
        })
        .unwrap();
    assert_eq!(report.attempted, 1);
    assert_eq!(report.status, RunStatus::Cancelled);
    assert_eq!(count, 1);
    m.steps = 0;
    assert_eq!(m.run(None).unwrap().attempted, 0);
}
#[test]
fn moments_and_correlations() {
    use statistics::*;
    let mut a = Moments::default();
    let mut b = Moments::default();
    let mut all = Moments::default();
    let mut covariance = Covariance::default();
    for i in 0..100 {
        let x = 1e10 + i as f64;
        all.push(x).unwrap();
        if i < 50 {
            a.push(x).unwrap();
        } else {
            b.push(x).unwrap();
        }
        covariance.push(x, 2.0 * x).unwrap();
    }
    a.merge(&b).unwrap();
    assert_eq!(a, all);
    assert!((a.variance().unwrap() - 841.6666666667).abs() < 1e-6);
    assert!((covariance.covariance().unwrap() - 2.0 * a.variance().unwrap()).abs() < 1e-6);
    let backup = a;
    assert!(a.push(f64::NAN).is_err());
    assert_eq!(a, backup);
    assert!(correlated_estimate(&[1.0; 100]).unwrap().is_none());
    assert!(correlated_estimate(&[]).unwrap().is_none());
    let mut rng = StdRng::seed_from_u64(57);
    for rho in [0.0, 0.8, -0.5] {
        let mut x = 0.0;
        let mut samples = Vec::new();
        let mut block = Blocking::new(128).unwrap();
        for i in 0..41_000 {
            let noise = csta_montecarlo::gaussian(&mut rng, 0.0, 1.0).unwrap();
            x = rho * x + noise;
            if i >= 1000 {
                samples.push(x);
                block.push(x).unwrap();
            }
        }
        let estimate = correlated_estimate(&samples).unwrap().unwrap();
        let expected = ((1.0 + rho) / (1.0 - rho)).max(1.0);
        assert!(
            (estimate.tau - expected).abs() < expected * 0.35,
            "{} vs {}",
            estimate.tau,
            expected
        );
        assert!(estimate.estimate.mean.abs() < 6.0 * estimate.estimate.standard_error);
        assert!(block.estimate().is_some());
    }
    assert!(Blocking::new(0).is_err());
}
#[test]
fn thermal_and_models() {
    use models::*;
    use statistics::*;
    let mut t = Thermodynamics::default();
    assert!(t.summary(1.0, 1.0, 1).is_err());
    for m in [-2.0, 2.0] {
        t.push(-1.0, m).unwrap();
    }
    let s = t.summary(0.5, 1.0, 2).unwrap();
    assert_eq!(s.heat_capacity, 0.0);
    assert_eq!(s.susceptibility, 1.0);
    assert!((s.binder.unwrap() - 2.0 / 3.0).abs() < 1e-14);
    let mut rng = StdRng::seed_from_u64(4);
    for boundary in [Boundary::Open, Boundary::Periodic] {
        let mut m = Ising2D::new(3, 3, boundary, 1.2, 0.3, vec![1; 9]).unwrap();
        let bonds = if boundary == Boundary::Open {
            12.0
        } else {
            18.0
        };
        assert!((m.energy(&mut ()) - (-1.2 * bonds - 2.7)).abs() < 1e-12);
        for _ in 0..1000 {
            let old = m.energy(&mut ());
            let c = m.propose_change(&mut rng);
            m.apply_change(c);
            assert!((m.energy(&mut ()) - m.recompute_energy()).abs() < 1e-10);
            m.revert_change(c);
            assert_eq!(m.energy(&mut ()), old);
            m.apply_change(c);
        }
        assert!(m.wolff_step(1.0, &mut rng).is_err());
    }
    let mut m = Ising2D::aligned(3).unwrap();
    assert_eq!(m.wolff_step(0.0, &mut rng).unwrap(), 1);
    let mut m = Ising2D::aligned(3).unwrap();
    assert_eq!(m.wolff_step(1000.0, &mut rng).unwrap(), 9);
    assert_eq!(m.magnetization(), -9.0);
    assert!(Ising2D::aligned(2).is_err());
    for mu in [-2.0, 0.0, 2.0] {
        let gas = LatticeGas::new(2, 2, Boundary::Open, 0.0, mu, vec![false; 4]).unwrap();
        let mut sampler = Metropolis::with_all(gas, (), 1.0, 0, StdRng::seed_from_u64(42));
        let mut count = 0.0;
        for _ in 0..30_000 {
            sampler.try_step().unwrap();
            assert!(
                (sampler.state.energy(&mut ()) - sampler.state.recompute_energy()).abs() < 1e-10
            );
            count += sampler.state.particles() as f64 / 4.0;
        }
        assert!((count / 30_000.0 - 1.0 / (1.0 + (-mu).exp())).abs() < 0.02);
    }
}
#[test]
fn replica_exchange() {
    use replica::ReplicaExchange;
    let states = vec![
        Two {
            value: false,
            high: 1.0,
        },
        Two {
            value: true,
            high: 1.0,
        },
    ];
    let mut r = ReplicaExchange::new(states.clone(), 0.0, vec![0.0, 2.0], 1).unwrap();
    let mut s = ReplicaExchange::new(states, 0.0, vec![0.0, 2.0], 1).unwrap();
    r.run(100, 3, None).unwrap();
    s.run(100, 3, None).unwrap();
    assert_eq!(r.swap_accepts, s.swap_accepts);
    assert_eq!(r.swap_attempts, 50);
    for (a, b) in r.replicas().iter().zip(s.replicas()) {
        assert_eq!(a.state.value, b.state.value);
        assert_eq!(a.attempted_moves, 300);
    }
    assert!(r.exchange(1).is_err());
    assert!(r.run(1, 0, None).is_err());
    assert!(
        ReplicaExchange::new(
            vec![
                Two {
                    value: false,
                    high: 1.0
                };
                2
            ],
            0.0,
            vec![1.0, 1.0],
            0
        )
        .is_err()
    );
}

#[derive(Clone)]
struct Asymmetric {
    i: usize,
    correction: f64,
}
const Q: [[f64; 3]; 3] = [[0.5, 0.25, 0.25], [0.5, 0.25, 0.25], [0.25, 0.5, 0.25]];
impl State for Asymmetric {
    type Params = ();
    type Change = (usize, usize);
    fn energy(&self, _: &mut ()) -> f64 {
        self.i as f64
    }
    fn propose_change(&self, r: &mut impl RngExt) -> Self::Change {
        let u = r.random::<f64>();
        let j = if u < Q[self.i][0] {
            0
        } else if u < Q[self.i][0] + Q[self.i][1] {
            1
        } else {
            2
        };
        (self.i, j)
    }
    fn apply_change(&mut self, c: Self::Change) {
        self.i = c.1;
    }
    fn revert_change(&mut self, c: Self::Change) {
        self.i = c.0;
    }
    fn log_proposal_ratio(&self, c: &Self::Change) -> f64 {
        (Q[c.1][c.0] / Q[c.0][c.1]).ln() + self.correction
    }
}
#[test]
fn hastings_detailed_balance_and_stationary_distribution() {
    for beta in [0.0, 0.7, -0.7] {
        let weights: [f64; 3] = std::array::from_fn(|i| (-beta * i as f64).exp());
        let z = weights.iter().sum::<f64>();
        for i in 0..3 {
            for j in 0..3 {
                let p = Q[i][j] * ((weights[j] * Q[j][i]) / (weights[i] * Q[i][j])).min(1.0);
                let reverse = Q[j][i] * ((weights[i] * Q[i][j]) / (weights[j] * Q[j][i])).min(1.0);
                assert!((weights[i] * p - weights[j] * reverse).abs() < 1e-14);
            }
        }
        let mut s = Metropolis::with_all(
            Asymmetric {
                i: 0,
                correction: 0.0,
            },
            (),
            beta,
            0,
            StdRng::seed_from_u64(8),
        );
        let mut counts = [0; 3];
        for _ in 0..80_000 {
            s.try_step().unwrap();
            counts[s.state.i] += 1;
        }
        for i in 0..3 {
            assert!((counts[i] as f64 / 80_000.0 - weights[i] / z).abs() < 0.015);
        }
    }
    let mut s = Metropolis::with_all(
        Asymmetric {
            i: 0,
            correction: f64::NEG_INFINITY,
        },
        (),
        0.0,
        0,
        StdRng::seed_from_u64(7),
    );
    for _ in 0..20 {
        assert!(!s.try_step().unwrap());
        assert_eq!(s.state.i, 0);
    }
    for bad in [f64::NAN, f64::INFINITY] {
        s.state.correction = bad;
        assert!(s.try_step().is_err());
        assert_eq!(s.attempted_moves, 20);
    }
}

#[test]
fn exact_lattice_sums_and_cluster_equilibrium() {
    use models::*;
    let beta = 0.3;
    let mut z = 0.0;
    let mut e = 0.0;
    let mut m2 = 0.0;
    for bits in 0..512 {
        let spins: Vec<i8> = (0..9)
            .map(|i| if bits & (1 << i) == 0 { -1 } else { 1 })
            .collect();
        // Independent positive-axis bond sum, not Lattice2D::bonds.
        let mut energy = 0.0;
        for y in 0..3 {
            for x in 0..3 {
                energy -= spins[y * 3 + x] as f64
                    * (spins[y * 3 + (x + 1) % 3] + spins[((y + 1) % 3) * 3 + x]) as f64;
            }
        }
        let model = Ising2D::new(3, 3, Boundary::Periodic, 1.0, 0.0, spins).unwrap();
        assert_eq!(model.recompute_energy(), energy);
        let w = (-beta * energy).exp();
        z += w;
        e += w * energy;
        m2 += w * model.magnetization().powi(2);
    }
    let mut model = Ising2D::aligned(3).unwrap();
    let mut rng = StdRng::seed_from_u64(77);
    let mut sampled_e = 0.0;
    let mut sampled_m2 = 0.0;
    for i in 0..61_000 {
        model.wolff_step(beta, &mut rng).unwrap();
        if i >= 1000 {
            sampled_e += model.energy(&mut ());
            sampled_m2 += model.magnetization().powi(2);
        }
    }
    assert!((sampled_e / 60_000.0 - e / z).abs() < 0.25);
    assert!((sampled_m2 / 60_000.0 - m2 / z).abs() < 0.8);
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(model.wolff_step(bad, &mut rng).is_err());
    }
    assert!(Ising2D::new(3, 3, Boundary::Periodic, 1.0, 0.0, vec![i8::MIN; 9]).is_err());
}

#[test]
fn interacting_grand_canonical_enumeration() {
    use models::*;
    let beta = 0.8;
    let mu = -0.4;
    let j = 0.7;
    let mut z = 0.0;
    let mut mean_n = 0.0;
    for bits in 0..16 {
        let occ: Vec<_> = (0..4).map(|i| bits & (1 << i) != 0).collect();
        let n = occ.iter().filter(|x| **x).count();
        let bonds = [(0, 1), (0, 2), (1, 3), (2, 3)]
            .iter()
            .filter(|(a, b)| occ[*a] && occ[*b])
            .count();
        let e = -j * bonds as f64 - mu * n as f64;
        let model = LatticeGas::new(2, 2, Boundary::Open, j, mu, occ).unwrap();
        assert!((model.energy(&mut ()) - e).abs() < 1e-14);
        let w = (-beta * e).exp();
        z += w;
        mean_n += w * n as f64;
    }
    let gas = LatticeGas::new(2, 2, Boundary::Open, j, mu, vec![false; 4]).unwrap();
    let mut m = Metropolis::with_all(gas, (), beta, 0, StdRng::seed_from_u64(11));
    let mut n = 0.0;
    for _ in 0..60_000 {
        m.try_step().unwrap();
        n += m.state.particles() as f64;
    }
    assert!((n / 60_000.0 - mean_n / z).abs() < 0.035);
}

#[test]
fn replica_canonical_marginals_and_rejection() {
    use replica::ReplicaExchange;
    let states = vec![
        Two {
            value: true,
            high: 1.0,
        },
        Two {
            value: false,
            high: 1.0,
        },
    ];
    let mut r = ReplicaExchange::new(states, 0.0, vec![0.0, 1000.0], 15).unwrap();
    assert!(!r.exchange(0).unwrap());
    assert!(r.replicas()[0].state.value);
    assert!(!r.replicas()[1].state.value);
    let mut r = ReplicaExchange::new(
        vec![
            Two {
                value: false,
                high: 1.0
            };
            2
        ],
        0.0,
        vec![0.1, 1.5],
        15,
    )
    .unwrap();
    let mut sums = [0.0; 2];
    for _ in 0..30_000 {
        r.run(1, 2, None).unwrap();
        for (sum, replica) in sums.iter_mut().zip(r.replicas()) {
            *sum += f64::from(replica.state.value);
        }
    }
    for (s, beta) in sums.into_iter().zip([0.1_f64, 1.5]) {
        assert!((s / 30_000.0 - 1.0 / (1.0 + beta.exp())).abs() < 0.02);
    }
    assert!(r.round_trips > 0);
    let flag = std::sync::atomic::AtomicBool::new(true);
    assert_eq!(r.run(1, 1, Some(&flag)).unwrap(), RunStatus::Cancelled);
    let mut single = ReplicaExchange::new(
        vec![Two {
            value: false,
            high: 1.0,
        }],
        0.0,
        vec![1.0],
        4,
    )
    .unwrap();
    single.run(10, 1, None).unwrap();
    assert_eq!(single.swap_attempts, 0);
    assert_eq!(single.round_trips, 0);
}
#[test]
fn block_uncertainty_and_insufficient_data() {
    use statistics::*;
    let mut complete = Blocking::new(8).unwrap();
    for i in 0..64 {
        complete.push(i as f64).unwrap();
    }
    let before = complete.estimate().unwrap();
    complete.push(1e6).unwrap();
    let after = complete.estimate().unwrap();
    assert_eq!(before.effective_samples, after.effective_samples);
    assert_eq!(before.standard_error, after.standard_error);
    let mut constant = vec![(0.0, 0.0); 64];
    constant.push((10.0, 10.0));
    assert!(thermal_errors(&constant, 1.0, 1.0, 1, 8).unwrap().is_none());

    assert!(thermal_errors(&[], 1.0, 1.0, 1, 8).unwrap().is_none());
    assert!(
        thermal_errors(&[(0.0, 0.0); 64], 1.0, 1.0, 1, 8)
            .unwrap()
            .is_none()
    );
    let samples: Vec<_> = (0..16)
        .flat_map(|i| vec![(i as f64, (i as f64 - 7.0)); 32])
        .collect();
    let blocked = thermal_errors(&samples, 1.0, 1.0, 1, 32).unwrap().unwrap();
    let iid = thermal_errors(&samples, 1.0, 1.0, 1, 1).unwrap().unwrap();
    assert!(blocked.energy > iid.energy * 5.0);
    assert!(blocked.binder.is_some());
    assert_eq!(blocked.used_samples, 512);
    let mut b = Blocking::new(8).unwrap();
    for i in 0..67 {
        b.push(i as f64).unwrap();
    }
    assert_eq!(b.incomplete_samples(), 3);
    assert_eq!(b.estimate().unwrap().samples, 64);
    let mut t = Thermodynamics::default();
    t.push(0.0, 0.0).unwrap();
    assert!(t.summary(0.0, 1.0, 1).unwrap().binder.is_none());
    assert!(t.summary(1.0, 1.0, 0).is_err());
    assert!(t.push(0.0, f64::MAX).is_err());
    assert_eq!(t.summary(1.0, 1.0, 1).unwrap().samples, 1);
}

#[derive(Clone)]
struct FixedRng(u64);
impl rand::TryRng for FixedRng {
    type Error = std::convert::Infallible;
    fn try_next_u64(&mut self) -> std::result::Result<u64, Self::Error> {
        Ok(self.0)
    }
    fn try_next_u32(&mut self) -> std::result::Result<u32, Self::Error> {
        Ok((self.0 >> 32) as u32)
    }
    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> std::result::Result<(), Self::Error> {
        for (i, b) in dst.iter_mut().enumerate() {
            *b = self.0.to_le_bytes()[i % 8];
        }
        Ok(())
    }
}
#[test]
fn acceptance_thresholds_and_symmetric_adapter() {
    let mut half = FixedRng(1 << 63);
    assert!(!accept_log(0.5_f64.ln(), &mut half));
    assert!(accept_log(0.5_f64.ln().next_up(), &mut half));
    assert!(!accept_log(f64::NEG_INFINITY, &mut FixedRng(0)));
    assert!(accept_log(0.0, &mut FixedRng(u64::MAX)));
    let mut a = sampler(0.7);
    let state = Hastings {
        state: Two {
            value: false,
            high: 1.0,
        },
        log_ratio: |_: &Two, _: &bool| 0.0,
    };
    let mut b = Metropolis::with_all(state, 0.0, 0.7, 100, StdRng::seed_from_u64(42));
    for _ in 0..100 {
        assert_eq!(a.try_step().unwrap(), b.try_step().unwrap());
        assert_eq!(a.state.value, b.state.state.value);
    }
    let mut m = sampler(0.0);
    m.steps = 3;
    let mut observations = Vec::new();
    m.run_observed(Schedule::default(), None, |s, _| observations.push(s.value))
        .unwrap();
    assert_eq!(observations, vec![true, false, true]);
    let mut count = 0;
    m.run_observed(
        Schedule {
            burn_in: 3,
            stride: 1,
        },
        None,
        |_, _| count += 1,
    )
    .unwrap();
    assert_eq!(count, 0);
    let empty: Vec<Vec<()>> = m.run_with_n(vec![]);
    assert!(empty.is_empty());
}
#[test]
fn iid_error_scaling_and_failure_atomicity() {
    use statistics::*;
    let mut r = StdRng::seed_from_u64(144);
    let samples: Vec<_> = (0..40_000)
        .map(|_| csta_montecarlo::gaussian(&mut r, 0.0, 1.0).unwrap())
        .collect();
    let small = correlated_estimate(&samples[..10_000]).unwrap().unwrap();
    let large = correlated_estimate(&samples).unwrap().unwrap();
    let ratio = small.estimate.standard_error / large.estimate.standard_error;
    assert!((ratio - 2.0).abs() < 0.3);
    let mut m = sampler(0.0);
    m.steps = 10;
    assert!(
        m.try_run_observed(Schedule::default(), None, |_, _| Err(Error(
            "measurement failed"
        )))
        .is_err()
    );
    assert_eq!(m.attempted_moves, 1);
    assert!(m.state.value);
    // An invalid proposed energy restores both state and a deliberately distinct cache.
    let mut m = sampler(1.0);
    m.state.high = f64::INFINITY;
    m.params = 72.0;
    assert!(m.try_step().is_err());
    assert_eq!(m.params, 72.0);
    assert!(!m.state.value);
    assert_eq!(m.attempted_moves, 0);
}
#[test]
fn cluster_bond_limits_and_grand_canonical_extremes() {
    use models::*;
    let mut all = Ising2D::new(2, 2, Boundary::Open, 1.0, 0.0, vec![1; 4]).unwrap();
    let mut none = all.clone();
    assert_eq!(all.wolff_step(0.5, &mut FixedRng(0)).unwrap(), 4);
    assert_eq!(none.wolff_step(0.5, &mut FixedRng(u64::MAX)).unwrap(), 1);
    let mut antiferro = Ising2D::new(2, 2, Boundary::Open, -1.0, 0.0, vec![1; 4]).unwrap();
    assert!(antiferro.wolff_step(0.5, &mut FixedRng(0)).is_err());
    assert_eq!(antiferro.spins(), &[1; 4]);
    for mu in [-1000.0, 1000.0] {
        let gas = LatticeGas::new(1, 1, Boundary::Open, 0.0, mu, vec![false]).unwrap();
        let mut m = Metropolis::with_all(gas, (), 1.0, 0, StdRng::seed_from_u64(42));
        for _ in 0..20 {
            m.try_step().unwrap();
            assert_eq!(m.state.particles(), usize::from(mu > 0.0));
        }
        m.beta = 0.0;
        let before = m.state.particles();
        assert!(m.try_step().unwrap());
        assert_eq!(m.state.particles(), 1 - before);
    }
    assert!(LatticeGas::new(0, 1, Boundary::Open, 0.0, 0.0, vec![]).is_err());
    assert!(LatticeGas::new(1, 1, Boundary::Open, 0.0, f64::NAN, vec![false]).is_err());
}
#[test]
fn replica_joint_distribution() {
    use replica::ReplicaExchange;
    let mut r = ReplicaExchange::new(
        vec![
            Two {
                value: false,
                high: 1.0
            };
            2
        ],
        0.0,
        vec![0.2, 1.0],
        64,
    )
    .unwrap();
    assert!(r.exchange(0).unwrap());
    let mut counts = [[0; 2]; 2];
    for _ in 0..40_000 {
        r.run(1, 3, None).unwrap();
        counts[usize::from(r.replicas()[0].state.value)]
            [usize::from(r.replicas()[1].state.value)] += 1;
    }
    let p = [1.0 / (1.0 + 0.2_f64.exp()), 1.0 / (1.0 + 1.0_f64.exp())];
    for (i, row) in counts.iter().enumerate() {
        for (j, count) in row.iter().enumerate() {
            let expected =
                (if i == 1 { p[0] } else { 1.0 - p[0] }) * (if j == 1 { p[1] } else { 1.0 - p[1] });
            assert!((*count as f64 / 40_000.0 - expected).abs() < 0.02);
        }
    }
}
