use crate::{analysis::*, joint::*, *};
#[test]
fn conditional_and_histogram_reweighting() {
    let grid = EnergyGrid::discrete(vec![0.0, 1.0]).unwrap();
    let dos = WLData::new(grid.clone(), vec![2_f64.ln(), 0.0], vec![2, 1]).unwrap();
    let mut c = ConditionalMoments::new(grid);
    assert!(c.evaluate(&dos, 0.0).is_err());
    c.record(0.0, -1.0).unwrap();
    c.record(0.0, 1.0).unwrap();
    c.record(1.0, 2.0).unwrap();
    let result = c.evaluate(&dos, 0.0).unwrap();
    assert!((result.mean - 2.0 / 3.0).abs() < 1e-14);
    assert!((result.second_moment - 2.0).abs() < 1e-14);
    assert!(c.record(1.0, f64::NAN).is_err());
    assert!(c.record(2.0, 0.0).is_err());
    let samples = vec![(0.0, -1.0), (0.0, 1.0), (1.0, 2.0)];
    assert!((reweight(&samples, 0.0, 0.0).unwrap().mean - result.mean).abs() < 1e-14);
    let b = 2.0_f64;
    let expected = 2.0 * (-b).exp() / (2.0 + (-b).exp());
    assert!((reweight(&samples, 0.0, b).unwrap().mean - expected).abs() < 1e-14);
    assert_eq!(
        reweight(&samples, 0.0, 1000.0).unwrap().overlap,
        Overlap::Low
    );
    assert!(reweight(&[], 0.0, 0.0).is_err());
}
#[test]
fn ensemble_offsets_and_uncertainty() {
    let data = |offset: f64| {
        WLData::new(
            EnergyGrid::discrete(vec![0.0, 1.0]).unwrap(),
            vec![offset, offset],
            vec![10; 2],
        )
        .unwrap()
    };
    let ensemble = DosEnsemble::new(
        vec![
            IndependentDos::from_data(1, data(0.0)).unwrap(),
            IndependentDos::from_data(2, data(1000.0)).unwrap(),
        ],
        2_f64.ln(),
    )
    .unwrap();
    let x = ensemble.energy(0.0).unwrap();
    assert!((x.mean - 0.5).abs() < 1e-12);
    assert!(x.standard_error.unwrap() < 1e-12);
    assert!(
        ensemble
            .log_dos_uncertainty()
            .unwrap()
            .iter()
            .all(Option::is_some)
    );
    assert!(
        DosEnsemble::new(
            vec![
                IndependentDos::from_data(1, data(0.0)).unwrap(),
                IndependentDos::from_data(1, data(0.0)).unwrap()
            ],
            0.0
        )
        .is_err()
    );
    let single =
        DosEnsemble::new(vec![IndependentDos::from_data(0, data(0.0)).unwrap()], 0.0).unwrap();
    assert!(single.energy(0.0).unwrap().standard_error.is_none());
}
#[test]
fn joint_exact_field_and_marginal() {
    // Two independent spins, E0=0; degeneracies M=-2,0,+2 are 1,2,1.
    let grid = JointGrid::new(vec![(0.0, -2.0), (0.0, 0.0), (0.0, 2.0)], 3).unwrap();
    let mut d = JointDos::new(grid, vec![0.0, 2_f64.ln(), 0.0]).unwrap();
    d.normalize_log_count(4_f64.ln()).unwrap();
    for h in [-1.0, 0.0, 1.0] {
        let r = d.evaluate(0.7, h).unwrap();
        assert!((r.magnetization - 2.0 * (0.7 * h).tanh()).abs() < 1e-13);
        assert!((r.energy + h * r.magnetization).abs() < 1e-13);
        assert!((r.probabilities.iter().sum::<f64>() - 1.0).abs() < 1e-13);
    }
    assert_eq!(
        d.energy_marginal().unwrap().energy_moments(1.0).unwrap(),
        (0.0, 0.0)
    );
    assert!(JointGrid::new(vec![(0.0, 0.0); 2], 2).is_err());
    assert!(JointGrid::new(vec![(0.0, 0.0)], 0).is_err());
    assert!(d.evaluate(f64::NAN, 0.0).is_err());
}

#[derive(Clone)]
struct JointSpin {
    spins: [i8; 2],
}
impl State for JointSpin {
    type Params = ();
    type Change = usize;
    fn energy(&self, _: &mut ()) -> f64 {
        0.0
    }
    fn propose_change(&self, r: &mut impl rand::RngExt) -> usize {
        r.random_range(0..2)
    }
    fn apply_change(&mut self, i: usize) {
        self.spins[i] *= -1;
    }
    fn revert_change(&mut self, i: usize) {
        self.spins[i] *= -1;
    }
}
#[test]
fn joint_sampling_and_partial_results() {
    use rand::{SeedableRng, rngs::StdRng};
    let grid = JointGrid::new(vec![(0.0, -2.0), (0.0, 0.0), (0.0, 2.0)], 3).unwrap();
    let config = Config {
        preliminary_stages: 2,
        min_stage_steps: 100,
        visits_per_bin: 20,
        sampling_steps: 30_000,
        max_steps: 100_000,
        ..Default::default()
    };
    let result = run_joint(
        JointSpin { spins: [1, 1] },
        (),
        StdRng::seed_from_u64(3),
        grid.clone(),
        |s| s.spins.iter().map(|x| *x as f64).sum(),
        config.clone(),
        None,
    )
    .unwrap();
    assert!(result.is_complete());
    let x = result.data.evaluate(1.0, 0.4).unwrap();
    assert!((x.magnetization - 2.0 * 0.4_f64.tanh()).abs() < 0.15);
    let cancel = std::sync::atomic::AtomicBool::new(true);
    let result = run_joint(
        JointSpin { spins: [1, 1] },
        (),
        StdRng::seed_from_u64(3),
        grid,
        |s| s.spins.iter().map(|x| *x as f64).sum(),
        config,
        Some(&cancel),
    )
    .unwrap();
    assert!(!result.is_complete());
    assert_eq!(result.diagnostics.proposals, 0);
}

#[derive(Clone)]
struct Asymmetric {
    i: usize,
}
impl State for Asymmetric {
    type Params = ();
    type Change = (usize, usize);
    fn energy(&self, _: &mut ()) -> f64 {
        self.i as f64
    }
    fn propose_change(&self, r: &mut impl rand::RngExt) -> Self::Change {
        let u = r.random::<f64>();
        (
            self.i,
            if u < 0.25 {
                0
            } else if u < 0.75 {
                1
            } else {
                2
            },
        )
    }
    fn apply_change(&mut self, c: Self::Change) {
        self.i = c.1;
    }
    fn revert_change(&mut self, c: Self::Change) {
        self.i = c.0;
    }
    fn log_proposal_ratio(&self, c: &Self::Change) -> f64 {
        let q: [f64; 3] = [0.25, 0.5, 0.25];
        (q[c.0] / q[c.1]).ln()
    }
}
#[test]
fn wl_hastings_frozen_dos() {
    use rand::{SeedableRng, rngs::StdRng};
    let data =
        RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0, 1.0, 2.0]).unwrap()).unwrap();
    let mut w = Walker::new(
        Asymmetric { i: 0 },
        (),
        StdRng::seed_from_u64(7),
        data,
        Config::default(),
    )
    .unwrap();
    let mut counts = [0; 3];
    for _ in 0..60_000 {
        w.transition(None).unwrap();
        counts[w.bin] += 1;
    }
    for n in counts {
        assert!((n as f64 / 60_000.0 - 1.0 / 3.0).abs() < 0.015);
    }
    assert_eq!(w.data.total_visits(), 0);
}
#[test]
fn continuous_mass_and_quadratic_microcanonical_stencil() {
    let grid = EnergyGrid::continuous(0.0, 4.0, 2).unwrap();
    let d = WLData::new(grid, vec![0.0, 0.0], vec![1, 1]).unwrap();
    assert!((d.log_partition_function(0.0).unwrap() - 2_f64.ln()).abs() < 1e-14);
    assert_eq!(d.energy_moments(0.0).unwrap(), (2.0, 5.0));
    for energies in [vec![0.0, 1.0, 2.0], vec![0.0, 1.0, 3.0]] {
        let dos = energies.iter().map(|e| e * e).collect();
        let d = WLData::new(EnergyGrid::discrete(energies).unwrap(), dos, vec![1; 3]).unwrap();
        let temperatures = d.microcanonical_temperature(1.0).unwrap();
        assert!((temperatures[1].unwrap() - 0.5).abs() < 1e-14);
    }
}
#[test]
fn analysis_invalid_support_and_uncertainty_propagation() {
    let grid = EnergyGrid::discrete(vec![0.0, 1.0]).unwrap();
    let data = |g: Vec<f64>| WLData::new(grid.clone(), g, vec![1; 2]).unwrap();
    let ensemble = DosEnsemble::new(
        vec![
            IndependentDos::from_data(0, data(vec![0.0, 0.0])).unwrap(),
            IndependentDos::from_data(1, data(vec![0.0, 3_f64.ln()])).unwrap(),
        ],
        0.0,
    )
    .unwrap();
    let result = ensemble.energy(0.0).unwrap();
    assert!((result.mean - 0.625).abs() < 1e-14);
    assert!((result.standard_error.unwrap() - 0.125).abs() < 1e-14);
    assert!(
        DosEnsemble::new(
            vec![
                IndependentDos::from_data(0, data(vec![0.0, 0.0])).unwrap(),
                IndependentDos::from_data(1, data(vec![0.0, f64::NEG_INFINITY])).unwrap()
            ],
            0.0
        )
        .is_err()
    );
    let mut conditional = ConditionalMoments::new(grid);
    conditional.record(0.0, 1.0).unwrap();
    let d = WLData::new(
        EnergyGrid::discrete(vec![0.0, 2.0]).unwrap(),
        vec![0.0; 2],
        vec![1; 2],
    )
    .unwrap();
    assert!(conditional.evaluate(&d, 1.0).is_err());
    assert!(reweight(&[(f64::MAX, 0.0), (-f64::MAX, 0.0)], 1.0, 2.0).is_err());
    let flag = std::sync::atomic::AtomicBool::new(true);
    let run = run(
        models::Ising::<4>::new([1; 4]).unwrap(),
        (),
        rand::rng(),
        RawWangLandauData::on_grid(models::Ising::<4>::grid().unwrap()).unwrap(),
        Config::default(),
        Some(&flag),
    )
    .unwrap();
    assert!(IndependentDos::from_run(9, run).is_err());
}
#[test]
fn session_observes_only_retained_production_states() {
    use rand::{SeedableRng, rngs::StdRng};
    let grid = models::Ising::<4>::grid().unwrap();
    let mut c = ConditionalMoments::new(grid.clone());
    let mut s = Session::new(
        models::Ising::<4>::new([1; 4]).unwrap(),
        (),
        StdRng::seed_from_u64(42),
        RawWangLandauData::on_grid(grid).unwrap(),
        Config {
            preliminary_stages: 1,
            min_stage_steps: 10,
            visits_per_bin: 3,
            sampling_steps: 2000,
            ..Default::default()
        },
    )
    .unwrap();
    let mut observations = 0;
    s.advance_observed(100_000, None, |state, energy| {
        observations += 1;
        c.record(energy, state.spins().iter().map(|x| *x as f64).sum())
    })
    .unwrap();
    assert_eq!(observations, 2000);
    assert_eq!(s.diagnostics().sampling_steps, 2000);
    let result = s.finish();
    assert!(result.is_complete());
    assert!(
        c.evaluate(&result.data.process_data().unwrap(), 0.3)
            .is_ok()
    );
}
#[test]
fn joint_offsets_and_extreme_fields() {
    let grid = JointGrid::new(vec![(0.0, -1.0), (0.0, 1.0)], 2).unwrap();
    let a = JointDos::new(grid.clone(), vec![0.0; 2]).unwrap();
    let b = JointDos::new(grid, vec![10000.0; 2]).unwrap();
    for beta in [-1e200, -0.4, 0.0, 0.4, 1e200] {
        let x = a.evaluate(beta, 1.0).unwrap();
        let y = b.evaluate(beta, 1.0).unwrap();
        assert!((x.magnetization - y.magnetization).abs() < 1e-12);
        assert!((x.probabilities.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }
    assert!(a.evaluate(f64::MAX, f64::MAX).is_err());
}
#[test]
fn disconnected_support_stays_explicit_and_bounded() {
    struct Fixed;
    impl State for Fixed {
        type Params = ();
        type Change = ();
        fn energy(&self, _: &mut ()) -> f64 {
            0.0
        }
        fn propose_change(&self, _: &mut impl rand::RngExt) {}
        fn apply_change(&mut self, _: ()) {}
        fn revert_change(&mut self, _: ()) {}
    }
    let result = run(
        Fixed,
        (),
        rand::rng(),
        RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0, 1.0]).unwrap()).unwrap(),
        Config {
            preliminary_stages: 1,
            min_stage_steps: 1,
            visits_per_bin: 1,
            sampling_steps: 100,
            max_steps: 10,
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(
        result.diagnostics.stop_reason,
        Some(StopReason::BudgetExhausted)
    );
    assert_eq!(result.data.accessible(), &[true, true]);
    assert_eq!(result.data.lifetime_bins(), &[10, 0]);
    assert!(!result.is_complete());
}
