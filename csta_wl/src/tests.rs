use std::{
    convert::Infallible,
    ops::{Deref, DerefMut},
};

use super::*;
use rand::{Rng, RngExt, SeedableRng, TryRng, rngs::StdRng};

struct StepRng(u64);

impl StepRng {
    fn new(value: u64, _: u64) -> Self {
        Self(value)
    }
}

impl TryRng for StepRng {
    type Error = Infallible;

    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> std::prelude::v1::Result<(), Self::Error> {
        for (i, b) in dst.iter_mut().enumerate() {
            *b = self.0.to_le_bytes()[i % 8];
        }
        Ok(())
    }

    fn try_next_u32(&mut self) -> std::prelude::v1::Result<u32, Self::Error> {
        Ok((self.0 >> 32) as u32)
    }

    fn try_next_u64(&mut self) -> std::prelude::v1::Result<u64, Self::Error> {
        Ok(self.0)
    }
}

// impl Deref for StepRng {
//     type Target = Self;
//
//     fn deref(&self) -> &Self::Target {
//         self
//     }
// }
//
// impl DerefMut for StepRng {
//     fn deref_mut(&mut self) -> &mut Self::Target {
//         self
//     }
// }

#[derive(Debug)]
pub(crate) struct CounterState(pub i32);
impl State for CounterState {
    type Params = i32;
    type Change = i32;
    fn energy(&self, cache: &mut i32) -> f64 {
        *cache = self.0;
        f64::from(self.0)
    }
    fn propose_change(&self, rng: &mut impl Rng) -> i32 {
        if rng.random_bool(0.5) { 1 } else { -1 }
    }
    fn apply_change(&mut self, change: i32) {
        self.0 += change;
    }
    fn revert_change(&mut self, change: i32) {
        self.0 -= change;
    }
}
pub(crate) fn simple_config(steps: u64) -> Config {
    Config {
        preliminary_stages: 0,
        min_stage_steps: 1,
        visits_per_bin: 1,
        sampling_steps: steps,
        max_steps: 10000,
        ..Config::default()
    }
}
fn grid() -> EnergyGrid {
    EnergyGrid::discrete(vec![0.0, 1.0, 2.0]).unwrap()
}
fn walker() -> Walker<CounterState, StdRng> {
    Walker::new(
        CounterState(1),
        0,
        StdRng::seed_from_u64(42),
        RawWangLandauData::on_grid(grid()).unwrap(),
        simple_config(20),
    )
    .unwrap()
}
#[test]
fn invalid_grids_and_boundary_neighbors() {
    for (a, b, n) in [
        (0.0, 1.0, 0),
        (1.0, 1.0, 1),
        (2.0, 1.0, 2),
        (f64::NAN, 1.0, 2),
        (0.0, f64::INFINITY, 2),
        (-f64::MAX, f64::MAX, 2),
        (0.0, f64::from_bits(1), 2),
    ] {
        assert!(EnergyGrid::continuous(a, b, n).is_err());
    }
    for e in [
        vec![],
        vec![1.0, 1.0],
        vec![2.0, 1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
    ] {
        assert!(EnergyGrid::discrete(e).is_err());
    }
    let g = EnergyGrid::continuous(0.0, 2.0, 2).unwrap();
    for (e, expected) in [
        (0.0, 0),
        (1_f64.next_down(), 0),
        (1.0, 1),
        (1_f64.next_up(), 1),
        (2.0, 1),
    ] {
        assert_eq!(g.bin(e).unwrap(), Some(expected));
    }
    for e in [0_f64.next_down(), 2_f64.next_up(), -100.0, 100.0] {
        assert_eq!(g.bin(e).unwrap(), None);
    }
    for e in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(g.bin(e).is_err());
    }
    let single = EnergyGrid::continuous(-1.0, 1.0, 1).unwrap();
    for e in [-1.0, 0.0, 1.0] {
        assert_eq!(single.bin(e).unwrap(), Some(0));
    }
    for (i, e) in grid().energies().iter().enumerate() {
        assert_eq!(grid().bin(*e).unwrap(), Some(i));
    }
    assert_eq!(grid().bin(0.5).unwrap(), None);
}
#[test]
fn support_dos_and_counter_validation() {
    for d in [
        vec![],
        vec![f64::NAN; 3],
        vec![f64::INFINITY; 3],
        vec![f64::NEG_INFINITY; 3],
    ] {
        assert!(RawWangLandauData::with_support(grid(), vec![true; 3], d).is_err());
    }
    let mut d = RawWangLandauData::with_support(
        grid(),
        vec![true, false, true],
        vec![0.0, f64::NEG_INFINITY, 0.0],
    )
    .unwrap();
    assert!(!d.is_flat());
    d.record(0, 1.0);
    assert!(!d.is_flat());
    d.record(2, 1.0);
    assert!(d.is_flat());
    assert_eq!(d.energy_to_bin(1.0).unwrap(), None);
    let dos = d.dos.clone();
    d.clear_hist();
    assert_eq!(d.dos, dos);
    assert_eq!(d.total_visits, 2);
    assert_eq!(d.visits, 0);
    assert!(!d.is_flat());
    d.total_visits = u64::MAX;
    assert_eq!(d.check_record(0, 1.0), Err(Error::CounterOverflow));
    let mut c = simple_config(1);
    c.visits_per_bin = u64::MAX;
    assert_eq!(c.validate(2), Err(Error::CounterOverflow));
}
#[test]
fn flatness_threshold_is_inclusive() {
    let mut d = RawWangLandauData::new(2, 0.0, 2.0).unwrap();
    d.bins = vec![8, 12];
    assert!(d.is_flat());
    d.bins = vec![7, 13];
    assert!(!d.is_flat());
    assert!(!d.flat_at(f64::NAN));
}
#[test]
fn acceptance_thresholds_and_extreme_logs() {
    assert!(accept(0.0, &mut StepRng::new(u64::MAX, 0)));
    assert!(accept(10000.0, &mut StepRng::new(u64::MAX, 0)));
    assert!(!accept(-10000.0, &mut StepRng::new(u64::MAX, 0)));
    assert!(accept(0.5_f64.ln(), &mut StepRng::new(0, 0)));
    assert!(!accept(0.5_f64.ln(), &mut StepRng::new(u64::MAX, 0)));
    assert!(!accept(0.5_f64.ln(), &mut StepRng::new(1_u64 << 63, 0)));
}
#[test]
fn each_transition_updates_exactly_one_retained_bin() {
    let mut w = walker();
    for n in 1..=100 {
        w.transition(Some(1.0)).unwrap();
        assert_eq!(w.data.bins.iter().sum::<u64>(), n);
        assert_eq!(w.data.visits, n);
        assert_eq!(w.data.total_visits, n);
        assert_eq!(w.params, w.state.0);
        assert_eq!(w.energy, f64::from(w.state.0));
        assert_eq!(w.data.energy_to_bin(w.energy).unwrap(), Some(w.bin));
    }
}
#[derive(Debug)]
struct Proposal {
    energy: f64,
    next: f64,
}
impl State for Proposal {
    type Params = f64;
    type Change = f64;
    fn energy(&self, cache: &mut f64) -> f64 {
        *cache = self.energy;
        self.energy
    }
    fn propose_change(&self, _: &mut impl Rng) -> f64 {
        self.energy
    }
    fn apply_change(&mut self, _: f64) {
        self.energy = self.next;
    }
    fn revert_change(&mut self, old: f64) {
        self.energy = old;
    }
}
#[test]
fn rejected_and_invalid_moves_restore_state_and_parameter_caches() {
    for next in [-1.0, 3.0, f64::NAN, f64::INFINITY] {
        let mut w = Walker::new(
            Proposal { energy: 1.0, next },
            0.0,
            StepRng::new(0, 0),
            RawWangLandauData::on_grid(grid()).unwrap(),
            simple_config(1),
        )
        .unwrap();
        let result = w.transition(Some(0.5));
        assert_eq!(
            (w.state.energy, w.energy, w.params, w.bin),
            (1.0, 1.0, 1.0, 1)
        );
        if next.is_finite() {
            assert!(!result.unwrap());
            assert_eq!(w.data.bins(), &[0, 1, 0]);
            assert_eq!(w.data.dos()[1], 0.5);
        } else {
            assert!(result.is_err());
            assert_eq!(w.data.total_visits(), 0);
        }
    }
    let mut w = Walker::new(
        Proposal {
            energy: 1.0,
            next: 2.0,
        },
        0.0,
        StepRng::new(u64::MAX, 0),
        RawWangLandauData::with_support(grid(), vec![true; 3], vec![0.0, 0.0, 1000.0]).unwrap(),
        simple_config(1),
    )
    .unwrap();
    assert!(!w.transition(Some(0.5)).unwrap());
    assert_eq!(w.params, 1.0);
    assert_eq!(w.state.energy, 1.0);
    let data = RawWangLandauData::new(1, 0.0, 2.0).unwrap();
    let mut same = Walker::new(
        Proposal {
            energy: 0.5,
            next: 1.5,
        },
        0.0,
        StepRng::new(0, 0),
        data,
        simple_config(1),
    )
    .unwrap();
    assert!(same.transition(Some(1.0)).unwrap());
    assert_eq!(same.data.bins(), &[1]);
    assert_eq!(same.energy, 1.5);
}
#[test]
fn initialization_and_invalid_controls() {
    let data = RawWangLandauData::on_grid(grid()).unwrap();
    assert!(matches!(
        Walker::new(
            CounterState(9),
            0,
            StdRng::seed_from_u64(1),
            data.clone(),
            simple_config(1)
        ),
        Err(Error::InitialStateOutside)
    ));
    assert!(sample_in_support::<Proposal>(&mut StdRng::seed_from_u64(1), &0.0, &data, 3).is_err());
    assert!(sample_in_support::<Proposal>(&mut StdRng::seed_from_u64(1), &0.0, &data, 0).is_err());
    for v in [-1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(target_steps(v, 3).is_err());
    }
    assert_eq!(target_steps(0.0, 3).unwrap(), 0);
    for i in 0..8 {
        let mut c = simple_config(1);
        match i {
            0 => c.min_stage_steps = 0,
            1 => c.visits_per_bin = 0,
            2 => c.flatness = 0.0,
            3 => c.flatness = f64::NAN,
            4 => c.t0 = Some(0.0),
            5 => c.t1 = Some(f64::INFINITY),
            6 => c.initial_ln_f = -1.0,
            _ => c.max_steps = 0,
        }
        assert!(c.validate(3).is_err());
    }
}
#[test]
fn stage_minima_schedule_and_exact_stopping() {
    let mut c = simple_config(3);
    c.preliminary_stages = 2;
    c.min_stage_steps = 2;
    let d = RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0]).unwrap()).unwrap();
    let r = run(CounterState(0), 0, StdRng::seed_from_u64(1), d, c, None).unwrap();
    assert!(r.is_complete());
    assert_eq!(r.diagnostics.proposals, 7);
    assert_eq!(r.diagnostics.completed_stages, 2);
    assert_eq!(r.data.bins(), &[3]);
    assert_eq!(r.data.total_visits(), 7);
    assert!((r.diagnostics.last_update - 1.0 / 12.0).abs() < 1e-12);
    let expected = 2.0 + 1.0 + 1.0 / 10.0 + 1.0 / 11.0 + 1.0 / 12.0;
    assert!((r.data.dos()[0] - expected).abs() < 1e-12);
    // Flatness cannot override a longer explicit minimum.
    let mut w = walker();
    w.config.preliminary_stages = 1;
    w.diagnostics.phase = Phase::Preliminary;
    w.minimum_visits = 10;
    w.data.bins = vec![3; 3];
    w.data.visits = 9;
    assert!(w.data.is_flat());
    assert_eq!(w.diagnostics.completed_stages, 0);
    w.advance(None).unwrap();
    assert_eq!(w.diagnostics.completed_stages, 1);
}
#[test]
fn cancellation_zero_target_budget_and_seeded_repetition() {
    let execute = |config, cancel| {
        run(
            CounterState(1),
            0,
            StdRng::seed_from_u64(123),
            RawWangLandauData::on_grid(grid()).unwrap(),
            config,
            cancel,
        )
        .unwrap()
    };
    let a = execute(simple_config(100), None);
    let b = execute(simple_config(100), None);
    assert_eq!(a.data.dos(), b.data.dos());
    assert_eq!(a.diagnostics, b.diagnostics);
    let flag = AtomicBool::new(true);
    let r = execute(simple_config(100), Some(&flag));
    assert_eq!(r.diagnostics.stop_reason, Some(StopReason::Cancelled));
    assert_eq!(r.diagnostics.proposals, 0);
    let mut c = simple_config(0);
    c.preliminary_stages = 5;
    let r = execute(c, None);
    assert!(r.is_complete());
    assert_eq!(r.data.total_visits(), 0);
    let mut c = simple_config(100);
    c.max_steps = 3;
    let r = execute(c, None);
    assert_eq!(r.diagnostics.stop_reason, Some(StopReason::BudgetExhausted));
    assert_eq!(r.diagnostics.proposals, 3);
    let mut c = simple_config(1);
    c.preliminary_stages = 1;
    c.max_steps = 5;
    let r = run(
        Proposal {
            energy: 1.0,
            next: 1.0,
        },
        0.0,
        StdRng::seed_from_u64(1),
        RawWangLandauData::on_grid(grid()).unwrap(),
        c,
        None,
    )
    .unwrap();
    assert_eq!(r.diagnostics.stop_reason, Some(StopReason::BudgetExhausted));
    assert_eq!(r.diagnostics.phase, Phase::Preliminary);
}
#[test]
fn ising_seeded_wl_smoke() {
    use models::Ising;
    let exact = Ising::<4>::exact_dos().unwrap();
    for seed in [7, 42, 987] {
        let c = Config {
            preliminary_stages: 3,
            min_stage_steps: 100,
            visits_per_bin: 20,
            sampling_steps: 30_000,
            max_steps: 100_000,
            ..Config::default()
        };
        let r = run(
            Ising::<4>::new([1; 4]).unwrap(),
            (),
            StdRng::seed_from_u64(seed),
            RawWangLandauData::on_grid(Ising::<4>::grid().unwrap()).unwrap(),
            c,
            None,
        )
        .unwrap();
        assert!(r.is_complete());
        let mut d = r.data.process_data().unwrap();
        d.normalize_log_count(16_f64.ln()).unwrap();
        for (g, reference) in d.dos().iter().zip(exact.dos()) {
            assert!(
                (g - reference).abs() < 0.25,
                "seed={seed}: {g} vs {reference}"
            );
        }
    }
}

impl Randomizable for Proposal {
    fn sample<R: Rng + ?Sized>(_: &mut R) -> Self {
        Self {
            energy: 9.0,
            next: 9.0,
        }
    }
}

#[test]
fn counter_and_dos_failure_roll_back_without_partial_updates() {
    for error in [
        Error::CounterOverflow,
        Error::Numerical("DOS update overflow"),
        Error::Invalid("DOS update must be finite and nonnegative"),
    ] {
        let mut w = Walker::new(
            Proposal {
                energy: 1.0,
                next: 2.0,
            },
            0.0,
            StepRng::new(0, 0),
            RawWangLandauData::on_grid(grid()).unwrap(),
            simple_config(1),
        )
        .unwrap();
        let delta = match error {
            Error::CounterOverflow => {
                w.data.total_visits = u64::MAX;
                1.0
            }
            Error::Numerical(_) => {
                w.data.dos.fill(f64::MAX);
                f64::MAX
            }
            _ => -1.0,
        };
        let dos = w.data.dos.clone();
        let bins = w.data.bins.clone();
        let visits = w.data.total_visits;
        assert_eq!(w.transition(Some(delta)), Err(error));
        assert_eq!(
            (w.state.energy, w.params, w.energy, w.bin),
            (1.0, 1.0, 1.0, 1)
        );
        assert_eq!(w.data.dos, dos);
        assert_eq!(w.data.bins, bins);
        assert_eq!(w.data.total_visits, visits);
    }
}
#[test]
fn phase_switch_retains_state_and_minimums_are_all_required() {
    let mut w = walker();
    w.diagnostics.phase = Phase::Preliminary;
    w.config.preliminary_stages = 1;
    w.minimum_visits = 2;
    w.data.visits = 10;
    w.data.bins = vec![10, 0, 0];
    w.advance(None).unwrap();
    assert_eq!(w.diagnostics.completed_stages, 0);
    let data = RawWangLandauData::on_grid(EnergyGrid::discrete(vec![1.0]).unwrap()).unwrap();
    let mut c = simple_config(1);
    c.preliminary_stages = 1;
    c.min_stage_steps = 2;
    let mut w = Walker::new(
        Proposal {
            energy: 1.0,
            next: 1.0,
        },
        0.0,
        StepRng::new(0, 0),
        data,
        c,
    )
    .unwrap();
    w.advance(None).unwrap();
    assert_eq!(w.diagnostics.completed_stages, 0);
    w.advance(None).unwrap();
    assert_eq!(w.diagnostics.phase, Phase::Sampling);
    assert_eq!(w.data.dos(), &[2.0]);
    assert_eq!(w.data.bins(), &[0]);
    assert_eq!(w.data.total_visits(), 2);
    assert_eq!((w.state.energy, w.energy, w.params), (1.0, 1.0, 1.0));
    w.advance(None).unwrap();
    assert_eq!(w.diagnostics.proposals, 3);
    assert_eq!(w.diagnostics.stop_reason, Some(StopReason::TargetReached));
}
#[test]
fn target_time_boundary_and_gamma_use_proposals() {
    for (time, expected) in [(1_f64.next_down(), 2), (1.0, 2), (1_f64.next_up(), 3)] {
        assert_eq!(target_steps(time, 2).unwrap(), expected);
    }
    let mut w = walker();
    w.advance(None).unwrap();
    assert!((w.diagnostics.last_update - 3.0 / 30.0).abs() < 1e-15);
    w.advance(None).unwrap();
    assert!((w.diagnostics.last_update - 3.0 / 31.0).abs() < 1e-15); // NOT 3/(30+1/3)
    let c = Config {
        initial_ln_f: f64::from_bits(1),
        preliminary_stages: 2,
        ..simple_config(1)
    };
    assert!(c.validate(1).is_err());
}
#[test]
fn oscillator_seeded_wl_matches_truncated_reference() {
    use models::Oscillators;
    let c = Config {
        preliminary_stages: 3,
        min_stage_steps: 100,
        visits_per_bin: 20,
        sampling_steps: 50_000,
        max_steps: 100_000,
        ..Config::default()
    };
    let r = run(
        Oscillators::<2>::new([0; 2]).unwrap(),
        (),
        StdRng::seed_from_u64(93),
        RawWangLandauData::on_grid(Oscillators::<2>::grid(5).unwrap()).unwrap(),
        c,
        None,
    )
    .unwrap();
    assert!(r.is_complete());
    let mut data = r.data.process_data().unwrap();
    let exact = Oscillators::<2>::exact_dos(5).unwrap();
    data.normalize_log_count(log_sum_exp(exact.dos()).unwrap())
        .unwrap();
    for (a, b) in data.dos().iter().zip(exact.dos()) {
        assert!((a - b).abs() < 0.25, "{a} != {b}");
    }
}

#[test]
#[ignore = "longer seeded convergence study; run explicitly with --ignored"]
fn ising_long_convergence_study() {
    use models::Ising;
    let exact = Ising::<8>::exact_dos().unwrap();
    for seed in 0..10 {
        let c = Config {
            sampling_steps: 1_000_000,
            max_steps: 2_000_000,
            ..Config::default()
        };
        let r = run(
            Ising::<8>::new([1; 8]).unwrap(),
            (),
            StdRng::seed_from_u64(seed),
            RawWangLandauData::on_grid(Ising::<8>::grid().unwrap()).unwrap(),
            c,
            None,
        )
        .unwrap();
        assert!(r.is_complete());
        let mut data = r.data.process_data().unwrap();
        data.normalize_log_count(256_f64.ln()).unwrap();
        for (a, b) in data.dos().iter().zip(exact.dos()) {
            assert!((a - b).abs() < 0.15, "seed={seed}: {a} != {b}");
        }
    }
}

#[test]
fn cancellation_during_run_stops_at_next_proposal_boundary() {
    use std::sync::Arc;
    struct CancelState {
        flag: Arc<AtomicBool>,
        changes: u64,
    }
    impl State for CancelState {
        type Params = ();
        type Change = ();
        fn energy(&self, _: &mut ()) -> f64 {
            0.0
        }
        fn propose_change(&self, _: &mut impl Rng) {}
        fn apply_change(&mut self, _: ()) {
            self.changes += 1;
            if self.changes == 3 {
                self.flag.store(true, Ordering::Relaxed);
            }
        }
        fn revert_change(&mut self, _: ()) {
            unreachable!()
        }
    }
    let flag = Arc::new(AtomicBool::new(false));
    let data = RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0]).unwrap()).unwrap();
    let r = run(
        CancelState {
            flag: flag.clone(),
            changes: 0,
        },
        (),
        StdRng::seed_from_u64(1),
        data,
        simple_config(100),
        Some(&flag),
    )
    .unwrap();
    assert_eq!(r.diagnostics.stop_reason, Some(StopReason::Cancelled));
    assert_eq!(r.diagnostics.proposals, 3);
    assert_eq!(r.data.total_visits(), 3);
}
#[test]
fn invalid_extreme_allocations_and_unrepresentable_updates() {
    assert!(EnergyGrid::continuous(0.0, 1.0, usize::MAX).is_err());
    let d = RawWangLandauData::with_support(grid(), vec![true; 3], vec![1e300; 3]).unwrap();
    assert!(matches!(d.check_record(0, 0.1), Err(Error::Numerical(_))));
    assert!(
        RawWangLandauData::with_support(grid(), vec![false; 3], vec![f64::NEG_INFINITY; 3])
            .is_err()
    );
    let mut d = RawWangLandauData::on_grid(grid()).unwrap();
    d.bins[1] = u64::MAX;
    assert_eq!(d.check_record(1, 0.1), Err(Error::CounterOverflow));
    d.bins[1] = 0;
    d.lifetime_bins[1] = u64::MAX;
    assert_eq!(d.check_record(1, 0.1), Err(Error::CounterOverflow));
    d.lifetime_bins[1] = 0;
    d.visits = u64::MAX;
    assert_eq!(d.check_record(1, 0.1), Err(Error::CounterOverflow));
}
