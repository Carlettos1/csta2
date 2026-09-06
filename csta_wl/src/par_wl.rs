//! Replica exchange and derivative-matched DOS stitching, following
//! arXiv:2103.15028v2 section 2, combined with the serial SAMC schedule.
use crate::*;
use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{
    ops::Range,
    panic::{AssertUnwindSafe, catch_unwind},
};

#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[derive(Clone, Debug, PartialEq)]
pub struct ParallelConfig {
    /// Ordered half-open GLOBAL BIN ranges on a shared grid. Cover the grid;
    /// adjacent windows must overlap in at least two bins and extend coverage.
    pub windows: Vec<Range<usize>>,
    /// Independent DOS estimates, averaged after anchoring within each window.
    pub walkers_per_window: usize,
    pub exchange_every: u64,
    pub seed: u64,
    pub run: Config,
}
impl ParallelConfig {
    fn validate(&self, bins: usize) -> Result<()> {
        if self.windows.is_empty() || self.walkers_per_window == 0 || self.exchange_every == 0 {
            return Err(Error::Invalid(
                "need windows, walkers and positive exchange interval",
            ));
        }
        self.windows
            .len()
            .checked_mul(self.walkers_per_window)
            .ok_or(Error::CounterOverflow)?;
        for (i, w) in self.windows.iter().enumerate() {
            if w.start >= w.end || w.end > bins {
                return Err(Error::Invalid("invalid window range"));
            }
            self.run.validate(w.len())?;
            if i == 0 && w.start != 0 {
                return Err(Error::Invalid("first window must start at bin zero"));
            }
            if i > 0 {
                let prev = &self.windows[i - 1];
                if w.start <= prev.start
                    || w.end <= prev.end
                    || w.start >= prev.end
                    || prev.end - w.start < 2
                {
                    return Err(Error::Invalid(
                        "windows must extend coverage with at least two overlap bins",
                    ));
                }
            }
        }
        if self.windows.last().unwrap().end != bins {
            return Err(Error::Invalid("windows do not cover the grid"));
        }
        Ok(())
    }
}
pub struct ParallelResult<S: State> {
    pub walkers: Vec<RunResult<S>>,
    /// Present only when every walker reaches its target and stitching succeeds.
    pub merged: Option<WLData>,
    pub stop_reason: StopReason,
    /// Eligible exchange trials; not included in histogram/proposal counters.
    pub exchange_attempts: u64,
    pub accepted_exchanges: u64,
}

/// Initial states are window-major, then replica-major. All replicas use the same
/// Hamiltonian parameters (cloned from `params`). State and parameter caches
/// travel together on exchange; RNG, DOS and adaptation stay in their windows.
/// Fixed seed/config/initial states are reproducible independent of Rayon thread
/// count, on the same rand version/platform. Model methods must terminate and
/// must not share mutable state across clones. Panics become WorkerFailed.
/// A failed/cancelled/budget-exhausted run never returns a merged DOS.
pub fn run_parallel<S>(
    grid: EnergyGrid,
    states: Vec<S>,
    params: S::Params,
    config: ParallelConfig,
    cancel: Option<&AtomicBool>,
) -> Result<ParallelResult<S>>
where
    S: State + Send,
    S::Params: Clone + Send,
{
    let mut session = ParallelSession::<S, StdRng>::new(grid, states, params, config)?;
    session.advance(u64::MAX, cancel)?;
    session.finish()
}

/// Incremental parallel execution. Checkpoints are taken between complete chunks,
/// after deterministic exchanges; default RNG matches `run_parallel`.
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "checkpoint",
    serde(bound(
        serialize = "S: serde::Serialize, S::Params: serde::Serialize, R: serde::Serialize",
        deserialize = "S: serde::Deserialize<'de>, S::Params: serde::Deserialize<'de>, R: serde::Deserialize<'de>"
    ))
)]
pub struct ParallelSession<S: State, R = StdRng> {
    walkers: Vec<Walker<S, R>>,
    exchange_rng: R,
    config: ParallelConfig,
    parity: usize,
    exchange_attempts: u64,
    accepted_exchanges: u64,
    stop_reason: Option<StopReason>,
    failed: bool,
}
impl<S: State + Send, R: RngExt + SeedableRng + Send> ParallelSession<S, R>
where
    S::Params: Clone + Send,
{
    pub fn new(
        grid: EnergyGrid,
        states: Vec<S>,
        params: S::Params,
        config: ParallelConfig,
    ) -> Result<Self> {
        config.validate(grid.len())?;
        let count = config.windows.len() * config.walkers_per_window;
        if states.len() != count {
            return Err(Error::Invalid("one initial state is required per walker"));
        }
        let mut seed_rng = R::seed_from_u64(config.seed);
        let mut walkers = Vec::with_capacity(count);
        for (i, state) in states.into_iter().enumerate() {
            let range = &config.windows[i / config.walkers_per_window];
            let mask: Vec<_> = (0..grid.len()).map(|b| range.contains(&b)).collect();
            let dos = mask
                .iter()
                .map(|a| if *a { 0.0 } else { f64::NEG_INFINITY })
                .collect();
            let data = RawWangLandauData::with_support(grid.clone(), mask, dos)?;
            let rng = R::seed_from_u64(seed_rng.random());
            walkers.push(
                catch_unwind(AssertUnwindSafe(|| {
                    Walker::new(state, params.clone(), rng, data, config.run.clone())
                }))
                .map_err(|_| Error::WorkerFailed)??,
            );
        }
        Ok(Self {
            walkers,
            exchange_rng: R::seed_from_u64(seed_rng.random()),
            config,
            parity: 0,
            exchange_attempts: 0,
            accepted_exchanges: 0,
            stop_reason: None,
            failed: false,
        })
    }
    pub fn diagnostics(&self) -> Vec<&Diagnostics> {
        self.walkers.iter().map(|w| &w.diagnostics).collect()
    }
    pub fn advance(
        &mut self,
        chunks: u64,
        cancel: Option<&AtomicBool>,
    ) -> Result<Option<StopReason>> {
        if self.failed {
            return Err(Error::WorkerFailed);
        }
        if chunks == 0 {
            return Ok(self.stop_reason);
        }
        if self.stop_reason == Some(StopReason::Cancelled) {
            self.stop_reason = None;
            for w in &mut self.walkers {
                if w.diagnostics.stop_reason == Some(StopReason::Cancelled) {
                    w.diagnostics.stop_reason = None;
                }
            }
        }
        if self.stop_reason.is_some() {
            return Ok(self.stop_reason);
        }
        for _ in 0..chunks {
            match self.chunk(cancel) {
                Ok(Some(reason)) => {
                    self.stop_reason = Some(reason);
                    break;
                }
                Ok(None) => {}
                Err(e) => {
                    self.failed = true;
                    return Err(e);
                }
            }
        }
        Ok(self.stop_reason)
    }
    fn chunk(&mut self, cancel: Option<&AtomicBool>) -> Result<Option<StopReason>> {
        if cancelled(cancel) {
            return Ok(Some(StopReason::Cancelled));
        }
        if self
            .walkers
            .iter()
            .all(|w| w.diagnostics.stop_reason == Some(StopReason::TargetReached))
        {
            return Ok(Some(StopReason::TargetReached));
        }
        let failed = AtomicBool::new(false);
        let results: Vec<Result<()>> = self
            .walkers
            .par_iter_mut()
            .map(|w| {
                let outcome = catch_unwind(AssertUnwindSafe(|| {
                    for _ in 0..self.config.exchange_every.min(self.config.run.max_steps) {
                        if cancelled(cancel) || failed.load(Ordering::Relaxed) {
                            break;
                        }
                        match w.diagnostics.stop_reason {
                            Some(StopReason::TargetReached) => {
                                let next = w
                                    .diagnostics
                                    .mixing_steps
                                    .checked_add(1)
                                    .ok_or(Error::CounterOverflow)?;
                                w.transition(None)?;
                                w.diagnostics.mixing_steps = next;
                            }
                            Some(_) => break,
                            None => {
                                w.advance(cancel)?;
                                if w.diagnostics.stop_reason.is_some() {
                                    break;
                                }
                            }
                        }
                    }
                    Ok(())
                }))
                .map_err(|_| Error::WorkerFailed)
                .and_then(|r| r);
                if outcome.is_err() {
                    failed.store(true, Ordering::Relaxed);
                }
                outcome
            })
            .collect();
        for r in results {
            r?;
        }
        if cancelled(cancel) {
            return Ok(Some(StopReason::Cancelled));
        }
        if self
            .walkers
            .iter()
            .any(|w| w.diagnostics.stop_reason == Some(StopReason::BudgetExhausted))
        {
            return Ok(Some(StopReason::BudgetExhausted));
        }
        if self
            .walkers
            .iter()
            .all(|w| w.diagnostics.stop_reason == Some(StopReason::TargetReached))
        {
            return Ok(Some(StopReason::TargetReached));
        }
        for i in (self.parity..self.config.windows.len().saturating_sub(1)).step_by(2) {
            for replica in 0..self.config.walkers_per_window {
                if self.exchange_attempts == u64::MAX || self.accepted_exchanges == u64::MAX {
                    return Err(Error::CounterOverflow);
                }
                let a = i * self.config.walkers_per_window + replica;
                let b = (i + 1) * self.config.walkers_per_window + replica;
                let (left, right) = self.walkers.split_at_mut(b);
                if let Some(accepted) =
                    exchange(&mut left[a], &mut right[0], &mut self.exchange_rng)?
                {
                    self.exchange_attempts += 1;
                    self.accepted_exchanges += u64::from(accepted);
                }
            }
        }
        self.parity ^= 1;
        Ok(None)
    }
    pub fn finish(self) -> Result<ParallelResult<S>> {
        if self.failed {
            return Err(Error::WorkerFailed);
        }
        let stop_reason = self.stop_reason.unwrap_or(StopReason::Cancelled);
        let mut runs: Vec<_> = self.walkers.into_iter().map(Walker::finish).collect();
        for r in &mut runs {
            if r.diagnostics.stop_reason.is_none() {
                r.diagnostics.stop_reason = Some(stop_reason);
            }
        }
        let merged =
            if stop_reason == StopReason::TargetReached && self.config.run.sampling_steps > 0 {
                let mut fragments = Vec::new();
                for group in runs.chunks(self.config.walkers_per_window) {
                    fragments.push(average_replicas(
                        &group.iter().map(|r| &r.data).collect::<Vec<_>>(),
                    )?);
                }
                Some(merge_windows(&fragments)?)
            } else {
                None
            };
        Ok(ParallelResult {
            walkers: runs,
            merged,
            stop_reason,
            exchange_attempts: self.exchange_attempts,
            accepted_exchanges: self.accepted_exchanges,
        })
    }
}
#[cfg(feature = "checkpoint")]
impl<S: State + Send, R: RngExt + SeedableRng + Send> ParallelSession<S, R>
where
    S: serde::Serialize + serde::de::DeserializeOwned,
    S::Params: Clone + Send + serde::Serialize + serde::de::DeserializeOwned,
    R: serde::Serialize + serde::de::DeserializeOwned,
{
    fn validate(&self) -> Result<()> {
        let first = self
            .walkers
            .first()
            .ok_or(Error::Invalid("empty checkpoint"))?;
        self.config.validate(first.data.grid.len())?;
        if self.failed
            || self.parity > 1
            || self.accepted_exchanges > self.exchange_attempts
            || self.walkers.len() != self.config.windows.len() * self.config.walkers_per_window
        {
            return Err(Error::Invalid("invalid parallel checkpoint"));
        }
        for (i, w) in self.walkers.iter().enumerate() {
            crate::session::validate_walker(w, true)?;
            let range = &self.config.windows[i / self.config.walkers_per_window];
            if w.config != self.config.run
                || w.data.grid != first.data.grid
                || w.data
                    .accessible
                    .iter()
                    .enumerate()
                    .any(|(b, a)| *a != range.contains(&b))
            {
                return Err(Error::Invalid("parallel checkpoint support mismatch"));
            }
        }
        if self.stop_reason == Some(StopReason::TargetReached)
            && self
                .walkers
                .iter()
                .any(|w| w.diagnostics.stop_reason != Some(StopReason::TargetReached))
        {
            return Err(Error::Invalid("incomplete parallel checkpoint"));
        }
        Ok(())
    }
    pub fn save(
        &self,
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.validate()?;
        csta_core::checkpoint::save(path, model, self)
    }
    pub fn load(
        path: impl AsRef<std::path::Path>,
        model: &str,
    ) -> std::result::Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let s: Self = csta_core::checkpoint::load(path, model)?;
        s.validate()?;
        Ok(s)
    }
}

fn exchange<S: State, R: RngExt>(
    a: &mut Walker<S, R>,
    b: &mut Walker<S, R>,
    rng: &mut impl RngExt,
) -> Result<Option<bool>> {
    if a.data.grid != b.data.grid {
        return Err(Error::Invalid("exchange grids differ"));
    }
    if !a.data.accessible[b.bin] || !b.data.accessible[a.bin] {
        return Ok(None);
    }
    // Eq. (4), evaluated as differences to remove arbitrary window offsets.
    let ratio = (a.data.dos[a.bin] - a.data.dos[b.bin]) + (b.data.dos[b.bin] - b.data.dos[a.bin]);
    if ratio.is_nan() {
        return Err(Error::Numerical("indeterminate exchange ratio"));
    }
    let accepted = accept(ratio, rng);
    if accepted {
        std::mem::swap(&mut a.state, &mut b.state);
        std::mem::swap(&mut a.params, &mut b.params);
        std::mem::swap(&mut a.energy, &mut b.energy);
        std::mem::swap(&mut a.bin, &mut b.bin);
    }
    Ok(Some(accepted))
}

fn average_replicas(data: &[&RawWangLandauData]) -> Result<RawWangLandauData> {
    let first = data.first().ok_or(Error::Invalid("no replicas"))?;
    let mut out = (*first).clone();
    out.bins.fill(0);
    out.lifetime_bins.fill(0);
    out.visits = 0;
    out.total_visits = 0;
    let anchor = out.accessible.iter().position(|v| *v).unwrap();
    for i in 0..out.dos.len() {
        if !out.accessible[i] {
            continue;
        }
        out.dos[i] = 0.0;
        for d in data {
            if d.grid != out.grid || d.accessible != out.accessible {
                return Err(Error::Invalid("replica supports differ"));
            }
            if d.lifetime_bins[i] == 0 {
                return Err(Error::InsufficientOverlap);
            }
            out.dos[i] += (d.dos[i] - d.dos[anchor]) / data.len() as f64;
            out.bins[i] = out.bins[i]
                .checked_add(d.bins[i])
                .ok_or(Error::CounterOverflow)?;
            out.lifetime_bins[i] = out.lifetime_bins[i]
                .checked_add(d.lifetime_bins[i])
                .ok_or(Error::CounterOverflow)?;
        }
        if !out.dos[i].is_finite() {
            return Err(Error::Numerical("replica DOS averaging overflow"));
        }
    }
    for d in data {
        out.visits = out
            .visits
            .checked_add(d.visits)
            .ok_or(Error::CounterOverflow)?;
        out.total_visits = out
            .total_visits
            .checked_add(d.total_visits)
            .ok_or(Error::CounterOverflow)?;
    }
    Ok(out)
}

/// Paste at the overlap point with smallest entropy-slope mismatch, then shift
/// the right fragment to agree there. This follows the paper/rwl strategy, not
/// a simple average across windows. All included levels must have been visited.
/// Histograms sum actual sampling visits across windows, regardless of paste.
pub fn merge_windows(fragments: &[RawWangLandauData]) -> Result<WLData> {
    let first = fragments
        .first()
        .ok_or(Error::Invalid("no DOS fragments"))?;
    let grid = first.grid.clone();
    let mut ranges = Vec::new();
    for f in fragments {
        if f.grid != grid {
            return Err(Error::Invalid("fragment grids differ"));
        }
        let start = f
            .accessible
            .iter()
            .position(|v| *v)
            .ok_or(Error::InsufficientOverlap)?;
        let end = f.accessible.iter().rposition(|v| *v).unwrap() + 1;
        if f.accessible[start..end].iter().any(|a| !a) || f.lifetime_bins[start..end].contains(&0) {
            return Err(Error::InsufficientOverlap);
        }
        ranges.push(start..end);
    }
    ParallelConfig {
        windows: ranges.clone(),
        walkers_per_window: 1,
        exchange_every: 1,
        seed: 0,
        run: Config::default(),
    }
    .validate(grid.len())?;
    let mut dos = first.dos.clone();
    let mut bins = first.bins.clone();
    let mut end = ranges[0].end;
    for (f, range) in fragments.iter().zip(&ranges).skip(1) {
        let e = grid.energies();
        let mut best: Option<(usize, f64)> = None;
        for i in range.start..end {
            let a = wl_data::derivative(&e[..end], &dos[..end], i)?;
            let b = wl_data::derivative(&e[range.clone()], &f.dos[range.clone()], i - range.start)?;
            let mismatch = (a - b).abs();
            if mismatch.is_finite() && best.is_none_or(|(_, score)| mismatch < score) {
                best = Some((i, mismatch));
            }
        }
        let (paste, _) = best.ok_or(Error::InsufficientOverlap)?;
        let anchor = dos[paste];
        for (i, value) in dos.iter_mut().enumerate().take(range.end).skip(paste) {
            *value = (f.dos[i] - f.dos[paste]) + anchor;
            if !value.is_finite() {
                return Err(Error::Numerical("merged DOS overflow"));
            }
        }
        for (i, value) in bins.iter_mut().enumerate() {
            *value = value.checked_add(f.bins[i]).ok_or(Error::CounterOverflow)?;
        }
        end = range.end;
    }
    let mut merged = WLData::new(grid, dos, bins)?;
    merged.anchor(0, 0.0)?;
    Ok(merged)
}

/// The old single-walker `par_wl` prototype is replaced by this explicit config API.
pub use run_parallel as par_wl;

#[cfg(test)]
mod tests {
    #![allow(clippy::single_range_in_vec_init)]
    use super::*;
    use crate::tests::{CounterState, simple_config};
    fn fragment(range: Range<usize>, offset: f64) -> RawWangLandauData {
        let grid = EnergyGrid::discrete((0..6).map(f64::from).collect()).unwrap();
        let mask: Vec<_> = (0..6).map(|i| range.contains(&i)).collect();
        let dos = (0..6)
            .map(|i| {
                if mask[i] {
                    i as f64 * 2.0 + offset
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let mut d = RawWangLandauData::with_support(grid, mask, dos).unwrap();
        for i in range {
            d.record(i, 0.0);
        }
        d
    }
    #[test]
    fn merges_known_offsets_and_checks_coverage() {
        let a = fragment(0..4, 100.0);
        let b = fragment(2..6, -10.0);
        assert_eq!(
            merge_windows(&[a.clone(), b.clone()]).unwrap().dos(),
            &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
        );
        let mut unvisited = b;
        unvisited.lifetime_bins[2] = 0;
        assert!(merge_windows(&[a, unvisited]).is_err());
        assert!(merge_windows(&[fragment(0..2, 0.0), fragment(3..6, 0.0)]).is_err());
        assert!(merge_windows(&[]).is_err());
    }
    fn config() -> ParallelConfig {
        ParallelConfig {
            windows: vec![0..4, 2..6],
            walkers_per_window: 1,
            exchange_every: 10,
            seed: 42,
            run: simple_config(500),
        }
    }
    #[test]
    fn invalid_window_configs() {
        let mut c = config();
        assert!(c.validate(6).is_ok());
        for windows in [
            vec![],
            vec![1..6],
            vec![0..4, 4..6],
            vec![0..5, 2..4],
            vec![0..3, 2..6],
            vec![0..7],
        ] {
            c.windows = windows;
            assert!(c.validate(6).is_err());
        }
        c = config();
        c.walkers_per_window = 0;
        assert!(c.validate(6).is_err());
        c = config();
        c.exchange_every = 0;
        assert!(c.validate(6).is_err());
    }
    #[test]
    fn exchanges_respect_support_and_keep_caches_with_state() {
        let mut a = Walker::new(
            CounterState(2),
            0,
            StdRng::seed_from_u64(1),
            fragment(0..4, 0.0),
            simple_config(2),
        );
        // Walker requires fresh counters; construct equivalent supports.
        assert!(a.is_err());
        let fresh = |r| {
            let d = fragment(r, 0.0);
            RawWangLandauData::with_support(d.grid.clone(), d.accessible.clone(), d.dos.clone())
                .unwrap()
        };
        a = Walker::new(
            CounterState(2),
            0,
            StdRng::seed_from_u64(1),
            fresh(0..4),
            simple_config(2),
        );
        let mut a = a.unwrap();
        let mut b = Walker::new(
            CounterState(3),
            0,
            StdRng::seed_from_u64(2),
            fresh(2..6),
            simple_config(2),
        )
        .unwrap();
        let mut rng = StdRng::seed_from_u64(5);
        assert_eq!(exchange(&mut a, &mut b, &mut rng).unwrap(), Some(true));
        assert_eq!((a.state.0, a.energy, a.params), (3, 3.0, 3));
        assert_eq!((b.state.0, b.energy, b.params), (2, 2.0, 2));
        assert_eq!(a.data.total_visits, 0);
        b.bin = 5;
        b.energy = 5.0;
        b.state.0 = 5;
        assert_eq!(exchange(&mut a, &mut b, &mut rng).unwrap(), None);
    }
    #[test]
    fn single_window_matches_serial_and_thread_count_is_reproducible() {
        let grid = EnergyGrid::discrete((0..6).map(f64::from).collect()).unwrap();
        let mut c = config();
        c.windows = vec![0..6];
        let mut seeder = StdRng::seed_from_u64(c.seed);
        let serial = run(
            CounterState(2),
            0,
            StdRng::seed_from_u64(seeder.random()),
            RawWangLandauData::on_grid(grid.clone()).unwrap(),
            c.run.clone(),
            None,
        )
        .unwrap();
        let parallel = run_parallel(grid.clone(), vec![CounterState(2)], 0, c, None).unwrap();
        assert_eq!(serial.data.dos(), parallel.walkers[0].data.dos());
        assert_eq!(serial.data.bins(), parallel.walkers[0].data.bins());
        let run_with = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    run_parallel(
                        grid.clone(),
                        vec![CounterState(1), CounterState(4)],
                        0,
                        config(),
                        None,
                    )
                    .unwrap()
                })
        };
        let a = run_with(1);
        let b = run_with(3);
        assert_eq!(a.merged, b.merged);
        assert_eq!(a.exchange_attempts, b.exchange_attempts);
        for (a, b) in a.walkers.iter().zip(&b.walkers) {
            assert_eq!(a.data.dos(), b.data.dos());
            assert_eq!(a.diagnostics, b.diagnostics);
            assert_eq!(a.data.total_visits(), 500);
        }
    }
    #[test]
    fn cancellation_and_budget_never_merge_partial_results() {
        let grid = EnergyGrid::discrete((0..6).map(f64::from).collect()).unwrap();
        let flag = AtomicBool::new(true);
        let r = run_parallel(
            grid.clone(),
            vec![CounterState(1), CounterState(4)],
            0,
            config(),
            Some(&flag),
        )
        .unwrap();
        assert_eq!(r.stop_reason, StopReason::Cancelled);
        assert!(r.merged.is_none());
        let mut c = config();
        c.run.max_steps = 3;
        let r = run_parallel(grid, vec![CounterState(1), CounterState(4)], 0, c, None).unwrap();
        assert_eq!(r.stop_reason, StopReason::BudgetExhausted);
        assert!(r.merged.is_none());
        assert!(r.walkers.iter().all(|w| w.diagnostics.proposals == 3));
    }
    #[test]
    fn worker_panic_is_an_error() {
        struct Broken;
        impl State for Broken {
            type Params = ();
            type Change = ();
            fn energy(&self, _: &mut ()) -> f64 {
                0.0
            }
            fn propose_change(&self, _: &mut impl RngExt) {
                panic!("mock model failure");
            }
            fn apply_change(&mut self, _: ()) {}
            fn revert_change(&mut self, _: ()) {}
        }
        let mut c = config();
        c.windows = vec![0..1];
        let r = run_parallel(
            EnergyGrid::discrete(vec![0.0]).unwrap(),
            vec![Broken],
            (),
            c,
            None,
        );
        assert!(matches!(r, Err(Error::WorkerFailed)));
    }
    #[test]
    fn rejected_exchange_leaves_both_walkers_unchanged() {
        let make = |state| {
            Walker::new(
                CounterState(state),
                0,
                StdRng::seed_from_u64(1),
                RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0, 1.0]).unwrap()).unwrap(),
                simple_config(10),
            )
            .unwrap()
        };
        let mut a = make(0);
        let mut b = make(1);
        a.data.dos = vec![0.0, 1000.0];
        b.data.dos = vec![1000.0, 0.0];
        assert_eq!(
            exchange(&mut a, &mut b, &mut StdRng::seed_from_u64(7)).unwrap(),
            Some(false)
        );
        assert_eq!((a.state.0, a.bin, a.params), (0, 0, 0));
        assert_eq!((b.state.0, b.bin, b.params), (1, 1, 1));
    }
    #[test]
    fn multiple_replicas_merge_and_count_visits() {
        let grid = EnergyGrid::discrete((0..6).map(f64::from).collect()).unwrap();
        let mut c = config();
        c.walkers_per_window = 2;
        let r = run_parallel(
            grid,
            vec![
                CounterState(1),
                CounterState(2),
                CounterState(3),
                CounterState(4),
            ],
            0,
            c,
            None,
        )
        .unwrap();
        assert_eq!(r.stop_reason, StopReason::TargetReached);
        assert_eq!(r.walkers.len(), 4);
        assert_eq!(r.merged.unwrap().bins().iter().sum::<u64>(), 2000);
        assert!(r.exchange_attempts > 0);
        assert!(r.accepted_exchanges <= r.exchange_attempts);
    }
    #[test]
    fn finished_window_mixes_without_adapting() {
        // The one-bin window completes immediately while the larger neighbor
        // must gather visits. Exercise the same frozen transition used by workers.
        let mut w = Walker::new(
            CounterState(1),
            0,
            StdRng::seed_from_u64(2),
            RawWangLandauData::on_grid(EnergyGrid::discrete(vec![0.0, 1.0, 2.0]).unwrap()).unwrap(),
            simple_config(1),
        )
        .unwrap();
        w.advance(None).unwrap();
        let dos = w.data.dos.clone();
        let hist = w.data.bins.clone();
        for _ in 0..100 {
            w.transition(None).unwrap();
        }
        assert_eq!(w.data.dos, dos);
        assert_eq!(w.data.bins, hist);
        assert_eq!(w.data.total_visits, 1);
        assert_eq!(w.params, w.state.0);
    }
    #[test]
    fn huge_exchange_interval_still_obeys_budget() {
        let grid = EnergyGrid::discrete((0..6).map(f64::from).collect()).unwrap();
        let mut c = config();
        c.exchange_every = u64::MAX;
        c.run.max_steps = 3;
        let r = run_parallel(grid, vec![CounterState(1), CounterState(4)], 0, c, None).unwrap();
        assert_eq!(r.stop_reason, StopReason::BudgetExhausted);
        assert!(r.walkers.iter().all(|w| w.diagnostics.proposals <= 3));
    }
    #[test]
    fn parallel_seeded_ising_matches_exact_reference() {
        use crate::models::Ising;
        let grid = Ising::<8>::grid().unwrap();
        let c = ParallelConfig {
            windows: vec![0..4, 2..5],
            walkers_per_window: 1,
            exchange_every: 100,
            seed: 42,
            run: Config {
                preliminary_stages: 3,
                min_stage_steps: 100,
                visits_per_bin: 20,
                sampling_steps: 50_000,
                max_steps: 200_000,
                ..Config::default()
            },
        };
        let r = run_parallel(
            grid,
            vec![
                Ising::<8>::new([1; 8]).unwrap(),
                Ising::<8>::new([1, -1, 1, -1, 1, -1, 1, -1]).unwrap(),
            ],
            (),
            c,
            None,
        )
        .unwrap();
        assert_eq!(r.stop_reason, StopReason::TargetReached);
        let mut d = r.merged.unwrap();
        d.normalize_log_count(256_f64.ln()).unwrap();
        let exact = Ising::<8>::exact_dos().unwrap();
        for (a, b) in d.dos().iter().zip(exact.dos()) {
            assert!((a - b).abs() < 0.3, "{a} != {b}");
        }
    }
}
