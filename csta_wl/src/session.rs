//! Incremental WL/SAMC sessions, retaining RNG and adaptation state between calls.
use crate::*;
#[cfg_attr(feature = "checkpoint", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "checkpoint",
    serde(bound(
        serialize = "S: serde::Serialize, S::Params: serde::Serialize, R: serde::Serialize",
        deserialize = "S: serde::Deserialize<'de>, S::Params: serde::Deserialize<'de>, R: serde::Deserialize<'de>"
    ))
)]
pub struct Session<S: State, R> {
    walker: Walker<S, R>,
}
impl<S: State, R: RngExt> Session<S, R>
where
    S::Params: Clone,
{
    pub fn new(
        state: S,
        params: S::Params,
        rng: R,
        data: RawWangLandauData,
        config: Config,
    ) -> Result<Self> {
        Ok(Self {
            walker: Walker::new(state, params, rng, data, config)?,
        })
    }
    /// Resume after cancellation at the next proposal boundary. Total budget and
    /// target remain fixed. `moves=0` leaves the snapshot unchanged.
    pub fn advance(&mut self, moves: u64, cancel: Option<&AtomicBool>) -> Result<&Diagnostics> {
        self.advance_observed(moves, cancel, |_, _| Ok(()))
    }
    /// Observe retained production states, including rejected and self-loop moves.
    /// The callback gets the cached energy and can record conditional moments.
    /// Warm-up is excluded. Callback errors stop after the move has committed.
    pub fn advance_observed(
        &mut self,
        moves: u64,
        cancel: Option<&AtomicBool>,
        mut observe: impl FnMut(&S, f64) -> Result<()>,
    ) -> Result<&Diagnostics> {
        if moves == 0 {
            return Ok(&self.walker.diagnostics);
        }
        if self.walker.diagnostics.stop_reason == Some(StopReason::Cancelled) {
            self.walker.diagnostics.stop_reason = None;
        }
        for _ in 0..moves {
            let before = self.walker.diagnostics.sampling_steps;
            self.walker.advance(cancel)?;
            if self.walker.diagnostics.sampling_steps > before {
                observe(&self.walker.state, self.walker.energy)?;
            }
            if self.walker.diagnostics.stop_reason.is_some() {
                break;
            }
        }
        Ok(&self.walker.diagnostics)
    }
    pub fn data(&self) -> &RawWangLandauData {
        &self.walker.data
    }
    pub fn diagnostics(&self) -> &Diagnostics {
        &self.walker.diagnostics
    }
    pub fn state(&self) -> &S {
        &self.walker.state
    }
    /// A deliberately unfinished session has an explicit partial stop reason.
    pub fn finish(mut self) -> RunResult<S> {
        if self.walker.diagnostics.stop_reason.is_none() {
            self.walker.diagnostics.stop_reason = Some(StopReason::Cancelled);
        }
        self.walker.finish()
    }
    #[cfg(feature = "checkpoint")]
    fn validate(&self) -> Result<()> {
        validate_walker(&self.walker, false)
    }
}
#[cfg(feature = "checkpoint")]
pub(crate) fn validate_walker<S: State, R: RngExt>(
    w: &Walker<S, R>,
    allow_mixing: bool,
) -> Result<()>
where
    S::Params: Clone,
{
    w.data.grid.validate()?;
    RawWangLandauData::with_support(
        w.data.grid.clone(),
        w.data.accessible.clone(),
        w.data.dos.clone(),
    )?;
    let n = w.data.grid.len();
    let checked_sum = |v: &[u64]| {
        v.iter()
            .try_fold(0u64, |a, b| a.checked_add(*b))
            .ok_or(Error::CounterOverflow)
    };
    if w.data.bins.len() != n
        || w.data.lifetime_bins.len() != n
        || checked_sum(&w.data.bins)? != w.data.visits
        || checked_sum(&w.data.lifetime_bins)? != w.data.total_visits
        || w.data.total_visits != w.diagnostics.proposals
        || w.data
            .bins
            .iter()
            .zip(&w.data.lifetime_bins)
            .any(|(a, b)| a > b)
        || w.data
            .lifetime_bins
            .iter()
            .zip(&w.data.accessible)
            .any(|(n, a)| !a && *n != 0)
    {
        return Err(Error::Invalid("invalid snapshot counters"));
    }
    let (minimum, t0, t1) = w.config.validate(w.data.active_bins())?;
    if minimum != w.minimum_visits
        || t0 != w.t0
        || t1 != w.t1
        || !w.ln_f.is_finite()
        || w.ln_f < 0.0
        || w.diagnostics.completed_stages > w.config.preliminary_stages
        || w.diagnostics.accepted > w.diagnostics.proposals
        || w.diagnostics.proposals > w.config.max_steps
        || w.diagnostics.sampling_steps > w.config.sampling_steps
        || w.diagnostics.sampling_steps > w.diagnostics.proposals
        || (!allow_mixing && w.diagnostics.mixing_steps != 0)
        || !w.diagnostics.last_update.is_finite()
    {
        return Err(Error::Invalid("invalid snapshot schedule"));
    }
    let expected = w.config.initial_ln_f * 0.5_f64.powi(w.diagnostics.completed_stages as i32);
    if w.ln_f != expected
        || (w.diagnostics.phase == Phase::Preliminary && w.diagnostics.sampling_steps != 0)
        || (w.diagnostics.phase == Phase::Sampling
            && w.diagnostics.completed_stages != w.config.preliminary_stages)
    {
        return Err(Error::Invalid("inconsistent snapshot phase"));
    }
    if w.diagnostics.stop_reason == Some(StopReason::TargetReached)
        && w.diagnostics.sampling_steps != w.config.sampling_steps
    {
        return Err(Error::Invalid("false completed snapshot"));
    }
    if (w.diagnostics.stop_reason == Some(StopReason::BudgetExhausted)
        && w.diagnostics.proposals != w.config.max_steps)
        || (w.diagnostics.stop_reason.is_none() && w.diagnostics.proposals == w.config.max_steps)
        || (w.diagnostics.phase == Phase::Sampling && w.data.visits != w.diagnostics.sampling_steps)
        || (w.diagnostics.phase == Phase::Preliminary
            && w.diagnostics.completed_stages == w.config.preliminary_stages)
    {
        return Err(Error::Invalid(
            "inconsistent snapshot completion or visit clock",
        ));
    }
    if !w.state.valid_state(&w.params) {
        return Err(Error::Invalid("invalid snapshot model"));
    }
    let energy = w.state.energy(&mut w.params.clone());
    if energy != w.energy || w.data.energy_to_bin(energy)? != Some(w.bin) {
        return Err(Error::Invalid("snapshot state/grid/cache mismatch"));
    }
    Ok(())
}
#[cfg(feature = "checkpoint")]
impl<S: State, R: RngExt> Session<S, R>
where
    S: serde::Serialize + serde::de::DeserializeOwned,
    S::Params: Clone + serde::Serialize + serde::de::DeserializeOwned,
    R: serde::Serialize + serde::de::DeserializeOwned,
{
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
