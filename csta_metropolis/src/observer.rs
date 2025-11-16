use crate::State;

pub trait Observer<S: State> {
    /// for energy, this is f64
    type Observation;

    /// Measure the state, returns an observation.
    fn measure(state: &S, params: &S::Params) -> Self::Observation;

    /// every nth step will measure.
    fn every() -> usize;

    /// only after nth step will start to measure.
    fn after() -> usize;
}

/// Not different from above, just has &self so it can be built into a dyn Trait
pub trait DynObserver<S: State> {
    type Observation;
    fn measure(&self, state: &S, params: &S::Params) -> Self::Observation;
    fn every(&self) -> usize;
    fn after(&self) -> usize;
}

/// All Observers are DynObservers
impl<O, S: State> DynObserver<S> for O
where
    O: Observer<S>,
{
    type Observation = O::Observation;

    fn after(&self) -> usize {
        O::after()
    }

    fn every(&self) -> usize {
        O::every()
    }

    fn measure(&self, state: &S, params: &<S as State>::Params) -> Self::Observation {
        O::measure(state, params)
    }
}
