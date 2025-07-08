pub trait State {
    type Params;

    fn energy(&self, params: &mut Self::Params) -> f64;
}

pub struct Metropolis<S: State> {
    pub init_state: S,
    pub state: S,
    pub beta: f64,
    pub steps: usize,
}
