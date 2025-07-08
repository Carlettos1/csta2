use csta::Vec2f64;
use serde_derive::{Deserialize, Serialize};

fn main() {
    let vec = Vec2f64::new(1.0, 1.0);
}

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
pub struct A {
    vec: Vec2f64,
}

pub enum State {
    Up,
    Down,
}

pub struct System {
    states: Vec<Vec<State>>,
}
