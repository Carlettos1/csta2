use csta::{Vec2f64, csta_derive::Randomizable};

fn main() {
    let ising = Ising::new(10, 10);
}

#[derive(Randomizable, Debug, Clone, Copy)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Randomizable, Debug)]
pub struct Ising {
    states: Vec<Spin>,
    w: usize,
    h: usize,
}

impl Spin {
    fn flip(&mut self) {
        match self {
            Spin::Down => *self = Spin::Up,
            Spin::Up => *self = Spin::Down,
        }
    }
}

impl Ising {
    fn new(w: usize, h: usize) -> Ising {
        Ising {
            states: Vec::with_capacity(w * h),
            w,
            h,
        }
    }

    fn flip(&mut self, x: usize, y: usize) {
        self.states[x + y * self.w].flip();
    }
}
