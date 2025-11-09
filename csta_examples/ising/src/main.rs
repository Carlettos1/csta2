use csta::csta_derive::Randomizable;
use rand::Rng;

fn main() {
    let ising = Ising::new(10, 10);
    println!("{:?} {:?} {:?}", ising.states, ising.w, ising.h);
}

#[derive(Randomizable, Debug, Clone, Copy)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Randomizable, Debug)]
pub struct Ising {
    #[csta(default = 10)]
    w: usize,
    #[csta(range(5..10))]
    h: usize,
    #[csta(len(w * h))]
    states: Vec<Spin>,
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
