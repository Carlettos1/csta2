//! A repeatable local microbenchmark, not a timing-sensitive unit test.
//! cargo run --release --example profile_energy
use csta_wl::{State, models::Ising};
use std::{hint::black_box, time::Instant};

fn main() {
    let state = Ising::<4096>::new([1; 4096]).unwrap();
    let iterations = 20_000;
    let start = Instant::now();
    for _ in 0..iterations {
        black_box(black_box(&state).recompute_energy());
    }
    let full = start.elapsed();
    let start = Instant::now();
    for _ in 0..iterations {
        black_box(black_box(&state).energy(&mut ()));
    }
    let cached = start.elapsed();
    println!("{iterations} energy queries on 4096 spins: full={full:?}, cached={cached:?}");
    assert_eq!(state.energy(&mut ()), state.recompute_energy());
}
