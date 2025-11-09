//!
//! This are just compilation and macro expansion tests for csta::csta_derive
//!
#![allow(dead_code)]
#![allow(clippy::disallowed_names)]

use csta::{Vec3f64, csta_derive::Randomizable};

#[derive(Randomizable)]
enum Spin {
    Up,
    Down,
}

#[derive(Randomizable)]
struct Ising {
    #[csta(default = 10)]
    w: usize,
    #[csta(range(5..10))]
    h: usize,
    #[csta(len(w * h))]
    states: Vec<Spin>,
}

#[derive(Randomizable)]
enum A {
    Part1,
    Part2,
    Part3,
}

#[derive(Randomizable)]
struct Nothing;

/// Energy: E
/// Kinetic energy: T
/// Potential energy: V
/// Velocity: v
/// mass: m
/// pos: r
///
/// E = T + V
///
/// T = m*|v|²/2, so
/// |v|² = 2*T/m
#[derive(Randomizable)]
struct Particle {
    #[csta(range(5.0..10.0))]
    mass: f64,
    #[csta(range(5.0..15.0))]
    kinetic_energy: f64,
    #[csta(default = 20.0)]
    energy: f64,
    // in a 10 units box
    #[csta(mul = 10.0)]
    pos: Vec3f64,
    // should respect kinetic energy and mass numbers
    #[csta(after(vel.normalize() * (kinetic_energy * 2.0_f64 / mass).sqrt()))]
    vel: Vec3f64,
    #[csta(default)]
    force: Vec3f64,
}

#[derive(Randomizable)]
struct After {
    #[csta(range(0.0..11.0))]
    padding: f64,
    #[csta(after(val1 + Vec3f64(0.0, 11.0, padding)))]
    val1: Vec3f64,
    #[csta(after(val2.normalize() * val1.len()))]
    val2: Vec3f64,
    #[csta(after(func1(val3 + val2)))]
    val3: Vec3f64,
}

fn func1(vec: Vec3f64) -> Vec3f64 {
    let mut foo = 1.0_f64;
    let mut x = vec.x();
    if x < 0.5 {
        x = x.sin();
    } else {
        foo -= 0.15;
        x -= x;
    }

    let mut y = vec.y();
    let mut z = vec.z();

    if ((y * z) * 1000.0) as i64 % 2 == 0 {
        y -= 1.0;
        foo = foo.powf(1.2);
    } else {
        std::mem::swap(&mut y, &mut z);
    }

    Vec3f64(x, y, z) * foo
}

#[derive(Randomizable)]
enum EEnum {
    #[csta(weight = 1.0)]
    Case1(After),
    #[csta(weight = 1.0)]
    Case2 {
        particle: Particle,
        #[csta(after(foo.normalize()))]
        foo: Vec3f64,
    },
    #[csta(weight = 1.0)]
    Case3(A),
}
