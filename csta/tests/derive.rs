use csta::{Randomizable, csta_derive::Randomizable as Derive};
use rand::{RngExt, SeedableRng, rngs::StdRng};
#[derive(Derive, Debug, PartialEq)]
struct Unit;
#[derive(Derive)]
struct Generic<T> {
    x: T,
}
#[derive(Derive)]
struct DefaultOnly<T> {
    #[csta(default)]
    x: T,
}
#[derive(Derive)]
struct Tuple(#[csta(range(2..=2))] u32, #[csta(default)] f64);
#[derive(Derive)]
struct Fields {
    #[csta(default=n+1)]
    next: usize,
    #[csta(default = 3)]
    n: usize,
    #[csta(len(n))]
    values: Vec<f64>,
    #[csta(after(raw*2.0))]
    raw: f64,
    #[csta(default=raw)]
    final_value: f64,
    #[csta(mul = 2.0, add = 1.0)]
    scaled: f64,
}
#[derive(Derive)]
enum Weighted {
    #[csta(weight = 0)]
    Never,
    #[csta(weight = 1)]
    A,
    #[csta(weight = 3)]
    B,
}
#[derive(Derive)]
enum Shapes {
    Unit,
    Tuple(f64),
    Named { x: f64 },
}
#[test]
fn macro_sampling_and_dependencies() {
    let _ = Weighted::Never;
    let mut r = StdRng::seed_from_u64(4);
    assert_eq!(Unit::sample(&mut r), Unit);
    assert!((0.0..1.0).contains(&Generic::<f64>::sample(&mut r).x));
    assert!(DefaultOnly::<Vec<i32>>::sample(&mut r).x.is_empty());
    let Tuple(a, b) = Tuple::sample(&mut r);
    assert_eq!((a, b), (2, 0.0));
    let mut reference = StdRng::seed_from_u64(7);
    let mut rng = StdRng::seed_from_u64(7);
    let raw = reference.random::<f64>();
    let values: Vec<f64> = (0..3).map(|_| reference.random()).collect();
    let scaled = reference.random::<f64>() * 2.0 + 1.0;
    let f = Fields::sample(&mut rng);
    assert_eq!(f.n, 3);
    assert_eq!(f.next, 4);
    assert_eq!(f.values, values);
    assert_eq!(f.raw, raw * 2.0);
    assert_eq!(f.final_value, f.raw);
    assert_eq!(f.scaled, scaled);
    let mut a = 0;
    for _ in 0..20_000 {
        match Weighted::sample(&mut r) {
            Weighted::Never => panic!("zero weight sampled"),
            Weighted::A => a += 1,
            Weighted::B => {}
        }
    }
    assert!((a as f64 / 20_000.0 - 0.25).abs() < 0.02);
    let mut seen = [false; 3];
    for _ in 0..100 {
        match Shapes::sample(&mut r) {
            Shapes::Unit => seen[0] = true,
            Shapes::Tuple(x) => {
                assert!((0.0..1.0).contains(&x));
                seen[1] = true
            }
            Shapes::Named { x } => {
                assert!((0.0..1.0).contains(&x));
                seen[2] = true
            }
        }
    }
    assert!(seen.into_iter().all(|x| x));
}
