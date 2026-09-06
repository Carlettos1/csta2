use std::marker::PhantomData;

use csta_core::{
    vec2::{Vec2f32, Vec2f64},
    vec3::{Vec3f32, Vec3f64},
    vec4::{Vec4f32, Vec4f64},
};
/// Differences with v1:
/// distr is removed, as it was almost exclusively used with StandardUniform
/// and csta_derive gives options as to mul, div, add, sub
///
/// MCIter is removed and now MonteCarlo is an iterator
/// as it was always used as such
///  
use rand::{RngExt, rngs::ThreadRng};

pub trait Randomizable {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self;
}

#[derive(Debug)]
pub struct MonteCarlo<T: Randomizable, R: RngExt> {
    rng: R,
    phantom: PhantomData<T>,
}

impl Randomizable for f64 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Randomizable for f32 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Randomizable for Vec2f64 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec2f64(rng.random(), rng.random())
    }
}

impl Randomizable for Vec2f32 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec2f32(rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f64 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec3f64(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f32 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec3f32(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f64 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec4f64(rng.random(), rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f32 {
    fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
        Vec4f32(rng.random(), rng.random(), rng.random(), rng.random())
    }
}

///
/// if two elements are randomizable, a tuple of both elements also will be
macro_rules! randomize_tuple {
    ($($t:tt),*) => {
        impl<$($t,)+> Randomizable for ($($t,)+)
        where
            $($t: Randomizable,)+
        {
            fn sample<R: RngExt + ?Sized>(rng: &mut R) -> Self {
                ( $( <$t>::sample(rng), )+ )
            }
        }
    };
}

randomize_tuple! {A, B}
randomize_tuple! {A, B, C}
randomize_tuple! {A, B, C, D}
randomize_tuple! {A, B, C, D, E}
randomize_tuple! {A, B, C, D, E, F}
randomize_tuple! {A, B, C, D, E, F, G}
randomize_tuple! {A, B, C, D, E, F, G, H}

impl<T: Randomizable, R: RngExt> MonteCarlo<T, R> {
    #[warn(unused_must_use)]
    pub fn new(rng: R) -> Self {
        Self {
            rng,
            phantom: PhantomData,
        }
    }
}

impl<T: Randomizable, R: RngExt> Iterator for MonteCarlo<T, R> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        Some(<T>::sample(&mut self.rng))
    }
}

impl<T> Default for MonteCarlo<T, ThreadRng>
where
    T: Randomizable,
{
    fn default() -> Self {
        Self::new(rand::rng())
    }
}

/// Normal sample via Box–Muller. Zero scale returns the mean without RNG draws.
pub fn gaussian(
    rng: &mut (impl RngExt + ?Sized),
    mean: f64,
    scale: f64,
) -> Result<f64, &'static str> {
    if !mean.is_finite() || !scale.is_finite() || scale < 0.0 {
        return Err("invalid Gaussian parameters");
    }
    if scale == 0.0 {
        return Ok(mean);
    }
    let radius = (-2.0 * (1.0 - rng.random::<f64>()).ln()).sqrt();
    let z = radius * (std::f64::consts::TAU * rng.random::<f64>()).cos();
    let x = scale.mul_add(z, mean);
    if x.is_finite() {
        Ok(x)
    } else {
        Err("Gaussian sample overflow")
    }
}
/// Uniform direction on S², with exactly two random draws and no rejection loop.
pub fn isotropic_direction(rng: &mut (impl RngExt + ?Sized)) -> Vec3f64 {
    let z = 2.0 * rng.random::<f64>() - 1.0;
    let phi = std::f64::consts::TAU * rng.random::<f64>();
    let r = (1.0 - z * z).max(0.0).sqrt();
    Vec3f64(r * phi.cos(), r * phi.sin(), z)
}

/// Portable PCG64 with serializable state, for opt-in exact checkpoints.
#[cfg(feature = "checkpoint")]
pub use rand_pcg::Pcg64 as CheckpointRng;
#[cfg(test)]
mod tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};
    #[test]
    fn initialization() {
        let mut a = StdRng::seed_from_u64(4);
        let mut b = StdRng::seed_from_u64(4);
        assert_eq!(gaussian(&mut a, 3.0, 0.0), Ok(3.0));
        assert_eq!(a.random::<u64>(), b.random::<u64>());
        for scale in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(gaussian(&mut a, 0.0, scale).is_err());
        }
        let mut means = [0.0; 3];
        let mut squares = [0.0; 3];
        let mut cross = 0.0;
        for _ in 0..40_000 {
            let v = isotropic_direction(&mut a);
            assert!((v.len() - 1.0).abs() < 1e-14);
            let x: [f64; 3] = v.into();
            for i in 0..3 {
                means[i] += x[i];
                squares[i] += x[i] * x[i];
            }
            cross += x[0] * x[1];
        }
        for i in 0..3 {
            assert!(means[i].abs() / 40_000.0 < 0.015);
            assert!((squares[i] / 40_000.0 - 1.0 / 3.0).abs() < 0.015);
        }
        assert!(cross.abs() / 40_000.0 < 0.015);
        let v: Vec<_> = MonteCarlo::<(f64, Vec3f64), _>::new(StdRng::seed_from_u64(9))
            .take(100)
            .collect();
        let w: Vec<_> = MonteCarlo::<(f64, Vec3f64), _>::new(StdRng::seed_from_u64(9))
            .take(100)
            .collect();
        assert_eq!(v, w);
        for (x, v) in v {
            assert!((0.0..1.0).contains(&x));
            for x in <[f64; 3]>::from(v) {
                assert!((0.0..1.0).contains(&x));
            }
        }
    }
}
#[cfg(test)]
mod support_tests {
    use super::*;
    use rand::{SeedableRng, rngs::StdRng};
    #[test]
    fn tuple_arities_and_no_draws_for_empty_take() {
        #[derive(Debug, PartialEq)]
        struct One;
        impl Randomizable for One {
            fn sample<R: RngExt + ?Sized>(_: &mut R) -> Self {
                Self
            }
        }
        let mut r = StdRng::seed_from_u64(12);
        let mut reference = StdRng::seed_from_u64(12);
        let _: Vec<f64> = MonteCarlo::new(&mut r).take(0).collect();
        assert_eq!(r.random::<u64>(), reference.random::<u64>());
        assert_eq!(<(One, One)>::sample(&mut r), (One, One));
        assert_eq!(<(One, One, One)>::sample(&mut r), (One, One, One));
        assert_eq!(<(One, One, One, One)>::sample(&mut r), (One, One, One, One));
        assert_eq!(
            <(One, One, One, One, One)>::sample(&mut r),
            (One, One, One, One, One)
        );
        assert_eq!(
            <(One, One, One, One, One, One)>::sample(&mut r),
            (One, One, One, One, One, One)
        );
        assert_eq!(
            <(One, One, One, One, One, One, One)>::sample(&mut r),
            (One, One, One, One, One, One, One)
        );
        assert_eq!(
            <(One, One, One, One, One, One, One, One)>::sample(&mut r),
            (One, One, One, One, One, One, One, One)
        );
    }
    #[test]
    fn uniform_and_gaussian_moments() {
        let mut rng = StdRng::seed_from_u64(31);
        let mut u = 0.0;
        let mut u2 = 0.0;
        let mut g = 0.0;
        let mut g2 = 0.0;
        for _ in 0..50_000 {
            let x = f64::sample(&mut rng);
            let y = gaussian(&mut rng, 2.0, 3.0).unwrap();
            u += x;
            u2 += x * x;
            g += y;
            g2 += y * y;
        }
        assert!((u / 50_000.0 - 0.5).abs() < 0.01);
        assert!((u2 / 50_000.0 - 1.0 / 3.0).abs() < 0.01);
        assert!((g / 50_000.0 - 2.0).abs() < 0.08);
        assert!((g2 / 50_000.0 - 13.0).abs() < 0.35);
    }
}
