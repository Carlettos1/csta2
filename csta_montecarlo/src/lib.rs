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
use rand::{Rng, rngs::ThreadRng};

pub trait Randomizable {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

#[derive(Debug)]
pub struct MonteCarlo<T: Randomizable, R: Rng> {
    rng: R,
    phantom: PhantomData<T>,
}

impl Randomizable for f64 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Randomizable for f32 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        rng.random()
    }
}

impl Randomizable for Vec2f64 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Vec2f64(rng.random(), rng.random())
    }
}

impl Randomizable for Vec2f32 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Vec2f32(rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f64 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Vec3f64(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f32 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Vec3f32(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f64 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        Vec4f64(rng.random(), rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f32 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
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
            fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
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

impl<T: Randomizable, R: Rng> MonteCarlo<T, R> {
    #[warn(unused_must_use)]
    pub fn new(rng: R) -> Self {
        Self {
            rng,
            phantom: PhantomData,
        }
    }
}

impl<T: Randomizable, R: Rng> Iterator for MonteCarlo<T, R> {
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
