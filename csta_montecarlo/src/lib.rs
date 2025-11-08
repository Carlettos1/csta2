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
    type Conf;
    fn sample<R: Rng + ?Sized>(rng: &mut R, config: &Self::Conf) -> Self;
}

#[derive(Debug)]
pub struct MonteCarlo<T: Randomizable, R: Rng> {
    rng: R,
    config: T::Conf,
}

impl Randomizable for f64 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        rng.random()
    }
}

impl Randomizable for f32 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        rng.random()
    }
}

impl Randomizable for Vec2f64 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec2f64(rng.random(), rng.random())
    }
}

impl Randomizable for Vec2f32 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec2f32(rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f64 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec3f64(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec3f32 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec3f32(rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f64 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec4f64(rng.random(), rng.random(), rng.random(), rng.random())
    }
}

impl Randomizable for Vec4f32 {
    type Conf = ();
    fn sample<R: Rng + ?Sized>(rng: &mut R, _: &Self::Conf) -> Self {
        Vec4f32(rng.random(), rng.random(), rng.random(), rng.random())
    }
}

impl<T> Randomizable for Vec<T>
where
    T: Randomizable,
{
    type Conf = (usize, T::Conf);
    fn sample<R: Rng + ?Sized>(rng: &mut R, config: &Self::Conf) -> Self {
        (0..config.0).map(|_| T::sample(rng, &config.1)).collect()
    }
}

///
/// if two elements are randomizable, a tuple of both elements also will be
macro_rules! randomize_tuple {
    ($($t:tt),* | $($n:tt),*) => {
        impl<$($t,)+> Randomizable for ($($t,)+)
        where
            $($t: Randomizable,)+
        {
            type Conf = ($($t::Conf,)+);
            fn sample<R: Rng + ?Sized>(rng: &mut R, config: &Self::Conf) -> Self {
                ( $( <$t>::sample(rng, &config.$n), )+ )
            }
        }
    };
}

randomize_tuple! {A, B | 0, 1}
randomize_tuple! {A, B, C | 0, 1, 2}
randomize_tuple! {A, B, C, D | 0, 1, 2, 3}
randomize_tuple! {A, B, C, D, E | 0, 1, 2, 3, 4}
randomize_tuple! {A, B, C, D, E, F | 0, 1, 2, 3, 4, 5}
randomize_tuple! {A, B, C, D, E, F, G | 0, 1, 2, 3, 4, 5, 6}
randomize_tuple! {A, B, C, D, E, F, G, H | 0, 1, 2, 3, 4, 5, 6, 7}

impl<T: Randomizable, R: Rng> MonteCarlo<T, R> {
    #[warn(unused_must_use)]
    pub fn new(rng: R, config: T::Conf) -> Self {
        Self { rng, config }
    }
}

impl<T: Randomizable, R: Rng> Iterator for MonteCarlo<T, R> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        Some(<T>::sample(&mut self.rng, &self.config))
    }
}

impl<T> Default for MonteCarlo<T, ThreadRng>
where
    T: Randomizable,
    T::Conf: Default,
{
    fn default() -> Self {
        Self::with_config(T::Conf::default())
    }
}

impl<T: Randomizable> MonteCarlo<T, ThreadRng> {
    pub fn with_config(config: T::Conf) -> Self {
        Self {
            rng: rand::rng(),
            config,
        }
    }
}
