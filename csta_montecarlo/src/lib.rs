use std::marker::PhantomData;

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
    _phantom: PhantomData<T>,
}

impl Randomizable for f64 {
    fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
        rng.random()
    }
}

///
/// if two elements are randomizable, a tuple of both elements also will be
macro_rules! randomize_tuple {
    ($($t:tt)+) => {
        impl<$($t: Randomizable,)+> Randomizable for ($($t,)+) {
            fn sample<R: Rng + ?Sized>(rng: &mut R) -> Self {
                ( $( <$t>::sample(rng), )+ )
            }
        }
    };
}

randomize_tuple! {A B}
randomize_tuple! {A B C}
randomize_tuple! {A B C E}
randomize_tuple! {A B C E F}
randomize_tuple! {A B C E F G}
randomize_tuple! {A B C E F G H}
randomize_tuple! {A B C E F G H I}

impl<T: Randomizable, R: Rng> MonteCarlo<T, R> {
    #[warn(unused_must_use)]
    pub fn new(rng: R) -> Self {
        Self {
            rng,
            _phantom: PhantomData,
        }
    }
}

impl<T: Randomizable, R: Rng> Iterator for MonteCarlo<T, R> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        Some(<T>::sample(&mut self.rng))
    }
}

impl<T: Randomizable> Default for MonteCarlo<T, ThreadRng> {
    fn default() -> Self {
        Self::new(rand::rng())
    }
}
