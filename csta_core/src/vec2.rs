//! A 2D vector representation.
//! This module provides a simple 2D vector type with basic operations.
//!

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec2f64(pub f64, pub f64);

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec2f32(pub f32, pub f32);

impl Vec2f64 {
    pub fn new(x: f64, y: f64) -> Self {
        Vec2f64(x, y)
    }

    pub fn x(&self) -> f64 {
        self.0
    }

    pub fn y(&self) -> f64 {
        self.1
    }

    pub fn len(&self) -> f64 {
        self.len_squared().sqrt()
    }

    pub fn len_squared(&self) -> f64 {
        self.dot(self)
    }

    pub fn normalize(&self) -> Self {
        let len = self.len();
        if len == 0.0 {
            Vec2f64(0.0, 0.0)
        } else {
            Vec2f64(self.0 / len, self.1 / len)
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.0 * other.0 + self.1 * other.1
    }
}

impl Vec2f32 {
    pub fn new(x: f32, y: f32) -> Self {
        Vec2f32(x, y)
    }

    pub fn x(&self) -> f32 {
        self.0
    }

    pub fn y(&self) -> f32 {
        self.1
    }

    pub fn len(&self) -> f32 {
        self.len_squared().sqrt()
    }

    pub fn len_squared(&self) -> f32 {
        self.dot(self)
    }

    pub fn normalize(&self) -> Self {
        let len = self.len();
        if len == 0.0 {
            Vec2f32(0.0, 0.0)
        } else {
            Vec2f32(self.0 / len, self.1 / len)
        }
    }

    pub fn dot(&self, other: &Self) -> f32 {
        self.0 * other.0 + self.1 * other.1
    }
}

macro_rules! impl_from_tuple {
    ($vec:ident, $tuple:ty) => {
        impl From<$tuple> for $vec {
            fn from(tuple: $tuple) -> Self {
                $vec(tuple.0, tuple.1)
            }
        }
    };
}

macro_rules! impl_tuple_from {
    ($vec:ty, $tuple:ty) => {
        impl From<$vec> for $tuple {
            fn from(vec: $vec) -> Self {
                (vec.0, vec.1)
            }
        }
    };
}

macro_rules! impl_from_array {
    ($vec:ident, $array:ty) => {
        impl From<$array> for $vec {
            fn from(arr: $array) -> Self {
                $vec(arr[0], arr[1])
            }
        }
    };
}

macro_rules! impl_array_from {
    ($vec:ty, $array:ty) => {
        impl From<$vec> for $array {
            fn from(vec: $vec) -> Self {
                [vec.0, vec.1]
            }
        }
    };
}

macro_rules! impl_add {
    ($vec:ident, $vec1:ty, $vec2:ty) => {
        impl std::ops::Add<$vec2> for $vec1 {
            type Output = $vec;

            fn add(self, other: $vec2) -> Self::Output {
                $vec(self.0 + other.0, self.1 + other.1)
            }
        }
    };
}

macro_rules! impl_add_assign {
    ($vec:ident, $other:ty) => {
        impl std::ops::AddAssign<$other> for $vec {
            fn add_assign(&mut self, other: $other) {
                self.0 += other.0;
                self.1 += other.1;
            }
        }
    };
}

macro_rules! impl_sub {
    ($vec:ident, $vec1:ty, $vec2:ty) => {
        impl std::ops::Sub<$vec2> for $vec1 {
            type Output = $vec;

            fn sub(self, other: $vec2) -> Self::Output {
                $vec(self.0 - other.0, self.1 - other.1)
            }
        }
    };
}

macro_rules! impl_sub_assign {
    ($vec:ident, $other:ty) => {
        impl std::ops::SubAssign<$other> for $vec {
            fn sub_assign(&mut self, other: $other) {
                self.0 -= other.0;
                self.1 -= other.1;
            }
        }
    };
}

macro_rules! impl_mul {
    ($vec:ident, $vec1:ty, $float:ty) => {
        impl std::ops::Mul<$float> for $vec1 {
            type Output = $vec;

            fn mul(self, scalar: $float) -> Self::Output {
                $vec(self.0 * scalar, self.1 * scalar)
            }
        }
    };
}

macro_rules! impl_mul_assign {
    ($vec:ident, $float:ident) => {
        impl std::ops::MulAssign<$float> for $vec {
            fn mul_assign(&mut self, scalar: $float) {
                self.0 *= scalar;
                self.1 *= scalar;
            }
        }
    };
}

macro_rules! impl_div {
    ($vec:ident, $vec1:ty, $float:ty) => {
        impl std::ops::Div<$float> for $vec1 {
            type Output = $vec;

            fn div(self, scalar: $float) -> Self::Output {
                if scalar.is_normal() {
                    $vec(self.0 / scalar, self.1 / scalar)
                } else {
                    $vec(self.0 / scalar, self.1 / scalar)
                }
            }
        }
    };
}

macro_rules! impl_div_assign {
    ($vec:ident, $float:ident) => {
        impl std::ops::DivAssign<$float> for $vec {
            fn div_assign(&mut self, scalar: $float) {
                if scalar.is_normal() {
                    self.0 /= scalar;
                    self.1 /= scalar;
                } else {
                    self.0 = 0.0;
                    self.1 = 0.0;
                }
            }
        }
    };
}

macro_rules! impl_neg {
    ($vec:ident, $vec1:ty) => {
        impl std::ops::Neg for $vec1 {
            type Output = $vec;

            fn neg(self) -> Self::Output {
                $vec(-self.0, -self.1)
            }
        }
    };
}

impl_from_tuple!(Vec2f64, (f64, f64));
impl_from_tuple!(Vec2f64, &(f64, f64));

impl_tuple_from!(Vec2f64, (f64, f64));
impl_tuple_from!(&Vec2f64, (f64, f64));

impl_from_array!(Vec2f64, [f64; 2]);
impl_from_array!(Vec2f64, &[f64; 2]);

impl_array_from!(Vec2f64, [f64; 2]);
impl_array_from!(&Vec2f64, [f64; 2]);

impl_add!(Vec2f64, Vec2f64, Vec2f64);
impl_add!(Vec2f64, Vec2f64, &Vec2f64);
impl_add!(Vec2f64, &Vec2f64, Vec2f64);
impl_add!(Vec2f64, &Vec2f64, &Vec2f64);

impl_add_assign!(Vec2f64, Vec2f64);
impl_add_assign!(Vec2f64, &Vec2f64);

impl_sub!(Vec2f64, Vec2f64, Vec2f64);
impl_sub!(Vec2f64, Vec2f64, &Vec2f64);
impl_sub!(Vec2f64, &Vec2f64, Vec2f64);
impl_sub!(Vec2f64, &Vec2f64, &Vec2f64);

impl_sub_assign!(Vec2f64, Vec2f64);
impl_sub_assign!(Vec2f64, &Vec2f64);

impl_mul!(Vec2f64, Vec2f64, f64);
impl_mul!(Vec2f64, &Vec2f64, f64);
impl_mul_assign!(Vec2f64, f64);

impl_div!(Vec2f64, Vec2f64, f64);
impl_div!(Vec2f64, &Vec2f64, f64);

impl_div_assign!(Vec2f64, f64);
impl_neg!(Vec2f64, Vec2f64);
impl_neg!(Vec2f64, &Vec2f64);

impl_from_tuple!(Vec2f32, (f32, f32));
impl_from_tuple!(Vec2f32, &(f32, f32));

impl_tuple_from!(Vec2f32, (f32, f32));
impl_tuple_from!(&Vec2f32, (f32, f32));

impl_from_array!(Vec2f32, [f32; 2]);
impl_from_array!(Vec2f32, &[f32; 2]);

impl_array_from!(Vec2f32, [f32; 2]);
impl_array_from!(&Vec2f32, [f32; 2]);

impl_add!(Vec2f32, Vec2f32, Vec2f32);
impl_add!(Vec2f32, Vec2f32, &Vec2f32);
impl_add!(Vec2f32, &Vec2f32, Vec2f32);
impl_add!(Vec2f32, &Vec2f32, &Vec2f32);

impl_add_assign!(Vec2f32, Vec2f32);
impl_add_assign!(Vec2f32, &Vec2f32);

impl_sub!(Vec2f32, Vec2f32, Vec2f32);
impl_sub!(Vec2f32, Vec2f32, &Vec2f32);
impl_sub!(Vec2f32, &Vec2f32, Vec2f32);
impl_sub!(Vec2f32, &Vec2f32, &Vec2f32);

impl_sub_assign!(Vec2f32, Vec2f32);
impl_sub_assign!(Vec2f32, &Vec2f32);

impl_mul!(Vec2f32, Vec2f32, f32);
impl_mul!(Vec2f32, &Vec2f32, f32);
impl_mul_assign!(Vec2f32, f32);

impl_div!(Vec2f32, Vec2f32, f32);
impl_div!(Vec2f32, &Vec2f32, f32);
impl_div_assign!(Vec2f32, f32);

impl_neg!(Vec2f32, Vec2f32);
impl_neg!(Vec2f32, &Vec2f32);
