//! A 3D vector representation.
//! This module provides a simple 3D vector type with basic operations.
//!

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3f64(pub f64, pub f64, pub f64);

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3f32(pub f32, pub f32, pub f32);

impl Vec3f64 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3f64(x, y, z)
    }

    pub fn x(&self) -> f64 {
        self.0
    }

    pub fn y(&self) -> f64 {
        self.1
    }

    pub fn z(&self) -> f64 {
        self.2
    }

    pub fn len(&self) -> f64 {
        self.len_squared().sqrt()
    }

    pub fn len_squared(&self) -> f64 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }

    pub fn normalize(&self) -> Self {
        let len = self.len();
        if len == 0.0 {
            Vec3f64(0.0, 0.0, 0.0)
        } else {
            Vec3f64(self.0 / len, self.1 / len, self.2 / len)
        }
    }
}

impl Vec3f32 {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Vec3f32(x, y, z)
    }

    pub fn x(&self) -> f32 {
        self.0
    }

    pub fn y(&self) -> f32 {
        self.1
    }

    pub fn z(&self) -> f32 {
        self.2
    }

    pub fn len(&self) -> f32 {
        self.len_squared().sqrt()
    }

    pub fn len_squared(&self) -> f32 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }

    pub fn normalize(&self) -> Self {
        let len = self.len();
        if len == 0.0 {
            Vec3f32(0.0, 0.0, 0.0)
        } else {
            Vec3f32(self.0 / len, self.1 / len, self.2 / len)
        }
    }
}

macro_rules! impl_from_tuple {
    ($vec:ident, $tuple:ty) => {
        impl From<$tuple> for $vec {
            fn from(tuple: $tuple) -> Self {
                $vec(tuple.0, tuple.1, tuple.2)
            }
        }
    };
}

macro_rules! impl_tuple_from {
    ($vec:ty, $tuple:ty) => {
        impl From<$vec> for $tuple {
            fn from(vec: $vec) -> Self {
                (vec.0, vec.1, vec.2)
            }
        }
    };
}

macro_rules! impl_from_array {
    ($vec:ident, $array:ty) => {
        impl From<$array> for $vec {
            fn from(arr: $array) -> Self {
                $vec(arr[0], arr[1], arr[2])
            }
        }
    };
}

macro_rules! impl_array_from {
    ($vec:ty, $array:ty) => {
        impl From<$vec> for $array {
            fn from(vec: $vec) -> Self {
                [vec.0, vec.1, vec.2]
            }
        }
    };
}

macro_rules! impl_add {
    ($vec:ident, $vec1:ty, $Vec3:ty) => {
        impl std::ops::Add<$Vec3> for $vec1 {
            type Output = $vec;

            fn add(self, other: $Vec3) -> Self::Output {
                $vec(self.0 + other.0, self.1 + other.1, self.2 + other.2)
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
                self.2 += other.2;
            }
        }
    };
}

macro_rules! impl_sub {
    ($vec:ident, $vec1:ty, $Vec3:ty) => {
        impl std::ops::Sub<$Vec3> for $vec1 {
            type Output = $vec;

            fn sub(self, other: $Vec3) -> Self::Output {
                $vec(self.0 - other.0, self.1 - other.1, self.2 - other.2)
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
                self.2 -= other.2;
            }
        }
    };
}

macro_rules! impl_mul {
    ($vec:ident, $vec1:ty, $float:ty) => {
        impl std::ops::Mul<$float> for $vec1 {
            type Output = $vec;

            fn mul(self, scalar: $float) -> Self::Output {
                $vec(self.0 * scalar, self.1 * scalar, self.2 * scalar)
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
                self.2 *= scalar;
            }
        }
    };
}

macro_rules! impl_div {
    ($vec:ident, $vec1:ty, $float:ty) => {
        impl std::ops::Div<$float> for $vec1 {
            type Output = $vec;

            fn div(self, scalar: $float) -> Self::Output {
                $vec(self.0 / scalar, self.1 / scalar, self.2 / scalar)
            }
        }
    };
}

macro_rules! impl_div_assign {
    ($vec:ident, $float:ident) => {
        impl std::ops::DivAssign<$float> for $vec {
            fn div_assign(&mut self, scalar: $float) {
                self.0 /= scalar;
                self.1 /= scalar;
                self.2 /= scalar;
            }
        }
    };
}

macro_rules! impl_neg {
    ($vec:ident, $vec1:ty) => {
        impl std::ops::Neg for $vec1 {
            type Output = $vec;

            fn neg(self) -> Self::Output {
                $vec(-self.0, -self.1, -self.2)
            }
        }
    };
}

impl_from_tuple!(Vec3f64, (f64, f64, f64));
impl_from_tuple!(Vec3f64, &(f64, f64, f64));

impl_tuple_from!(Vec3f64, (f64, f64, f64));
impl_tuple_from!(&Vec3f64, (f64, f64, f64));

impl_from_array!(Vec3f64, [f64; 3]);
impl_from_array!(Vec3f64, &[f64; 3]);

impl_array_from!(Vec3f64, [f64; 3]);
impl_array_from!(&Vec3f64, [f64; 3]);

impl_add!(Vec3f64, Vec3f64, Vec3f64);
impl_add!(Vec3f64, Vec3f64, &Vec3f64);
impl_add!(Vec3f64, &Vec3f64, Vec3f64);
impl_add!(Vec3f64, &Vec3f64, &Vec3f64);

impl_add_assign!(Vec3f64, Vec3f64);
impl_add_assign!(Vec3f64, &Vec3f64);

impl_sub!(Vec3f64, Vec3f64, Vec3f64);
impl_sub!(Vec3f64, Vec3f64, &Vec3f64);
impl_sub!(Vec3f64, &Vec3f64, Vec3f64);
impl_sub!(Vec3f64, &Vec3f64, &Vec3f64);

impl_sub_assign!(Vec3f64, Vec3f64);
impl_sub_assign!(Vec3f64, &Vec3f64);

impl_mul!(Vec3f64, Vec3f64, f64);
impl_mul!(Vec3f64, &Vec3f64, f64);
impl_mul_assign!(Vec3f64, f64);

impl_div!(Vec3f64, Vec3f64, f64);
impl_div!(Vec3f64, &Vec3f64, f64);

impl_div_assign!(Vec3f64, f64);
impl_neg!(Vec3f64, Vec3f64);
impl_neg!(Vec3f64, &Vec3f64);

impl_from_tuple!(Vec3f32, (f32, f32, f32));
impl_from_tuple!(Vec3f32, &(f32, f32, f32));

impl_tuple_from!(Vec3f32, (f32, f32, f32));
impl_tuple_from!(&Vec3f32, (f32, f32, f32));

impl_from_array!(Vec3f32, [f32; 3]);
impl_from_array!(Vec3f32, &[f32; 3]);

impl_array_from!(Vec3f32, [f32; 3]);
impl_array_from!(&Vec3f32, [f32; 3]);

impl_add!(Vec3f32, Vec3f32, Vec3f32);
impl_add!(Vec3f32, Vec3f32, &Vec3f32);
impl_add!(Vec3f32, &Vec3f32, Vec3f32);
impl_add!(Vec3f32, &Vec3f32, &Vec3f32);

impl_add_assign!(Vec3f32, Vec3f32);
impl_add_assign!(Vec3f32, &Vec3f32);

impl_sub!(Vec3f32, Vec3f32, Vec3f32);
impl_sub!(Vec3f32, Vec3f32, &Vec3f32);
impl_sub!(Vec3f32, &Vec3f32, Vec3f32);
impl_sub!(Vec3f32, &Vec3f32, &Vec3f32);

impl_sub_assign!(Vec3f32, Vec3f32);
impl_sub_assign!(Vec3f32, &Vec3f32);

impl_mul!(Vec3f32, Vec3f32, f32);
impl_mul!(Vec3f32, &Vec3f32, f32);
impl_mul_assign!(Vec3f32, f32);

impl_div!(Vec3f32, Vec3f32, f32);
impl_div!(Vec3f32, &Vec3f32, f32);
impl_div_assign!(Vec3f32, f32);

impl_neg!(Vec3f32, Vec3f32);
impl_neg!(Vec3f32, &Vec3f32);
