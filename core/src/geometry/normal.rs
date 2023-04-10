//! 3-D normals

use super::common::*;
use crate::geometry::*;
use crate::pbrt::*;
use num_traits::{Num, Zero};
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

/// A 3-D normal containing numeric values.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Normal3<T> {
    /// X-coordinate.
    pub x: T,

    /// Y-coordinate.
    pub y: T,

    /// Y-coordinate.
    pub z: T,
}

/// 3-D normal containing `Float` values.
pub type Normal3f = Normal3<Float>;
impl Normal3f {
    /// Zero normal.
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };
}

impl<T: Num> Normal3<T> {
    /// Creates a new 3-D normal.
    ///
    /// * `x` - X-coordinate.
    /// * `y` - Y-coordinate.
    /// * `z` - Z-coordinate.
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Creates a new 3-D zero normal.
    pub fn zero() -> Self
    where
        T: Zero,
    {
        Self::new(T::zero(), T::zero(), T::zero())
    }

    /// Returns true if either coordinate is NaN.
    pub fn has_nans(&self) -> bool
    where
        T: num_traits::Float,
    {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }

    /// Returns the square of the normal's length.
    pub fn length_squared(&self) -> T
    where
        T: Mul<Output = T> + Add<Output = T> + Copy,
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the normal's length.
    pub fn length(&self) -> T
    where
        T: num_traits::Float,
    {
        self.length_squared().sqrt()
    }

    /// Returns the unit normal.
    pub fn normalize(&self) -> Self
    where
        T: num_traits::Float,
    {
        *self / self.length()
    }

    /// Returns a new normal containing absolute values of the components.
    pub fn abs(&self) -> Self
    where
        T: Neg<Output = T> + PartialOrd + Copy,
    {
        Self::new(abs(self.x), abs(self.y), abs(self.z))
    }
}

impl<T: Num + Neg<Output = T> + PartialOrd + Copy> Dot<Normal3<T>> for Normal3<T> {
    type Output = T;

    /// Returns the dot product with another normal.
    ///
    /// * `other` - The other normal.
    fn dot(&self, other: &Self) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl<T: Num + Neg<Output = T> + PartialOrd + Copy> Dot<Vector3<T>> for Normal3<T> {
    type Output = T;

    /// Returns the dot product with another vector.
    ///
    /// * `other` - The other vector.
    fn dot(&self, other: &Vector3<T>) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl<T: Num + Copy> Cross<Normal3<T>> for Normal3<T> {
    type Output = Self;

    /// Returns the cross product with another normal.
    ///
    /// * `other` - The other vector.
    fn cross(&self, other: &Normal3<T>) -> Self::Output {
        Self::new(
            (self.y * other.z) - (self.z * other.y),
            (self.z * other.x) - (self.x * other.z),
            (self.x * other.y) - (self.y * other.x),
        )
    }
}

impl<T: Num + Copy> Cross<Vector3<T>> for Normal3<T> {
    type Output = Self;

    /// Returns the cross product with another vector.
    ///
    /// * `other` - The other vector.
    fn cross(&self, other: &Vector3<T>) -> Self::Output {
        Self::new(
            (self.y * other.z) - (self.z * other.y),
            (self.z * other.x) - (self.x * other.z),
            (self.x * other.y) - (self.y * other.x),
        )
    }
}

// Implement FaceForward trait which allows pointing vectors in the same hemisphere as another normal/vector.
impl<T: Num + Neg<Output = T> + PartialOrd + Copy> FaceForward<T, Vector3<T>> for Normal3<T> {}
impl<T: Num + Neg<Output = T> + PartialOrd + Copy> FaceForward<T, Normal3<T>> for Normal3<T> {}

impl<T: Num> Add for Normal3<T> {
    type Output = Self;

    /// Adds the given normal and returns the result.
    ///
    /// * `other` - The normal to add.
    fn add(self, other: Self) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Num + Copy> Add<&Normal3<T>> for Normal3<T> {
    type Output = Self;

    /// Adds the given normal and returns the result.
    ///
    /// * `other` - The normal to add.
    fn add(self, other: &Normal3<T>) -> Self::Output {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Num + Copy> AddAssign for Normal3<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The normal to add.
    fn add_assign(&mut self, other: Self) {
        *self = Self::new(self.x + other.x, self.y + other.y, self.z + other.z);
    }
}

impl<T: Num + Copy> AddAssign<&Normal3<T>> for Normal3<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The normal to add.
    fn add_assign(&mut self, other: &Normal3<T>) {
        *self = Self::new(self.x + other.x, self.y + other.y, self.z + other.z);
    }
}

impl<T: Num> Sub for Normal3<T> {
    type Output = Normal3<T>;

    /// Subtracts the given normal and returns the result.
    ///
    /// * `other` - The normal to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Num + Copy> Sub<&Normal3<T>> for Normal3<T> {
    type Output = Normal3<T>;

    /// Subtracts the given normal and returns the result.
    ///
    /// * `other` - The normal to subtract.
    fn sub(self, other: &Normal3<T>) -> Self::Output {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Num + Copy> SubAssign for Normal3<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The normal to subtract.
    fn sub_assign(&mut self, other: Self) {
        *self = Self::new(self.x - other.x, self.y - other.y, self.z - other.z);
    }
}

impl<T: Num + Copy> SubAssign<&Normal3<T>> for Normal3<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The normal to subtract.
    fn sub_assign(&mut self, other: &Normal3<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y, self.z - other.z);
    }
}

impl<T: Num + Copy> Mul<T> for Normal3<T> {
    type Output = Self;

    /// Scale the vector.
    ///
    /// * `f` - The scaling factor.
    fn mul(self, f: T) -> Self::Output {
        Self::new(f * self.x, f * self.y, f * self.z)
    }
}

macro_rules! premul {
    ($t: ty) => {
        impl Mul<Normal3<$t>> for $t {
            type Output = Normal3<$t>;
            /// Scale the normal.
            ///
            /// * `n` - The normal.
            fn mul(self, n: Normal3<$t>) -> Normal3<$t> {
                Normal3::new(self * n.x, self * n.y, self * n.z)
            }
        }

        impl Mul<&Normal3<$t>> for $t {
            type Output = Normal3<$t>;
            /// Scale the normal.
            ///
            /// * `n` - The normal.
            fn mul(self, n: &Normal3<$t>) -> Normal3<$t> {
                Normal3::new(self * n.x, self * n.y, self * n.z)
            }
        }
    };
}

premul!(f32);
premul!(f64);
premul!(i8);
premul!(i16);
premul!(i32);
premul!(i64);
premul!(u8);
premul!(u16);
premul!(u32);
premul!(u64);

impl<T: Num + Copy> MulAssign<T> for Normal3<T> {
    /// Scale and assign the result to the vector.
    ///
    /// * `f` - The scaling factor.
    fn mul_assign(&mut self, f: T) {
        *self = Self::new(f * self.x, f * self.y, f * self.z);
    }
}

impl<T: Num + Copy> Div<T> for Normal3<T> {
    type Output = Self;

    /// Scale the vector by 1/f.
    ///
    /// * `f` - The scaling factor.
    fn div(self, f: T) -> Self::Output {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        Self::Output::new(inv * self.x, inv * self.y, inv * self.z)
    }
}

impl<T: Num + Copy> DivAssign<T> for Normal3<T> {
    /// Scale the vector by 1/f and assign the result to the vector.
    ///
    /// * `f` - The scaling factor.
    fn div_assign(&mut self, f: T) {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        *self = Self::new(inv * self.x, inv * self.y, inv * self.z);
    }
}

impl<T: Num + Neg<Output = T>> Neg for Normal3<T> {
    type Output = Self;

    /// Flip the vector's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y, -self.z)
    }
}

impl<T: Num + Neg<Output = T> + Copy> Neg for &Normal3<T> {
    type Output = Normal3<T>;

    /// Flip the vector's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y, -self.z)
    }
}

impl<T> Index<Axis> for Normal3<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` - A 3-D coordinate axis.
    fn index(&self, axis: Axis) -> &Self::Output {
        match axis {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            Axis::Z => &self.z,
        }
    }
}

impl<T> Index<usize> for Normal3<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index(&self, axis: usize) -> &Self::Output {
        &self[Axis::from(axis)]
    }
}

impl<T> IndexMut<Axis> for Normal3<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` - A 3-D coordinate axis.
    fn index_mut(&mut self, axis: Axis) -> &mut Self::Output {
        match axis {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            Axis::Z => &mut self.z,
        }
    }
}

impl<T> IndexMut<usize> for Normal3<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index_mut(&mut self, axis: usize) -> &mut Self::Output {
        &mut self[Axis::from(axis)]
    }
}

impl<T> From<Vector3<T>> for Normal3<T> {
    /// Convert a 3-D vector to a 3-D normal.
    ///
    /// * `v` - 3-D vector.
    fn from(v: Vector3<T>) -> Self {
        Self { x: v.x, y: v.y, z: v.z }
    }
}

impl<T: Copy> From<&Vector3<T>> for Normal3<T> {
    /// Convert a 3-D vector to a 3-D normal.
    ///
    /// * `v` - 3-D vector.
    fn from(v: &Vector3<T>) -> Self {
        Self { x: v.x, y: v.y, z: v.z }
    }
}

impl<T: fmt::Display> fmt::Display for Normal3<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "<{}, {}, {}>", self.x, self.y, self.z)
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn zero_vector() {
        assert!(Normal3::new(0, 0, 0) == Normal3::zero());
        assert!(Normal3::new(0.0, 0.0, 0.0) == Normal3::zero());
    }

    #[test]
    fn has_nans() {
        assert!(!Normal3::new(0.0, 0.0, 0.0).has_nans());
        assert!(Normal3::new(f32::NAN, f32::NAN, f32::NAN).has_nans());
        assert!(Normal3::new(f64::NAN, f64::NAN, f64::NAN).has_nans());
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f64() {
        Normal3::<f64>::zero().normalize();
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f32() {
        Normal3::<f32>::zero().normalize();
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_i64() {
        Normal3::<i64>::zero() / 0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_f64() {
        Normal3::<f64>::new(1.0, 1.0, 1.0) / 0.0;
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

    prop_normal3!(normal3_i32, i32, -100..100i32, -100..100i32, -100..100i32);
    prop_normal3!(normal3_f32, f32, -100.0..100.0f32, -100.0..100.0f32, -100.0..100.0f32);

    proptest! {
        #[test]
        fn length_squared_i32(n in normal3_i32()) {
            prop_assert_eq!(n.length_squared(), n.x * n.x + n.y * n.y + n.z * n.z);
        }

        #[test]
        fn length_squared_f32(n in normal3_f32()) {
            prop_assert_eq!(n.length_squared(), n.x * n.x + n.y * n.y + n.z * n.z);
        }

        #[test]
        fn length_f32(n in normal3_f32()) {
            prop_assert_eq!(n.length(), (n.x * n.x + n.y * n.y + n.z * n.z).sqrt());
        }

        #[test]
        fn normalize_f32(n in normal3_f32()) {
            // Since we do 1.0 / l in implementation we have to do the
            // same here. Doing Normal3::new(x / l, y / l) will not work
            // for some of the floating point values due to precision
            // errors.
            let f = 1.0 / (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
            prop_assert_eq!(n.normalize(), Normal3::new(n.x * f, n.y * f, n.z * f));
        }

        #[test]
        fn dot_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            prop_assert_eq!(n1.dot(&n2), n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
        }

        #[test]
        fn abs_dot_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            prop_assert_eq!(n1.abs_dot(&n2), (n1.x * n2.x + n1.y * n2.y + n1.z * n2.z).abs());
        }

        #[test]
        fn add_i32(n1 in normal3_i32(), n2 in normal3_i32()) {
            prop_assert_eq!(n1 + n2, Normal3::new(n1.x + n2.x, n1.y + n2.y, n1.z + n2.z));
        }

        #[test]
        fn add_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            prop_assert_eq!(n1 + n2, Normal3::new(n1.x + n2.x, n1.y + n2.y, n1.z + n2.z));
        }

        #[test]
        fn add_assign_i32(n1 in normal3_i32(), n2 in normal3_i32()) {
            let mut n = n1;
            n += n2;
            prop_assert_eq!(n, Normal3::new(n1.x + n2.x, n1.y + n2.y, n1.z + n2.z));
        }

        #[test]
        fn add_assign_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            let mut n = n1;
            n += n2;
            prop_assert_eq!(n, Normal3::new(n1.x + n2.x, n1.y + n2.y, n1.z + n2.z));
        }

        #[test]
        fn sub_i32(n1 in normal3_i32(), n2 in normal3_i32()) {
            prop_assert_eq!(n1 - n2, Normal3::new(n1.x - n2.x, n1.y - n2.y, n1.z - n2.z));
        }

        #[test]
        fn sub_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            prop_assert_eq!(n1 - n2, Normal3::new(n1.x - n2.x, n1.y - n2.y, n1.z - n2.z));
        }

        #[test]
        fn sub_assign_i32(n1 in normal3_i32(), n2 in normal3_i32()) {
            let mut n = n1;
            n -= n2;
            prop_assert_eq!(n, Normal3::new(n1.x - n2.x, n1.y - n2.y, n1.z - n2.z));
        }

        #[test]
        fn sub_assign_f32(n1 in normal3_f32(), n2 in normal3_f32()) {
            let mut n = n1;
            n -= n2;
            prop_assert_eq!(n, Normal3::new(n1.x - n2.x, n1.y - n2.y, n1.z - n2.z));
        }

        #[test]
        fn mul_i32(n in normal3_i32(), f in range_i32()) {
            let expected = Normal3::new(n.x * f, n.y * f, n.z * f);
            prop_assert_eq!(n * f, expected);
            prop_assert_eq!(f * n, expected);
        }

        #[test]
        fn mul_f32(n in normal3_f32(), f in range_f32()) {
            let expected = Normal3::new(n.x * f, n.y * f, n.z * f);
            prop_assert_eq!(n * f, expected);
            prop_assert_eq!(f * n, expected);
        }

        #[test]
        fn mul_assign_i32(n in normal3_i32(), f in range_i32()) {
            let mut n1 = n;
            n1 *= f;
            prop_assert_eq!(n1, Normal3::new(n.x * f, n.y * f, n.z * f));
        }

        #[test]
        fn mul_assign_f32(n in normal3_f32(), f in range_f32()) {
            let mut n1 = n;
            n1 *= f;
            prop_assert_eq!(n1, Normal3::new(n.x * f, n.y * f, n.z * f));
        }

        #[test]
        fn div_i32(
            n in normal3_i32(),
            f in (-100..100i32).prop_filter("non-zero", |x| *x != 0)
        ) {
            let s = 1 / f;
            prop_assert_eq!(n / f, Normal3::new(n.x * s, n.y * s, n.z * s));
        }

        #[test]
        fn div_f32(n in normal3_f32(), f in non_zero_f32()) {
            let s = 1.0 / f;
            prop_assert_eq!(n / f, Normal3::new(n.x * s, n.y * s, n.z * s));
        }

        #[test]
        fn div_assign_i32(n in normal3_i32(), f in non_zero_i32()) {
            let mut n1 = n;
            n1 /= f;

            let s = 1 / f;
            prop_assert_eq!(n1, Normal3::new(n.x * s, n.y * s, n.z * s));
        }

        #[test]
        fn div_assign_f32(n in normal3_f32(), f in non_zero_f32()) {
            let mut n1 = n;
            n1 /= f;

            let s = 1.0 / f;
            prop_assert_eq!(n1, Normal3::new(n.x * s, n.y * s, n.z * s));
        }

        #[test]
        fn neg_i32(n in normal3_i32()) {
            prop_assert_eq!(-n, Normal3::new(-n.x, -n.y, -n.z));
            prop_assert_eq!(--n, n);
        }

        #[test]
        fn neg_f32(n in normal3_f32()) {
            prop_assert_eq!(-n, Normal3::new(-n.x, -n.y, -n.z));
            prop_assert_eq!(--n, n);
        }

        #[test]
        fn index_i32(n in normal3_i32()) {
            prop_assert_eq!(n[Axis::X], n.x);
            prop_assert_eq!(n[Axis::Y], n.y);
        }

        #[test]
        fn index_f32(n in normal3_f32()) {
            prop_assert_eq!(n[Axis::X], n.x);
            prop_assert_eq!(n[Axis::Y], n.y);
        }

        #[test]
        fn index_mut_i32(n in normal3_i32()) {
            let mut n1 = Normal3::new(-200, 200, -200);
            n1[Axis::X] = n.x;
            n1[Axis::Y] = n.y;
            n1[Axis::Z] = n.z;
            prop_assert_eq!(n1, n);
        }

        #[test]
        fn index_mut_f32(n in normal3_f32()) {
            let mut n1 = Normal3::new(-200.0, 200.0, -200.0);
            n1[Axis::X] = n.x;
            n1[Axis::Y] = n.y;
            n1[Axis::Z] = n.z;
            prop_assert_eq!(n1, n);
        }
    }
}
