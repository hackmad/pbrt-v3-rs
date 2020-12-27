//! 2-D Vectors

#![allow(dead_code)]
use super::common::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use num_traits::{Num, Zero};
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A 2-D vector containing numeric values.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vector2<T> {
    /// X-coordinate.
    pub x: T,

    /// Y-coordinate.
    pub y: T,
}

/// 2-D vector containing `Float` values.
pub type Vector2f = Vector2<Float>;

/// 2-D vector containing `Int` values.
pub type Vector2i = Vector2<Int>;

impl<T: Num> Vector2<T> {
    /// Creates a new 2-D vector.
    ///
    /// * `x` - X-coordinate.
    /// * `y` - Y-coordinate.
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    /// Creates a new 2-D zero vector.
    pub fn zero() -> Self
    where
        T: Zero,
    {
        Self::new(T::zero(), T::zero())
    }

    /// Returns true if either coordinate is NaN.
    pub fn has_nans(&self) -> bool
    where
        T: num_traits::Float,
    {
        self.x.is_nan() || self.y.is_nan()
    }

    /// Returns the square of the vector's length.
    pub fn length_squared(&self) -> T
    where
        T: Mul<Output = T> + Add<Output = T> + Copy,
    {
        self.x * self.x + self.y * self.y
    }

    /// Returns the vector's length.
    pub fn length(&self) -> T
    where
        T: num_traits::Float,
    {
        self.length_squared().sqrt()
    }

    /// Returns the unit vector.
    pub fn normalize(&self) -> Self
    where
        T: num_traits::Float,
    {
        *self / self.length()
    }

    /// Returns a new vector containing absolute values of the components.
    pub fn abs(&self) -> Self
    where
        T: Neg<Output = T> + PartialOrd + Copy,
    {
        Self::new(abs(self.x), abs(self.y))
    }

    /// Returns the smallest coordinate value.
    pub fn min_component(&self) -> T
    where
        T: PartialOrd + Copy,
    {
        if self.x.lt(&self.y) {
            self.x
        } else {
            self.y
        }
    }

    /// Returns the largest coordinate value.
    pub fn max_component(&self) -> T
    where
        T: PartialOrd + Copy,
    {
        if self.x.gt(&self.y) {
            self.x
        } else {
            self.y
        }
    }

    /// Returns the axis with largest coordinate value.
    pub fn max_dimension(&self) -> Axis
    where
        T: PartialOrd + Copy,
    {
        if self.x.gt(&self.y) {
            Axis::X
        } else {
            Axis::Y
        }
    }

    /// Return the component-wise minimum coordinate values with another vector.
    ///
    /// * `other` -  The other vector.
    pub fn min(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(min(self.x, other.x), min(self.y, other.y))
    }

    /// Return the component-wise maximum coordinate values with another vector.
    ///
    /// * `other` -  The other vector.
    pub fn max(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(max(self.x, other.x), max(self.y, other.y))
    }

    /// Returns a new vector with permuted coordinates according to given axes.
    ///
    /// * `x` -  Axis to use for the x-coordinate of returned vector.
    /// * `y` -  Axis to use for the y-coordinate of returned vector.
    pub fn permute(&self, x: Axis, y: Axis) -> Self
    where
        T: Copy,
    {
        Self::new(self[x], self[y])
    }
}

impl<T: Num + Neg<Output = T> + PartialOrd + Copy> Dot<Vector2<T>> for Vector2<T> {
    type Output = T;

    /// Returns the dot product with another vector.
    ///
    /// * `other` -  The other vector.
    fn dot(&self, other: &Vector2<T>) -> T {
        self.x * other.x + self.y * other.y
    }
}

impl<T: Num> Add for Vector2<T> {
    type Output = Self;

    /// Adds the given vector and returns the result.
    ///
    /// * `other` -  The vector to add.
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Num + Copy> AddAssign for Vector2<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` -  The vector to add.
    fn add_assign(&mut self, other: Self) {
        *self = Self::new(self.x + other.x, self.y + other.y);
    }
}

impl<T: Num> Sub for Vector2<T> {
    type Output = Self;

    /// Subtracts the given vector and returns the result.
    ///
    /// * `other` -  The vector to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Self::Output::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Num + Copy> SubAssign for Vector2<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` -  The vector to subtract.
    fn sub_assign(&mut self, other: Self) {
        *self = Self::new(self.x - other.x, self.y - other.y);
    }
}

impl<T: Num + Copy> Mul<T> for Vector2<T> {
    type Output = Self;

    /// Scale the vector.
    ///
    /// * `f` -  The scaling factor.
    fn mul(self, f: T) -> Self::Output {
        Self::Output::new(f * self.x, f * self.y)
    }
}

macro_rules! premul {
    ($t: ty) => {
        impl Mul<Vector2<$t>> for $t {
            type Output = Vector2<$t>;
            /// Scale the vector.
            ///
            /// * `v` -  The vector.
            fn mul(self, v: Vector2<$t>) -> Vector2<$t> {
                Vector2::<$t>::new(self * v.x, self * v.y)
            }
        }

        impl Mul<&Vector2<$t>> for $t {
            type Output = Vector2<$t>;
            /// Scale the vector.
            ///
            /// * `v` -  The vector.
            fn mul(self, v: &Vector2<$t>) -> Vector2<$t> {
                Vector2::<$t>::new(self * v.x, self * v.y)
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

impl<T: Num + Copy> MulAssign<T> for Vector2<T> {
    /// Scale and assign the result to the vector.
    ///
    /// * `f` -  The scaling factor.
    fn mul_assign(&mut self, f: T) {
        *self = Self::new(f * self.x, f * self.y)
    }
}

impl<T: Num + Copy> Div<T> for Vector2<T> {
    type Output = Self;

    /// Scale the vector by 1/f.
    ///
    /// * `f` -  The scaling factor.
    fn div(self, f: T) -> Self::Output {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        Self::Output::new(inv * self.x, inv * self.y)
    }
}

impl<T: Num + Copy> DivAssign<T> for Vector2<T> {
    /// Scale the vector by 1/f and assign the result to the vector.
    ///
    /// * `f` -  The scaling factor.
    fn div_assign(&mut self, f: T) {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        *self = Self::new(inv * self.x, inv * self.y);
    }
}

impl<T: Num + Neg<Output = T>> Neg for Vector2<T> {
    type Output = Self;

    /// Flip the vector's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y)
    }
}

impl<T> Index<Axis> for Vector2<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 2-D coordinate axis.
    fn index(&self, axis: Axis) -> &Self::Output {
        match axis {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            _ => panic!("Invalid axis for std::Index on Vector2<T>"),
        }
    }
}

impl<T> Index<usize> for Vector2<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 2-D coordinate axis.
    fn index(&self, axis: usize) -> &Self::Output {
        &self[Axis::from(axis)]
    }
}

impl<T> IndexMut<Axis> for Vector2<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` - A 2-D coordinate axis.
    fn index_mut(&mut self, axis: Axis) -> &mut Self::Output {
        match axis {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            _ => panic!("Invalid axis for std::IndexMut on Vector2<T>"),
        }
    }
}

impl<T> IndexMut<usize> for Vector2<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 2-D coordinate axis.
    fn index_mut(&mut self, axis: usize) -> &mut Self::Output {
        &mut self[Axis::from(axis)]
    }
}

impl<T> From<Point2<T>> for Vector2<T> {
    /// Convert a 2-D point to a 2-D vector.
    ///
    /// * `p` - 2-D point.
    fn from(p: Point2<T>) -> Self {
        Self { x: p.x, y: p.y }
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
        assert!(Vector2::new(0, 0) == Vector2::zero());
        assert!(Vector2::new(0.0, 0.0) == Vector2::zero());
    }

    #[test]
    fn has_nans() {
        assert!(!Vector2::new(0.0, 0.0).has_nans());
        assert!(Vector2::new(f32::NAN, f32::NAN).has_nans());
        assert!(Vector2::new(f64::NAN, f64::NAN).has_nans());
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f64() {
        Vector2::<f64>::zero().normalize();
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f32() {
        Vector2::<f32>::zero().normalize();
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_i64() {
        Vector2::<i64>::zero() / 0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_f64() {
        Vector2::<f64>::new(1.0, 1.0) / 0.0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn invalid_index() {
        let z = Vector2::<i64>::zero()[Axis::Z];
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn invalid_index_mut() {
        let mut v = Vector2::<i64>::zero();
        v[Axis::Z] = 1;
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

    prop_vector2!(vector2_i32, i32, -100..100i32, -100..100i32);
    prop_vector2!(vector2_f32, f32, -100.0..100.0f32, -100.0..100.0f32);

    proptest! {
        #[test]
        fn length_squared_i32(v in vector2_i32()) {
            prop_assert_eq!(v.length_squared(), v.x * v.x + v.y * v.y);
        }

        #[test]
        fn length_squared_f32(v in vector2_f32()) {
            prop_assert_eq!(v.length_squared(), v.x * v.x + v.y * v.y);
        }

        #[test]
        fn length_f32(v in vector2_f32()) {
            prop_assert_eq!(v.length(), (v.x * v.x + v.y * v.y).sqrt());
        }

        #[test]
        fn normalize_f32(v in vector2_f32()) {
            // Since we do 1.0 / l in implementation we have to do the
            // same here. Doing Vector2::new(x / l, y / l) will not work
            // for some of the floating point values due to precision
            // errors.
            let f = 1.0 / (v.x * v.x + v.y * v.y).sqrt();
            prop_assert_eq!(v.normalize(), Vector2::new(v.x * f, v.y * f));
        }

        #[test]
        fn abs_i32(v in vector2_i32()) {
            prop_assert_eq!(v.abs(), Vector2::new(abs(v.x), abs(v.y)));
        }

        #[test]
        fn abs_f32(v in vector2_f32()) {
            prop_assert_eq!(v.abs(), Vector2::new(abs(v.x), abs(v.y)));
        }

        #[test]
        fn dot_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1.dot(&v2), v1.x * v2.x + v1.y * v2.y);
        }

        #[test]
        fn abs_dot_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1.abs_dot(&v2), (v1.x * v2.x + v1.y * v2.y).abs());
        }

        #[test]
        fn min_component_f32(v in vector2_f32()) {
            prop_assert_eq!(v.min_component(), v.x.min(v.y));
        }

        #[test]
        fn max_component_f32(v in vector2_f32()) {
            prop_assert_eq!(v.max_component(), v.x.max(v.y));
        }

        #[test]
        fn max_dimension_f32(v in vector2_f32()) {
            let dim = if v.x > v.y { Axis::X } else { Axis::Y };
            prop_assert_eq!(v.max_dimension(), dim);
        }

        #[test]
        fn min_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1.min(&v2), Vector2::new(v1.x.min(v2.x), v1.y.min(v2.y)));
        }

        #[test]
        fn max_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1.max(&v2), Vector2::new(v1.x.max(v2.x), v1.y.max(v2.y)));
        }

        #[test]
        fn permute_f32(v in vector2_f32(), a1 in axis_2d_strategy(), a2 in axis_2d_strategy()) {
            let permuted_v = v.permute(a1, a2);
            prop_assert_eq!(permuted_v.x, v[a1]);
            prop_assert_eq!(permuted_v.y, v[a2]);
        }

        #[test]
        fn add_i32(v1 in vector2_i32(), v2 in vector2_i32()) {
            prop_assert_eq!(v1 + v2, Vector2::new(v1.x + v2.x, v1.y + v2.y));
        }

        #[test]
        fn add_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1 + v2, Vector2::new(v1.x + v2.x, v1.y + v2.y));
        }

        #[test]
        fn add_assign_i32(v1 in vector2_i32(), v2 in vector2_i32()) {
            let mut v = v1;
            v += v2;
            prop_assert_eq!(v, Vector2::new(v1.x + v2.x, v1.y + v2.y));
        }

        #[test]
        fn add_assign_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            let mut v = v1;
            v += v2;
            prop_assert_eq!(v, Vector2::new(v1.x + v2.x, v1.y + v2.y));
        }

        #[test]
        fn sub_i32(v1 in vector2_i32(), v2 in vector2_i32()) {
            prop_assert_eq!(v1 - v2, Vector2::new(v1.x - v2.x, v1.y - v2.y));
        }

        #[test]
        fn sub_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            prop_assert_eq!(v1 - v2, Vector2::new(v1.x - v2.x, v1.y - v2.y));
        }

        #[test]
        fn sub_assign_i32(v1 in vector2_i32(), v2 in vector2_i32()) {
            let mut v = v1;
            v -= v2;
            prop_assert_eq!(v, Vector2::new(v1.x - v2.x, v1.y - v2.y));
        }

        #[test]
        fn sub_assign_f32(v1 in vector2_f32(), v2 in vector2_f32()) {
            let mut v = v1;
            v -= v2;
            prop_assert_eq!(v, Vector2::new(v1.x - v2.x, v1.y - v2.y));
        }

        #[test]
        fn mul_i32(v in vector2_i32(), f in range_i32()) {
            let expected = Vector2::new(v.x * f, v.y * f);
            prop_assert_eq!(v * f, expected);
            prop_assert_eq!(f * v, expected);
        }

        #[test]
        fn mul_f32(v in vector2_f32(), f in range_f32()) {
            let expected = Vector2::new(v.x * f, v.y * f);
            prop_assert_eq!(v * f, expected);
            prop_assert_eq!(f * v, expected);
        }

        #[test]
        fn mul_assign_i32(v in vector2_i32(), f in range_i32()) {
            let mut v1 = v;
            v1 *= f;
            prop_assert_eq!(v1, Vector2::new(v.x * f, v.y * f));
        }

        #[test]
        fn mul_assign_f32(v in vector2_f32(), f in range_f32()) {
            let mut v1 = v;
            v1 *= f;
            prop_assert_eq!(v1, Vector2::new(v.x * f, v.y * f));
        }

        #[test]
        fn div_i32(
            v in vector2_i32(),
            f in (-100..100i32).prop_filter("non-zero", |x| *x != 0)
        ) {
            let s = 1 / f;
            prop_assert_eq!(v / f, Vector2::new(v.x * s, v.y * s));
        }

        #[test]
        fn div_f32(v in vector2_f32(), f in non_zero_f32()) {
            let s = 1.0 / f;
            prop_assert_eq!(v / f, Vector2::new(v.x * s, v.y * s));
        }

        #[test]
        fn div_assign_i32(v in vector2_i32(), f in non_zero_i32()) {
            let mut v1 = v;
            v1 /= f;

            let s = 1 / f;
            prop_assert_eq!(v1, Vector2::new(v.x * s, v.y * s));
        }

        #[test]
        fn div_assign_f32(v in vector2_f32(), f in non_zero_f32()) {
            let mut v1 = v;
            v1 /= f;

            let s = 1.0 / f;
            prop_assert_eq!(v1, Vector2::new(v.x * s, v.y * s));
        }

        #[test]
        fn neg_i32(v in vector2_i32()) {
            prop_assert_eq!(-v, Vector2::new(-v.x, -v.y));
            prop_assert_eq!(--v, v);
        }

        #[test]
        fn neg_f32(v in vector2_f32()) {
            prop_assert_eq!(-v, Vector2::new(-v.x, -v.y));
            prop_assert_eq!(--v, v);
        }

        #[test]
        fn index_i32(v in vector2_i32()) {
            prop_assert_eq!(v[Axis::X], v.x);
            prop_assert_eq!(v[Axis::Y], v.y);
        }

        #[test]
        fn index_f32(v in vector2_f32()) {
            prop_assert_eq!(v[Axis::X], v.x);
            prop_assert_eq!(v[Axis::Y], v.y);
        }

        #[test]
        fn index_mut_i32(v in vector2_i32()) {
            let mut v1 = Vector2::new(-200, 200);
            v1[Axis::X] = v.x;
            v1[Axis::Y] = v.y;
            prop_assert_eq!(v1, v);
        }

        #[test]
        fn index_mut_f32(v in vector2_f32()) {
            let mut v1 = Vector2::new(-200.0, 200.0);
            v1[Axis::X] = v.x;
            v1[Axis::Y] = v.y;
            prop_assert_eq!(v1, v);
        }
    }
}
