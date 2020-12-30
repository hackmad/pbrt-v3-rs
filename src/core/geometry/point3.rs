//! 3-D Points

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use num_traits::{Num, Zero};
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A 3-D point containing numeric values.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Point3<T> {
    /// X-coordinate.
    pub x: T,

    /// Y-coordinate.
    pub y: T,

    /// Z-coordinate.
    pub z: T,
}

/// 3-D point containing `Float` values.
pub type Point3f = Point3<Float>;

/// 3-D point containing `Int` values.
pub type Point3i = Point3<Int>;

impl<T: Num> Point3<T> {
    /// Creates a new 3-D point.
    ///
    /// * `x` - X-coordinate.
    /// * `y` - Y-coordinate.
    /// * `z` - Z-coordinate.
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Creates a new 3-D zero point.
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
        self.x.is_nan() || self.y.is_nan()
    }

    /// Returns a new point containing absolute values of the components.
    pub fn abs(&self) -> Self
    where
        T: Neg<Output = T> + PartialOrd + Copy,
    {
        Self::new(abs(self.x), abs(self.y), abs(self.z))
    }

    /// Returns a new point containing floor of values of the components.
    pub fn floor(&self) -> Self
    where
        T: num_traits::Float,
    {
        Self::new(self.x.floor(), self.y.floor(), self.z.floor())
    }

    /// Returns a new point containing ceil of values of the components.
    pub fn ceil(&self) -> Self
    where
        T: num_traits::Float,
    {
        Self::new(self.x.ceil(), self.y.ceil(), self.z.ceil())
    }

    /// Return the component-wise minimum coordinate values with another point.
    ///
    /// * `other` - The other point.
    pub fn min(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(
            min(self.x, other.x),
            min(self.y, other.y),
            min(self.z, other.z),
        )
    }

    /// Return the component-wise maximum coordinate values with another point.
    ///
    /// * `other` - The other point.
    pub fn max(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(
            max(self.x, other.x),
            max(self.y, other.y),
            max(self.z, other.z),
        )
    }

    /// Returns a new point with permuted coordinates according to given axes.
    ///
    /// * `x` - Axis to use for the x-coordinate of returned point.
    /// * `y` - Axis to use for the y-coordinate of returned point.
    /// * `z` - Axis to use for the z-coordinate of returned point.
    pub fn permute(&self, x: Axis, y: Axis, z: Axis) -> Self
    where
        T: Copy,
    {
        Self::new(self[x], self[y], self[z])
    }

    /// Returns the distance to another point.
    ///
    /// * `other` - The other point.
    pub fn distance(self, other: Self) -> T
    where
        T: num_traits::Float,
    {
        (self - other).length()
    }

    /// Returns the square of the istance to another point.
    ///
    /// * `other` - The other point.
    pub fn distance_squared(self, other: Self) -> T
    where
        T: num_traits::Float,
    {
        (self - other).length_squared()
    }
}

impl<T: Num> Add for Point3<T> {
    type Output = Self;

    /// Adds the given point and returns the result.
    ///
    /// * `other` - The point to add.
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Num + Copy> AddAssign for Point3<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The point to add.
    fn add_assign(&mut self, other: Self) {
        *self = Self::new(self.x + other.x, self.y + other.y, self.z + other.z);
    }
}

impl<T: Num> Add<Vector3<T>> for Point3<T> {
    type Output = Self;

    /// Offsets the point by the given vector.
    ///
    /// * `other` - The vector to add.
    fn add(self, other: Vector3<T>) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Num + Copy> AddAssign<Vector3<T>> for Point3<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The vector to add.
    fn add_assign(&mut self, other: Vector3<T>) {
        *self = Self::new(self.x + other.x, self.y + other.y, self.z + other.z);
    }
}

impl<T: Num> Sub for Point3<T> {
    type Output = Vector3<T>;

    /// Subtracts the given point and returns the vector towards that point.
    ///
    /// * `other` - The point to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Vector3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Num> Sub<Vector3<T>> for Point3<T> {
    type Output = Self;

    /// Subtracts the given vector and returns the result.
    ///
    /// * `other` - The point to subtract.
    fn sub(self, other: Vector3<T>) -> Self::Output {
        Self::Output::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Num + Copy> SubAssign<Vector3<T>> for Point3<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The vector to subtract.
    fn sub_assign(&mut self, other: Vector3<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y, self.z - other.z);
    }
}

impl<T: Num + Copy> Mul<T> for Point3<T> {
    type Output = Self;

    /// Scale the point.
    ///
    /// * `f` - The scaling factor.
    fn mul(self, f: T) -> Self::Output {
        Self::Output::new(f * self.x, f * self.y, f * self.z)
    }
}

macro_rules! premul {
    ($t: ty) => {
        impl Mul<Point3<$t>> for $t {
            type Output = Point3<$t>;
            /// Scale the point.
            ///
            /// * `p` - The point.
            fn mul(self, p: Point3<$t>) -> Point3<$t> {
                Point3::<$t>::new(self * p.x, self * p.y, self * p.z)
            }
        }

        impl Mul<&Point3<$t>> for $t {
            type Output = Point3<$t>;
            /// Scale the point.
            ///
            /// * `p` - The point.
            fn mul(self, p: &Point3<$t>) -> Point3<$t> {
                Point3::<$t>::new(self * p.x, self * p.y, self * p.z)
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

impl<T: Num + Copy> MulAssign<T> for Point3<T> {
    /// Scale and assign the result to the point.
    ///
    /// * `f` - The scaling factor.
    fn mul_assign(&mut self, f: T) {
        *self = Self::new(f * self.x, f * self.y, f * self.z);
    }
}

impl<T: Num + Copy> Div<T> for Point3<T> {
    type Output = Self;

    /// Scale the point by 1/f.
    ///
    /// * `f` - The scaling factor.
    fn div(self, f: T) -> Self::Output {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        Self::Output::new(inv * self.x, inv * self.y, inv * self.z)
    }
}

impl<T: Num + Copy> DivAssign<T> for Point3<T> {
    /// Scale the point by 1/f and assign the result to the point.
    ///
    /// * `f` - The scaling factor.
    fn div_assign(&mut self, f: T) {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        *self = Self::new(inv * self.x, inv * self.y, inv * self.z);
    }
}

impl<T: Num + Neg<Output = T>> Neg for Point3<T> {
    type Output = Self;

    /// Flip the point's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y, -self.z)
    }
}

impl<T> Index<Axis> for Point3<T> {
    type Output = T;

    /// Index the point by an axis to get the immutable coordinate axis value.
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

impl<T> Index<usize> for Point3<T> {
    type Output = T;

    /// Index the point by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index(&self, axis: usize) -> &Self::Output {
        &self[Axis::from(axis)]
    }
}

impl<T> IndexMut<Axis> for Point3<T> {
    /// Index the point by an axis to get a mutable coordinate axis value.
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

impl<T> IndexMut<usize> for Point3<T> {
    /// Index the point by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index_mut(&mut self, axis: usize) -> &mut Self::Output {
        &mut self[Axis::from(axis)]
    }
}

impl<T> From<Vector3<T>> for Point3<T> {
    /// Convert a 3-D vector to a 3-D point.
    ///
    /// * `v` - 3-D vector.
    fn from(v: Vector3<T>) -> Self {
        Self {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

impl From<Point3i> for Point3f {
    /// Convert a `Point3i` to `Point3f`.
    ///
    /// * `b` - The `Point3i` to convert.
    fn from(p: Point3i) -> Self {
        Self {
            x: p.x as Float,
            y: p.y as Float,
            z: p.z as Float,
        }
    }
}

impl From<Point3f> for Point3i {
    /// Convert a `Point3f` to `Point3i`.
    ///
    /// * `b` - The `Point3f` to convert.
    fn from(p: Point3f) -> Self {
        Self {
            x: p.x as Int,
            y: p.y as Int,
            z: p.z as Int,
        }
    }
}

impl<T: fmt::Display> fmt::Display for Point3<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.x, self.y, self.z)
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
    fn zero_point() {
        assert!(Point3::new(0, 0, 0) == Point3::zero());
        assert!(Point3::new(0.0, 0.0, 0.0) == Point3::zero());
    }

    #[test]
    fn has_nans() {
        assert!(!Point3::new(0.0, 0.0, 0.0).has_nans());
        assert!(Point3::new(f32::NAN, f32::NAN, f32::NAN).has_nans());
        assert!(Point3::new(f64::NAN, f64::NAN, f64::NAN).has_nans());
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_i64() {
        Point3::<i64>::zero() / 0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_f64() {
        Point3::<f64>::new(1.0, 1.0, 1.0) / 0.0;
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

    prop_point3!(point3_i32, i32, -100..100i32, -100..100i32, -100..100i32);
    prop_point3!(
        point3_f32,
        f32,
        -100.0..100.0f32,
        -100.0..100.0f32,
        -100.0..100.0f32
    );

    prop_vector3!(vector3_i32, i32, -100..100i32, -100..100i32, -100..100i32);
    prop_vector3!(
        vector3_f32,
        f32,
        -100.0..100.0f32,
        -100.0..100.0f32,
        -100.0..100.0f32
    );

    proptest! {
        #[test]
        fn distance_squared_f32(p1 in point3_f32(), p2 in point3_f32()) {
            let expected =
                (p1.x - p2.x) * (p1.x - p2.x) +
                (p1.y - p2.y) * (p1.y - p2.y) +
                (p1.z - p2.z) * (p1.z - p2.z);
            prop_assert_eq!(p1.distance_squared(p2), expected);
        }

        #[test]
        fn length_f32(p1 in point3_f32(), p2 in point3_f32()) {
            let expected =
                (p1.x - p2.x) * (p1.x - p2.x) +
                (p1.y - p2.y) * (p1.y - p2.y) +
                (p1.z - p2.z) * (p1.z - p2.z);
            let expected = expected.sqrt();
            prop_assert_eq!(p1.distance(p2), expected);
        }

        #[test]
        fn abs_i32(p in point3_i32()) {
            prop_assert_eq!(p.abs(), Point3::new(abs(p.x), abs(p.y), abs(p.z)));
        }

        #[test]
        fn floor_f32(p in point3_f32()) {
            prop_assert_eq!(p.floor(), Point3::new(p.x.floor(), p.y.floor(), p.z.floor()));
        }

        #[test]
        fn ceil_f32(p in point3_f32()) {
            prop_assert_eq!(p.ceil(), Point3::new(p.x.ceil(), p.y.ceil(), p.z.ceil()));
        }

        #[test]
        fn min_f32(p1 in point3_f32(), p2 in point3_f32()) {
            prop_assert_eq!(p1.min(&p2), Point3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z)));
        }

        #[test]
        fn max_f32(p1 in point3_f32(), p2 in point3_f32()) {
            prop_assert_eq!(p1.max(&p2), Point3::new(p1.x.max(p2.x), p1.y.max(p2.y), p1.z.max(p2.z)));
        }

        #[test]
        fn permute_f32(
            p in point3_f32(),
            a1 in axis_3d_strategy(),
            a2 in axis_3d_strategy(),
            a3 in axis_3d_strategy(),
        ) {
            let permuted_p = p.permute(a1, a2, a3);
            prop_assert_eq!(permuted_p.x, p[a1]);
            prop_assert_eq!(permuted_p.y, p[a2]);
            prop_assert_eq!(permuted_p.z, p[a3]);
        }

        #[test]
        fn lerp_edge_case_f32(p1 in point3_f32(), p2 in point3_f32()) {
            prop_assert_eq!(lerp(0.0, p1, p2), p1);
            prop_assert_eq!(lerp(1.0, p1, p2), p2);
        }

        #[test]
        fn lerp_f32(p1 in point3_f32(), p2 in point3_f32(), t in -2.0..2.0f32) {
            prop_assert_eq!(lerp(t, p1, p2), (1.0 - t) * p1 + t * p2);
        }

        #[test]
        fn add_point_i32(p1 in point3_i32(), p2 in point3_i32()) {
            prop_assert_eq!(p1 + p2, Point3::new(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z));
        }

        #[test]
        fn add_point_f32(p1 in point3_f32(), p2 in point3_f32()) {
            prop_assert_eq!(p1 + p2, Point3::new(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z));
        }

        #[test]
        fn add_vector_f32(p in point3_f32(), v in vector3_f32()) {
            prop_assert_eq!(p + v, Point3::new(p.x + v.x, p.y + v.y, p.z + v.z));
        }

        #[test]
        fn add_assign_vector_i32(p in point3_i32(), v in vector3_i32()) {
            let mut p1 = p;
            p1 += v;
            prop_assert_eq!(p1, Point3::new(p.x + v.x, p.y + v.y, p.z + v.z));
        }

        #[test]
        fn add_assign_vector_f32(p in point3_f32(), v in vector3_f32()) {
            let mut p1 = p;
            p1 += v;
            prop_assert_eq!(p1, Point3::new(p.x + v.x, p.y + v.y, p.z + v.z));
        }

        #[test]
        fn sub_point_i32(p1 in point3_i32(), p2 in point3_i32()) {
            prop_assert_eq!(p1 - p2, Vector3::new(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z));
        }

        #[test]
        fn sub_point_f32(p1 in point3_f32(), p2 in point3_f32()) {
            prop_assert_eq!(p1 - p2, Vector3::new(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z));
        }

        #[test]
        fn sub_vector_f32(p in point3_f32(), v in vector3_f32()) {
            prop_assert_eq!(p - v, Point3::new(p.x - v.x, p.y - v.y, p.z - v.z));
        }

        #[test]
        fn sub_assign_vector_i32(p in point3_i32(), v in vector3_i32()) {
            let mut p1 = p;
            p1 -= v;
            prop_assert_eq!(p1, Point3::new(p.x - v.x, p.y - v.y, p.z - v.z));
        }

        #[test]
        fn sub_assign_vector_f32(p in point3_f32(), v in vector3_f32()) {
            let mut p1 = p;
            p1 -= v;
            prop_assert_eq!(p1, Point3::new(p.x - v.x, p.y - v.y, p.z - v.z));
        }

        #[test]
        fn mul_i32(p in point3_i32(), f in range_i32()) {
            let expected = Point3::new(p.x * f, p.y * f, p.z * f);
            prop_assert_eq!(p * f, expected);
            prop_assert_eq!(f * p, expected);
        }

        #[test]
        fn mul_f32(p in point3_f32(), f in range_f32()) {
            let expected = Point3::new(p.x * f, p.y * f, p.z * f);
            prop_assert_eq!(p * f, expected);
            prop_assert_eq!(f * p, expected);
        }

        #[test]
        fn mul_assign_i32(p in point3_i32(), f in range_i32()) {
            let mut p1 = p;
            p1 *= f;
            prop_assert_eq!(p1, Point3::new(p.x * f, p.y * f, p.z * f));
        }

        #[test]
        fn mul_assign_f32(p in point3_f32(), f in range_f32()) {
            let mut p1 = p;
            p1 *= f;
            prop_assert_eq!(p1, Point3::new(p.x * f, p.y * f, p.z * f));
        }

        #[test]
        fn div_i32(
            p in point3_i32(),
            f in (-100..100i32).prop_filter("non-zero", |x| *x != 0)
        ) {
            let s = 1 / f;
            prop_assert_eq!(p / f, Point3::new(p.x * s, p.y * s, p.z * s));
        }

        #[test]
        fn div_f32(p in point3_f32(), f in non_zero_f32()) {
            let s = 1.0 / f;
            prop_assert_eq!(p / f, Point3::new(p.x * s, p.y * s, p.z * s));
        }

        #[test]
        fn div_assign_i32(p in point3_i32(), f in non_zero_i32()) {
            let mut p1 = p;
            p1 /= f;

            let s = 1 / f;
            prop_assert_eq!(p1, Point3::new(p.x * s, p.y * s, p.z * s));
        }

        #[test]
        fn div_assign_f32(p in point3_f32(), f in non_zero_f32()) {
            let mut p1 = p;
            p1 /= f;

            let s = 1.0 / f;
            prop_assert_eq!(p1, Point3::new(p.x * s, p.y * s, p.z * s));
        }

        #[test]
        fn neg_i32(p in point3_i32()) {
            prop_assert_eq!(-p, Point3::new(-p.x, -p.y, -p.z));
            prop_assert_eq!(--p, p);
        }

        #[test]
        fn neg_f32(p in point3_f32()) {
            prop_assert_eq!(-p, Point3::new(-p.x, -p.y, -p.z));
            prop_assert_eq!(--p, p);
        }

        #[test]
        fn index_i32(p in point3_i32()) {
            prop_assert_eq!(p[Axis::X], p.x);
            prop_assert_eq!(p[Axis::Y], p.y);
            prop_assert_eq!(p[Axis::Z], p.z);
        }

        #[test]
        fn index_f32(p in point3_f32()) {
            prop_assert_eq!(p[Axis::X], p.x);
            prop_assert_eq!(p[Axis::Y], p.y);
            prop_assert_eq!(p[Axis::Z], p.z);
        }

        #[test]
        fn index_mut_i32(p in point3_i32()) {
            let mut p1 = Point3::new(-200, 200, -200);
            p1[Axis::X] = p.x;
            p1[Axis::Y] = p.y;
            p1[Axis::Z] = p.z;
            prop_assert_eq!(p1, p);
        }

        #[test]
        fn index_mut_f32(p in point3_f32()) {
            let mut p1 = Point3::new(-200.0, 200.0, -200.0);
            p1[Axis::X] = p.x;
            p1[Axis::Y] = p.y;
            p1[Axis::Z] = p.z;
            prop_assert_eq!(p1, p);
        }
    }
}
