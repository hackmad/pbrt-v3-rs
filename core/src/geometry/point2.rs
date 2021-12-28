//! 2-D Points

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use num_traits::{Num, Zero};
use std::fmt;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A 2-D point containing numeric values.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Point2<T> {
    /// X-coordinate.
    pub x: T,

    /// Y-coordinate.
    pub y: T,
}

/// 2-D point containing `Float` values.
pub type Point2f = Point2<Float>;

/// 2-D point containing `Int` values.
pub type Point2i = Point2<Int>;

impl<T: Num> Point2<T> {
    /// Creates a new 2-D point.
    ///
    /// * `x` - X-coordinate.
    /// * `y` - Y-coordinate.
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }

    /// Creates a new 2-D zero point.
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

    /// Returns a new point containing absolute values of the components.
    pub fn abs(&self) -> Self
    where
        T: Neg<Output = T> + PartialOrd + Copy,
    {
        Self::new(abs(self.x), abs(self.y))
    }

    /// Returns a new point containing floor of values of the components.
    pub fn floor(&self) -> Self
    where
        T: num_traits::Float,
    {
        Self::new(self.x.floor(), self.y.floor())
    }

    /// Returns a new point containing ceil of values of the components.
    pub fn ceil(&self) -> Self
    where
        T: num_traits::Float,
    {
        Self::new(self.x.ceil(), self.y.ceil())
    }

    /// Return the component-wise minimum coordinate values with another point.
    ///
    /// * `other` - The other point.
    pub fn min(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(min(self.x, other.x), min(self.y, other.y))
    }

    /// Return the component-wise maximum coordinate values with another point.
    ///
    /// * `other` - The other point.
    pub fn max(&self, other: &Self) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self::new(max(self.x, other.x), max(self.y, other.y))
    }

    /// Returns a new point with permuted coordinates according to given axes.
    ///
    /// * `x` - Axis to use for the x-coordinate of returned point.
    /// * `y` - Axis to use for the y-coordinate of returned point.
    pub fn permute(&self, x: Axis, y: Axis) -> Self
    where
        T: Copy,
    {
        Self::new(self[x], self[y])
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

impl<T: Num> Add for Point2<T> {
    type Output = Self;

    /// Adds the given point and returns the result.
    ///
    /// * `other` - The point to add.
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Num + Copy> Add<&Point2<T>> for Point2<T> {
    type Output = Point2<T>;

    /// Adds the given point and returns the result.
    ///
    /// * `other` - The point to add.
    fn add(self, other: &Point2<T>) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Num> Add<Vector2<T>> for Point2<T> {
    type Output = Self;

    /// Offsets the point by the given vector.
    ///
    /// * `other` - The vector to add.
    fn add(self, other: Vector2<T>) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Num + Copy> Add<&Vector2<T>> for Point2<T> {
    type Output = Self;

    /// Offsets the point by the given vector.
    ///
    /// * `other` - The vector to add.
    fn add(self, other: &Vector2<T>) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y)
    }
}

impl<T: Num + Copy> AddAssign for Point2<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The point to add.
    fn add_assign(&mut self, other: Self) {
        *self = Point2::new(self.x + other.x, self.y + other.y);
    }
}

impl<T: Num + Copy> AddAssign<&Point2<T>> for Point2<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The point to add.
    fn add_assign(&mut self, other: &Point2<T>) {
        *self = Point2::new(self.x + other.x, self.y + other.y);
    }
}

impl<T: Num + Copy> AddAssign<Vector2<T>> for Point2<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The vector to add.
    fn add_assign(&mut self, other: Vector2<T>) {
        *self = Self::new(self.x + other.x, self.y + other.y);
    }
}

impl<T: Num + Copy> AddAssign<&Vector2<T>> for Point2<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` - The vector to add.
    fn add_assign(&mut self, other: &Vector2<T>) {
        *self = Self::new(self.x + other.x, self.y + other.y);
    }
}

impl<T: Num> Sub for Point2<T> {
    type Output = Vector2<T>;

    /// Subtracts the given point and returns the vector towards that point.
    ///
    /// * `other` - The point to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Vector2::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Num + Copy> Sub<&Point2<T>> for Point2<T> {
    type Output = Vector2<T>;

    /// Subtracts the given point and returns the vector towards that point.
    ///
    /// * `other` - The point to subtract.
    fn sub(self, other: &Point2<T>) -> Self::Output {
        Vector2::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Num + Copy> SubAssign<Point2<T>> for Point2<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The point to subtract.
    fn sub_assign(&mut self, other: Point2<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y);
    }
}

impl<T: Num + Copy> SubAssign<&Point2<T>> for Point2<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The point to subtract.
    fn sub_assign(&mut self, other: &Point2<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y);
    }
}

impl<T: Num> Sub<Vector2<T>> for Point2<T> {
    type Output = Self;

    /// Subtracts the given vector and returns the result.
    ///
    /// * `other` - The point to subtract.
    fn sub(self, other: Vector2<T>) -> Self::Output {
        Self::Output::new(self.x - other.x, self.y - other.y)
    }
}

impl<T: Num + Copy> SubAssign<Vector2<T>> for Point2<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The vector to subtract.
    fn sub_assign(&mut self, other: Vector2<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y);
    }
}

impl<T: Num + Copy> SubAssign<&Vector2<T>> for Point2<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` - The vector to subtract.
    fn sub_assign(&mut self, other: &Vector2<T>) {
        *self = Self::new(self.x - other.x, self.y - other.y);
    }
}

impl<T: Num + Copy> Mul<T> for Point2<T> {
    type Output = Self;

    /// Scale the point.
    ///
    /// * `f` - The scaling factor.
    fn mul(self, f: T) -> Self::Output {
        Self::Output::new(f * self.x, f * self.y)
    }
}

macro_rules! premul {
    ($t: ty) => {
        impl Mul<Point2<$t>> for $t {
            type Output = Point2<$t>;
            /// Scale the vector.
            ///
            /// * `p` - The point.
            fn mul(self, p: Point2<$t>) -> Point2<$t> {
                Point2::<$t>::new(self * p.x, self * p.y)
            }
        }

        impl Mul<&Point2<$t>> for $t {
            type Output = Point2<$t>;
            /// Scale the vector.
            ///
            /// * `p` - The point.
            fn mul(self, p: &Point2<$t>) -> Point2<$t> {
                Point2::<$t>::new(self * p.x, self * p.y)
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

impl<T: Num + Copy> MulAssign<T> for Point2<T> {
    /// Scale and assign the result to the point.
    ///
    /// * `f` - The scaling factor.
    fn mul_assign(&mut self, f: T) {
        *self = Self::new(f * self.x, f * self.y);
    }
}

impl<T: Num + Copy> Div<T> for Point2<T> {
    type Output = Self;

    /// Scale the point by 1/f.
    ///
    /// * `f` - The scaling factor.
    fn div(self, f: T) -> Self::Output {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        Self::Output::new(inv * self.x, inv * self.y)
    }
}

impl<T: Num + Copy> DivAssign<T> for Point2<T> {
    /// Scale the point by 1/f and assign the result to the point.
    ///
    /// * `f` - The scaling factor.
    fn div_assign(&mut self, f: T) {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        *self = Self::new(inv * self.x, inv * self.y);
    }
}

impl<T: Num + Neg<Output = T>> Neg for Point2<T> {
    type Output = Self;

    /// Flip the point's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y)
    }
}

impl<T> Index<Axis> for Point2<T> {
    type Output = T;

    /// Index the point by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` - A 2-D coordinate axis.
    fn index(&self, axis: Axis) -> &Self::Output {
        match axis {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            _ => panic!("Invalid axis for std::Index on Point2<T>"),
        }
    }
}

impl<T> Index<usize> for Point2<T> {
    type Output = T;

    /// Index the point by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 2-D coordinate axis.
    fn index(&self, axis: usize) -> &Self::Output {
        &self[Axis::from(axis)]
    }
}

impl<T> IndexMut<Axis> for Point2<T> {
    /// Index the point by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` - A 2-D coordinate axis.
    fn index_mut(&mut self, axis: Axis) -> &mut Self::Output {
        match axis {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            _ => panic!("Invalid axis for std::IndexMut on Point2<T>"),
        }
    }
}

impl<T> IndexMut<usize> for Point2<T> {
    /// Index the point by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 2-D coordinate axis.
    fn index_mut(&mut self, axis: usize) -> &mut Self::Output {
        &mut self[Axis::from(axis)]
    }
}

impl<T> From<Vector2<T>> for Point2<T> {
    /// Convert a 2-D vector to a 2-D point.
    ///
    /// * `v` - 2-D vector.
    fn from(v: Vector2<T>) -> Self {
        Self { x: v.x, y: v.y }
    }
}

impl<T> From<Point3<T>> for Point2<T> {
    /// Convert a 3-D point to a 2-D point by dropping the z-coordinate.
    ///
    /// * `p` - 3-D point.
    fn from(p: Point3<T>) -> Self {
        Self { x: p.x, y: p.y }
    }
}

impl From<Point2i> for Point2f {
    /// Convert a `Point2i` to `Point2f`.
    ///
    /// * `b` - The `Point2i` to convert.
    fn from(p: Point2i) -> Self {
        Self {
            x: p.x as Float,
            y: p.y as Float,
        }
    }
}

impl From<Point2f> for Point2i {
    /// Convert a `Point2f` to `Point2i`.
    ///
    /// * `b` - The `Point2f` to convert.
    fn from(p: Point2f) -> Self {
        Self {
            x: p.x as Int,
            y: p.y as Int,
        }
    }
}

impl<T: fmt::Display> fmt::Display for Point2<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
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
        assert!(Point2::new(0, 0) == Point2::zero());
        assert!(Point2::new(0.0, 0.0) == Point2::zero());
    }

    #[test]
    fn has_nans() {
        assert!(!Point2::new(0.0, 0.0).has_nans());
        assert!(Point2::new(f32::NAN, f32::NAN).has_nans());
        assert!(Point2::new(f64::NAN, f64::NAN).has_nans());
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_i64() {
        Point2::<i64>::zero() / 0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_f64() {
        Point2::<f64>::new(1.0, 1.0) / 0.0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn invalid_index() {
        let z = Point2::<i64>::zero()[Axis::Z];
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn invalid_index_mut() {
        let mut v = Point2::<i64>::zero();
        v[Axis::Z] = 1;
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

    prop_point2!(point2_i32, i32, -100..100i32, -100..100i32);
    prop_point2!(point2_f32, f32, -100.0..100.0f32, -100.0..100.0f32);

    prop_vector2!(vector2_i32, i32, -100..100i32, -100..100i32);
    prop_vector2!(vector2_f32, f32, -100.0..100.0f32, -100.0..100.0f32);

    proptest! {
        #[test]
        fn distance_squared_f32(p1 in point2_f32(), p2 in point2_f32()) {
            let expected = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
            prop_assert_eq!(p1.distance_squared(p2), expected);
        }

        #[test]
        fn length_f32(p1 in point2_f32(), p2 in point2_f32()) {
            let expected = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
            let expected = expected.sqrt();
            prop_assert_eq!(p1.distance(p2), expected);
        }

        #[test]
        fn abs_i32(p in point2_i32()) {
            prop_assert_eq!(p.abs(), Point2::new(abs(p.x), abs(p.y)));
        }

        #[test]
        fn floor_f32(p in point2_f32()) {
            prop_assert_eq!(p.floor(), Point2::new(p.x.floor(), p.y.floor()));
        }

        #[test]
        fn ceil_f32(p in point2_f32()) {
            prop_assert_eq!(p.ceil(), Point2::new(p.x.ceil(), p.y.ceil()));
        }

        #[test]
        fn min_f32(p1 in point2_f32(), p2 in point2_f32()) {
            prop_assert_eq!(p1.min(&p2), Point2::new(p1.x.min(p2.x), p1.y.min(p2.y)));
        }

        #[test]
        fn max_f32(p1 in point2_f32(), p2 in point2_f32()) {
            prop_assert_eq!(p1.max(&p2), Point2::new(p1.x.max(p2.x), p1.y.max(p2.y)));
        }

        #[test]
        fn permute_f32(p in point2_f32(), a1 in axis_2d_strategy(), a2 in axis_2d_strategy()) {
            let permuted_p = p.permute(a1, a2);
            prop_assert_eq!(permuted_p.x, p[a1]);
            prop_assert_eq!(permuted_p.y, p[a2]);
        }

        #[test]
        fn lerp_edge_case_f32(p1 in point2_f32(), p2 in point2_f32()) {
            prop_assert_eq!(lerp(0.0, p1, p2), p1);
            prop_assert_eq!(lerp(1.0, p1, p2), p2);
        }

        #[test]
        fn lerp_f32(p1 in point2_f32(), p2 in point2_f32(), t in -2.0..2.0f32) {
            prop_assert_eq!(lerp(t, p1, p2), (1.0 - t) * p1 + t * p2);
        }

        #[test]
        fn add_point_i32(p1 in point2_i32(), p2 in point2_i32()) {
            prop_assert_eq!(p1 + p2, Point2::new(p1.x + p2.x, p1.y + p2.y));
        }

        #[test]
        fn add_point_f32(p1 in point2_f32(), p2 in point2_f32()) {
            prop_assert_eq!(p1 + p2, Point2::new(p1.x + p2.x, p1.y + p2.y));
        }

        #[test]
        fn add_vector_f32(p in point2_f32(), v in vector2_f32()) {
            prop_assert_eq!(p + v, Point2::new(p.x + v.x, p.y + v.y));
        }

        #[test]
        fn add_assign_vector_i32(p in point2_i32(), v in vector2_i32()) {
            let mut p1 = p;
            p1 += v;
            prop_assert_eq!(p1, Point2::new(p.x + v.x, p.y + v.y));
        }

        #[test]
        fn add_assign_vector_f32(p in point2_f32(), v in vector2_f32()) {
            let mut p1 = p;
            p1 += v;
            prop_assert_eq!(p1, Point2::new(p.x + v.x, p.y + v.y));
        }

        #[test]
        fn sub_point_i32(p1 in point2_i32(), p2 in point2_i32()) {
            prop_assert_eq!(p1 - p2, Vector2::new(p1.x - p2.x, p1.y - p2.y));
        }

        #[test]
        fn sub_point_f32(p1 in point2_f32(), p2 in point2_f32()) {
            prop_assert_eq!(p1 - p2, Vector2::new(p1.x - p2.x, p1.y - p2.y));
        }

        #[test]
        fn sub_vector_f32(p in point2_f32(), v in vector2_f32()) {
            prop_assert_eq!(p - v, Point2::new(p.x - v.x, p.y - v.y));
        }

        #[test]
        fn sub_assign_vector_i32(p in point2_i32(), v in vector2_i32()) {
            let mut p1 = p;
            p1 -= v;
            prop_assert_eq!(p1, Point2::new(p.x - v.x, p.y - v.y));
        }

        #[test]
        fn sub_assign_vector_f32(p in point2_f32(), v in vector2_f32()) {
            let mut p1 = p;
            p1 -= v;
            prop_assert_eq!(p1, Point2::new(p.x - v.x, p.y - v.y));
        }

        #[test]
        fn mul_i32(p in point2_i32(), f in range_i32()) {
            let expected = Point2::new(p.x * f, p.y * f);
            prop_assert_eq!(p * f, expected);
            prop_assert_eq!(f * p, expected);
        }

        #[test]
        fn mul_f32(p in point2_f32(), f in range_f32()) {
            let expected = Point2::new(p.x * f, p.y * f);
            prop_assert_eq!(p * f, expected);
            prop_assert_eq!(f * p, expected);
        }

        #[test]
        fn mul_assign_i32(p in point2_i32(), f in range_i32()) {
            let mut p1 = p;
            p1 *= f;
            prop_assert_eq!(p1, Point2::new(p.x * f, p.y * f));
        }

        #[test]
        fn mul_assign_f32(p in point2_f32(), f in range_f32()) {
            let mut p1 = p;
            p1 *= f;
            prop_assert_eq!(p1, Point2::new(p.x * f, p.y * f));
        }

        #[test]
        fn div_i32(
            p in point2_i32(),
            f in (-100..100i32).prop_filter("non-zero", |x| *x != 0)
        ) {
            let s = 1 / f;
            prop_assert_eq!(p / f, Point2::new(p.x * s, p.y * s));
        }

        #[test]
        fn div_f32(p in point2_f32(), f in non_zero_f32()) {
            let s = 1.0 / f;
            prop_assert_eq!(p / f, Point2::new(p.x * s, p.y * s));
        }

        #[test]
        fn div_assign_i32(p in point2_i32(), f in non_zero_i32()) {
            let mut p1 = p;
            p1 /= f;

            let s = 1 / f;
            prop_assert_eq!(p1, Point2::new(p.x * s, p.y * s));
        }

        #[test]
        fn div_assign_f32(p in point2_f32(), f in non_zero_f32()) {
            let mut p1 = p;
            p1 /= f;

            let s = 1.0 / f;
            prop_assert_eq!(p1, Point2::new(p.x * s, p.y * s));
        }

        #[test]
        fn neg_i32(p in point2_i32()) {
            prop_assert_eq!(-p, Point2::new(-p.x, -p.y));
            prop_assert_eq!(--p, p);
        }

        #[test]
        fn neg_f32(p in point2_f32()) {
            prop_assert_eq!(-p, Point2::new(-p.x, -p.y));
            prop_assert_eq!(--p, p);
        }

        #[test]
        fn index_i32(p in point2_i32()) {
            prop_assert_eq!(p[Axis::X], p.x);
            prop_assert_eq!(p[Axis::Y], p.y);
        }

        #[test]
        fn index_f32(p in point2_f32()) {
            prop_assert_eq!(p[Axis::X], p.x);
            prop_assert_eq!(p[Axis::Y], p.y);
        }

        #[test]
        fn index_mut_i32(p in point2_i32()) {
            let mut p1 = Point2::new(-200, 200);
            p1[Axis::X] = p.x;
            p1[Axis::Y] = p.y;
            prop_assert_eq!(p1, p);
        }

        #[test]
        fn index_mut_f32(p in point2_f32()) {
            let mut p1 = Point2::new(-200.0, 200.0);
            p1[Axis::X] = p.x;
            p1[Axis::Y] = p.y;
            prop_assert_eq!(p1, p);
        }
    }
}
