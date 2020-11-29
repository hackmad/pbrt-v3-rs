//! 3-D Vectors

#![allow(dead_code)]
use super::common::*;
use super::{abs, max, min, Axis, Float, Int, Normal3, Point3};
use num_traits::{Num, Zero};
use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

/// A 3-D vector containing numeric values.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vector3<T> {
    /// X-coordinate.
    pub x: T,

    /// Y-coordinate.
    pub y: T,

    /// Y-coordinate.
    pub z: T,
}

/// 3-D vector containing `Float` values.
pub type Vector3f = Vector3<Float>;

/// 3-D vector containing `Int` values.
pub type Vector3i = Vector3<Int>;

impl<T: Num> Vector3<T> {
    /// Creates a new 3-D vector.
    ///
    /// * `x` - X-coordinate.
    /// * `y` - Y-coordinate.
    /// * `z` - Z-coordinate.
    pub fn new(x: T, y: T, z: T) -> Self {
        Self { x, y, z }
    }

    /// Creates a new 3-D zero vector.
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

    /// Returns the square of the vector's length.
    pub fn length_squared(&self) -> T
    where
        T: Mul<Output = T> + Add<Output = T> + Copy,
    {
        self.x * self.x + self.y * self.y + self.z * self.z
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
        Self::new(abs(self.x), abs(self.y), abs(self.z))
    }

    /// Returns the cross product with another vector.
    ///
    /// * `other` - The other vector.
    pub fn cross(&self, other: &Self) -> Self
    where
        T: Copy,
    {
        Self::new(
            (self.y * other.z) - (self.z * other.y),
            (self.z * other.x) - (self.x * other.z),
            (self.x * other.y) - (self.y * other.x),
        )
    }

    /// Returns the smallest coordinate value.
    pub fn min_component(&self) -> T
    where
        T: PartialOrd + Copy,
    {
        if self.x.lt(&self.y) {
            if self.x.lt(&self.z) {
                self.x
            } else {
                self.z
            }
        } else if self.y.lt(&self.z) {
            self.y
        } else {
            self.z
        }
    }

    /// Returns the largest coordinate value.
    pub fn max_component(&self) -> T
    where
        T: PartialOrd + Copy,
    {
        if self.x.gt(&self.y) {
            if self.x.gt(&self.z) {
                self.x
            } else {
                self.z
            }
        } else if self.y.gt(&self.z) {
            self.y
        } else {
            self.z
        }
    }

    /// Returns the axis with largest coordinate value.
    pub fn max_dimension(&self) -> Axis
    where
        T: PartialOrd + Copy,
    {
        if self.x.gt(&self.y) {
            if self.x.gt(&self.z) {
                Axis::X
            } else {
                Axis::Z
            }
        } else if self.y.gt(&self.z) {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    /// Return the component-wise minimum coordinate values with another vector.
    ///
    /// * `other` - The other vector.
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

    /// Return the component-wise maximum coordinate values with another vector.
    ///
    /// * `other` - The other vector.
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

    /// Returns a new vector with permuted coordinates according to given axes.
    ///
    /// * `x` - Axis to use for the x-coordinate of returned vector.
    /// * `y` - Axis to use for the y-coordinate of returned vector.
    /// * `z` - Axis to use for the z-coordinate of returned vector.
    pub fn permute(&self, x: Axis, y: Axis, z: Axis) -> Self
    where
        T: Copy,
    {
        Self::new(self[x], self[y], self[z])
    }
}

impl<T: Num + Neg<Output = T> + PartialOrd + Copy> Dot<Vector3<T>> for Vector3<T> {
    type Output = T;

    /// Returns the dot product with another vector.
    ///
    /// * `other` -  The other vector.
    fn dot(&self, other: &Vector3<T>) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl<T: Num + Neg<Output = T> + PartialOrd + Copy> Dot<Normal3<T>> for Vector3<T> {
    type Output = T;

    /// Returns the dot product with another normal.
    ///
    /// * `other` -  The other normal.
    fn dot(&self, other: &Normal3<T>) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

impl<T: Num> Add for Vector3<T> {
    type Output = Self;

    /// Adds the given vector and returns the result.
    ///
    /// * `other` -  The vector to add.
    fn add(self, other: Self) -> Self::Output {
        Self::Output::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl<T: Num + Copy> AddAssign for Vector3<T> {
    /// Performs the `+=` operation.
    ///
    /// * `other` -  The vector to add.
    fn add_assign(&mut self, other: Self) {
        *self = Self::new(self.x + other.x, self.y + other.y, self.z + other.z);
    }
}

impl<T: Num> Sub for Vector3<T> {
    type Output = Self;

    /// Subtracts the given vector and returns the result.
    ///
    /// * `other` -  The vector to subtract.
    fn sub(self, other: Self) -> Self::Output {
        Self::Output::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl<T: Num + Copy> SubAssign for Vector3<T> {
    /// Performs the `-=` operation.
    ///
    /// * `other` -  The vector to subtract.
    fn sub_assign(&mut self, other: Self) {
        *self = Self::new(self.x - other.x, self.y - other.y, self.z - other.z);
    }
}

impl<T: Num + Copy> Mul<T> for Vector3<T> {
    type Output = Vector3<T>;

    /// Scale the vector.
    ///
    /// * `f` -  The scaling factor.
    fn mul(self, f: T) -> Self::Output {
        Self::Output::new(f * self.x, f * self.y, f * self.z)
    }
}

macro_rules! premul {
    ($t: ty) => {
        impl Mul<Vector3<$t>> for $t {
            type Output = Vector3<$t>;
            /// Scale the vector.
            ///
            /// * `v` -  The vector.
            fn mul(self, v: Vector3<$t>) -> Vector3<$t> {
                Vector3::<$t>::new(self * v.x, self * v.y, self * v.z)
            }
        }

        impl Mul<&Vector3<$t>> for $t {
            type Output = Vector3<$t>;
            /// Scale the vector.
            ///
            /// * `v` -  The vector.
            fn mul(self, v: &Vector3<$t>) -> Vector3<$t> {
                Vector3::<$t>::new(self * v.x, self * v.y, self * v.z)
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

impl<T: Num + Copy> MulAssign<T> for Vector3<T> {
    /// Scale and assign the result to the vector.
    ///
    /// * `f` -  The scaling factor.
    fn mul_assign(&mut self, f: T) {
        *self = Self::new(f * self.x, f * self.y, f * self.z);
    }
}

impl<T: Num + Copy> Div<T> for Vector3<T> {
    type Output = Self;

    /// Scale the vector by 1/f.
    ///
    /// * `f` -  The scaling factor.
    fn div(self, f: T) -> Self::Output {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        Self::Output::new(inv * self.x, inv * self.y, inv * self.z)
    }
}

impl<T: Num + Copy> DivAssign<T> for Vector3<T> {
    /// Scale the vector by 1/f and assign the result to the vector.
    ///
    /// * `f` -  The scaling factor.
    fn div_assign(&mut self, f: T) {
        debug_assert!(!f.is_zero());

        let inv = T::one() / f;
        *self = Self::new(inv * self.x, inv * self.y, inv * self.z);
    }
}

impl<T: Num + Neg<Output = T>> Neg for Vector3<T> {
    type Output = Vector3<T>;

    /// Flip the vector's direction (scale by -1).
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.x, -self.y, -self.z)
    }
}

impl<T> Index<Axis> for Vector3<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index(&self, axis: Axis) -> &Self::Output {
        match axis {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            Axis::Z => &self.z,
        }
    }
}

impl<T> Index<usize> for Vector3<T> {
    type Output = T;

    /// Index the vector by an axis to get the immutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index(&self, axis: usize) -> &Self::Output {
        &self[Axis::from(axis)]
    }
}

impl<T> IndexMut<Axis> for Vector3<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index_mut(&mut self, axis: Axis) -> &mut Self::Output {
        match axis {
            Axis::X => &mut self.x,
            Axis::Y => &mut self.y,
            Axis::Z => &mut self.z,
        }
    }
}

impl<T> IndexMut<usize> for Vector3<T> {
    /// Index the vector by an axis to get a mutable coordinate axis value.
    ///
    /// * `axis` -  A 3-D coordinate axis.
    fn index_mut(&mut self, axis: usize) -> &mut Self::Output {
        &mut self[Axis::from(axis)]
    }
}

impl<T> From<Point3<T>> for Vector3<T> {
    /// Convert a 3-D point to a 3-D vector.
    ///
    /// * `p` -  3-D point.
    fn from(p: Point3<T>) -> Self {
        Self {
            x: p.x,
            y: p.y,
            z: p.z,
        }
    }
}

impl<T> From<Normal3<T>> for Vector3<T> {
    /// Convert a 3-D normal to a 3-D vector.
    ///
    /// * `n` -  3-D normal.
    fn from(n: Normal3<T>) -> Self {
        Self {
            x: n.x,
            y: n.y,
            z: n.z,
        }
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::super::axis_3d_strategy;
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn zero_vector() {
        assert!(Vector3::new(0, 0, 0) == Vector3::zero());
        assert!(Vector3::new(0.0, 0.0, 0.0) == Vector3::zero());
    }

    #[test]
    fn has_nans() {
        assert!(!Vector3::new(0.0, 0.0, 0.0).has_nans());
        assert!(Vector3::new(f32::NAN, f32::NAN, f32::NAN).has_nans());
        assert!(Vector3::new(f64::NAN, f64::NAN, f64::NAN).has_nans());
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f64() {
        Vector3::<f64>::zero().normalize();
    }

    #[test]
    #[should_panic]
    fn normalize_zero_f32() {
        Vector3::<f32>::zero().normalize();
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_i64() {
        Vector3::<i64>::zero() / 0;
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn div_zero_f64() {
        Vector3::<f64>::new(1.0, 1.0, 1.0) / 0.0;
    }

    #[test]
    fn cross_axis_f32() {
        let x_axis = Vector3::new(1.0, 0.0, 0.0);
        let y_axis = Vector3::new(0.0, 1.0, 0.0);
        let z_axis = Vector3::new(0.0, 0.0, 1.0);

        assert!(x_axis.cross(&y_axis) == z_axis);
        assert!(y_axis.cross(&x_axis) == -z_axis);

        assert!(y_axis.cross(&z_axis) == x_axis);
        assert!(z_axis.cross(&y_axis) == -x_axis);

        assert!(z_axis.cross(&x_axis) == y_axis);
        assert!(x_axis.cross(&z_axis) == -y_axis);
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

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
        fn length_squared_i32(v in vector3_i32()) {
            prop_assert_eq!(v.length_squared(), v.x * v.x + v.y * v.y + v.z * v.z);
        }

        #[test]
        fn length_squared_f32(v in vector3_f32()) {
            prop_assert_eq!(v.length_squared(), v.x * v.x + v.y * v.y + v.z * v.z);
        }

        #[test]
        fn length_f32(v in vector3_f32()) {
            prop_assert_eq!(v.length(), (v.x * v.x + v.y * v.y + v.z * v.z).sqrt());
        }

        #[test]
        fn normalize_f32(v in vector3_f32()) {
            // Since we do 1.0 / l in implementation we have to do the
            // same here. Doing Vector3::new(x / l, y / l) will not work
            // for some of the floating point values due to precision
            // errors.
            let f = 1.0 / (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
            prop_assert_eq!(v.normalize(), Vector3::new(v.x * f, v.y * f, v.z * f));
        }

        #[test]
        fn abs_i32(v in vector3_i32()) {
            prop_assert_eq!(v.abs(), Vector3::new(abs(v.x), abs(v.y), abs(v.z)));
        }

        #[test]
        fn abs_f32(v in vector3_f32()) {
            prop_assert_eq!(v.abs(), Vector3::new(abs(v.x), abs(v.y), abs(v.z)));
        }

        #[test]
        fn dot_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1.dot(&v2), v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
        }

        #[test]
        fn abs_dot_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1.abs_dot(&v2), (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z).abs());
        }

        #[test]
        fn cross_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            let expected = Vector3 {
                x: (v1.y * v2.z) - (v1.z * v2.y),
                y: (v1.z * v2.x) - (v1.x * v2.z),
                z: (v1.x * v2.y) - (v1.y * v2.x),
            };

            let c = v1.cross(&v2);
            prop_assert_eq!(c, expected);
        }

        #[test]
        fn cross_zero_f32(v in vector3_f32()) {
            let zero = Vector3::<f32>::zero();
            prop_assert_eq!(zero.cross(&v), zero);
            prop_assert_eq!(v.cross(&zero), zero);
            prop_assert_eq!(v.cross(&v), zero);
        }

        #[test]
        fn min_component_f32(v in vector3_f32()) {
            prop_assert_eq!(v.min_component(), v.x.min(v.y).min(v.z));
        }

        #[test]
        fn max_component_f32(v in vector3_f32()) {
            prop_assert_eq!(v.max_component(), v.x.max(v.y).max(v.z));
        }

        #[test]
        fn max_dimension_f32(v in vector3_f32()) {
            let dim = if v.x > v.y {
                if v.x > v.z  {
                    Axis::X
                } else {
                    Axis::Z
                }
            } else if v.y > v.z {
                Axis::Y
            } else {
                Axis::Z
            };
            prop_assert_eq!(v.max_dimension(), dim);
        }

        #[test]
        fn min_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1.min(&v2), Vector3::new(v1.x.min(v2.x), v1.y.min(v2.y), v1.z.min(v2.z)));
        }

        #[test]
        fn max_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1.max(&v2), Vector3::new(v1.x.max(v2.x), v1.y.max(v2.y), v1.z.max(v2.z)));
        }

        #[test]
        fn permute_f32(
            v in vector3_f32(),
            a1 in axis_3d_strategy(),
            a2 in axis_3d_strategy(),
            a3 in axis_3d_strategy()
        ) {
            let permuted_v = v.permute(a1, a2, a3);
            prop_assert_eq!(permuted_v.x, v[a1]);
            prop_assert_eq!(permuted_v.y, v[a2]);
            prop_assert_eq!(permuted_v.z, v[a3]);
        }

        #[test]
        fn add_i32(v1 in vector3_i32(), v2 in vector3_i32()) {
            prop_assert_eq!(v1 + v2, Vector3::new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
        }

        #[test]
        fn add_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1 + v2, Vector3::new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
        }

        #[test]
        fn add_assign_i32(v1 in vector3_i32(), v2 in vector3_i32()) {
            let mut v = v1;
            v += v2;
            prop_assert_eq!(v, Vector3::new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
        }

        #[test]
        fn add_assign_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            let mut v = v1;
            v += v2;
            prop_assert_eq!(v, Vector3::new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
        }

        #[test]
        fn sub_i32(v1 in vector3_i32(), v2 in vector3_i32()) {
            prop_assert_eq!(v1 - v2, Vector3::new(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
        }

        #[test]
        fn sub_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            prop_assert_eq!(v1 - v2, Vector3::new(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
        }

        #[test]
        fn sub_assign_i32(v1 in vector3_i32(), v2 in vector3_i32()) {
            let mut v = v1;
            v -= v2;
            prop_assert_eq!(v, Vector3::new(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
        }

        #[test]
        fn sub_assign_f32(v1 in vector3_f32(), v2 in vector3_f32()) {
            let mut v = v1;
            v -= v2;
            prop_assert_eq!(v, Vector3::new(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z));
        }

        #[test]
        fn mul_i32(v in vector3_i32(), f in range_i32()) {
            let expected = Vector3::new(v.x * f, v.y * f, v.z * f);
            prop_assert_eq!(v * f, expected);
            prop_assert_eq!(f * v, expected);
        }

        #[test]
        fn mul_f32(v in vector3_f32(), f in range_f32()) {
            let expected = Vector3::new(v.x * f, v.y * f, v.z * f);
            prop_assert_eq!(v * f, expected);
            prop_assert_eq!(f * v, expected);
        }

        #[test]
        fn mul_assign_i32(v in vector3_i32(), f in range_i32()) {
            let mut v1 = v;
            v1 *= f;
            prop_assert_eq!(v1, Vector3::new(v.x * f, v.y * f, v.z * f));
        }

        #[test]
        fn mul_assign_f32(v in vector3_f32(), f in range_f32()) {
            let mut v1 = v;
            v1 *= f;
            prop_assert_eq!(v1, Vector3::new(v.x * f, v.y * f, v.z * f));
        }

        #[test]
        fn div_i32(
            v in vector3_i32(),
            f in (-100..100i32).prop_filter("non-zero", |x| *x != 0)
        ) {
            let s = 1 / f;
            prop_assert_eq!(v / f, Vector3::new(v.x * s, v.y * s, v.z * s));
        }

        #[test]
        fn div_f32(v in vector3_f32(), f in non_zero_f32()) {
            let s = 1.0 / f;
            prop_assert_eq!(v / f, Vector3::new(v.x * s, v.y * s, v.z * s));
        }

        #[test]
        fn div_assign_i32(v in vector3_i32(), f in non_zero_i32()) {
            let mut v1 = v;
            v1 /= f;

            let s = 1 / f;
            prop_assert_eq!(v1, Vector3::new(v.x * s, v.y * s, v.z * s));
        }

        #[test]
        fn div_assign_f32(v in vector3_f32(), f in non_zero_f32()) {
            let mut v1 = v;
            v1 /= f;

            let s = 1.0 / f;
            prop_assert_eq!(v1, Vector3::new(v.x * s, v.y * s, v.z * s));
        }

        #[test]
        fn neg_i32(v in vector3_i32()) {
            prop_assert_eq!(-v, Vector3::new(-v.x, -v.y, -v.z));
            prop_assert_eq!(--v, v);
        }

        #[test]
        fn neg_f32(v in vector3_f32()) {
            prop_assert_eq!(-v, Vector3::new(-v.x, -v.y, -v.z));
            prop_assert_eq!(--v, v);
        }

        #[test]
        fn index_i32(v in vector3_i32()) {
            prop_assert_eq!(v[Axis::X], v.x);
            prop_assert_eq!(v[Axis::Y], v.y);
        }

        #[test]
        fn index_f32(v in vector3_f32()) {
            prop_assert_eq!(v[Axis::X], v.x);
            prop_assert_eq!(v[Axis::Y], v.y);
        }

        #[test]
        fn index_mut_i32(v in vector3_i32()) {
            let mut v1 = Vector3::new(-200, 200, -200);
            v1[Axis::X] = v.x;
            v1[Axis::Y] = v.y;
            v1[Axis::Z] = v.z;
            prop_assert_eq!(v1, v);
        }

        #[test]
        fn index_mut_f32(v in vector3_f32()) {
            let mut v1 = Vector3::new(-200.0, 200.0, -200.0);
            v1[Axis::X] = v.x;
            v1[Axis::Y] = v.y;
            v1[Axis::Z] = v.z;
            prop_assert_eq!(v1, v);
        }
    }
}
