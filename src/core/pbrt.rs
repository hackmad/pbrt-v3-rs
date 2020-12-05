//! PBRT common stuff

#![allow(dead_code)]

use num_traits::Num;
use std::ops::{Add, Mul, Neg};

/// Use 32-bit precision for floating point numbers.
pub type Float = f32;

/// Default signed integer to 32-bit.
pub type Int = i32;

/// Infinty (∞)
pub const INFINITY: Float = Float::INFINITY;

/// PI (π)
pub const PI: Float = std::f32::consts::PI;

/// PI/2 (π/2)
pub const PI_OVER_TWO: Float = PI * 0.5;

/// PI/4 (π/4)
pub const PI_OVER_FOUR: Float = PI * 0.25;

/// 2*PI (2π)
pub const TWO_PI: Float = PI * 2.0;

/// Machine Epsilon
pub const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;

/// Shadow Epsilon
pub const SHADOW_EPSILON: Float = 0.0001;

/// Axis enumeration
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Axis {
    X = 0,
    Y = 1,
    Z = 2,
}
impl From<usize> for Axis {
    fn from(i: usize) -> Self {
        match i {
            0 => Axis::X,
            1 => Axis::Y,
            2 => Axis::Z,
            _ => panic!("invalid axis value"),
        }
    }
}
impl From<u8> for Axis {
    fn from(i: u8) -> Self {
        match i {
            0 => Axis::X,
            1 => Axis::Y,
            2 => Axis::Z,
            _ => panic!("invalid axis value"),
        }
    }
}
impl Into<u8> for Axis {
    fn into(self) -> u8 {
        match self {
            Axis::X => 0_u8,
            Axis::Y => 1_u8,
            Axis::Z => 2_u8,
        }
    }
}
impl Into<usize> for Axis {
    fn into(self) -> usize {
        match self {
            Axis::X => 0_usize,
            Axis::Y => 1_usize,
            Axis::Z => 2_usize,
        }
    }
}
impl Add<usize> for Axis {
    type Output = Axis;
    fn add(self, i: usize) -> Self::Output {
        Axis::from((self as usize + i) % 3)
    }
}
impl Default for Axis {
    fn default() -> Self {
        Axis::X
    }
}

/// Returns the absolute value of a number.
///
/// * `n` - The number.
pub fn abs<T>(n: T) -> T
where
    T: Num + Neg<Output = T> + PartialOrd + Copy,
{
    if n < T::zero() {
        -n
    } else {
        n
    }
}

/// Returns the minimum of 2 numbers.
///
/// * `a` - First number.
/// * `b` - Second number.
pub fn min<T>(a: T, b: T) -> T
where
    T: Num + PartialOrd + Copy,
{
    if a < b {
        a
    } else {
        b
    }
}

/// Returns the maximum of 2 numbers.
///
/// * `a` - First number.
/// * `b` - Second number.
pub fn max<T>(a: T, b: T) -> T
where
    T: Num + PartialOrd + Copy,
{
    if a > b {
        a
    } else {
        b
    }
}

/// Clamps a value x to [min, max].
///
/// See https://github.com/rust-lang/rust/issues/44095
///
/// * `x` - The number to clamp.
/// * `min` - Minimum value.
/// * `max` - Maximum value.
pub fn clamp<T>(x: T, min: T, max: T) -> T
where
    T: Num + PartialOrd + Copy,
{
    if x < min {
        min
    } else if x > max {
        max
    } else {
        x
    }
}

/// Returns the error bound for adding n terms.
///
/// * `n` - Number of terms
pub fn gamma(n: Int) -> Float {
    (n as Float * MACHINE_EPSILON) / (1.0 - n as Float * MACHINE_EPSILON)
}

/// Linearly interpolate between two points for parameters in [0, 1] and
/// extrapolate for parameters outside that interval.
///
/// * `t` - Parameter.
/// * `p0` - Point at t=0.
/// * `p1` - Point at t=1.
pub fn lerp<P>(t: Float, p0: P, p1: P) -> P
where
    Float: Mul<P, Output = P>,
    P: Add<P, Output = P>,
{
    (1.0 - t) * p0 + t * p1
}

/// Convert a 32-bit floating point value to its constituent bits and
/// return the representation as 32-bit unsigned integer.
///
/// * `f` - The 32-bit floating point number.
pub fn float_to_bits(f: f32) -> u32 {
    let result: u32;
    unsafe {
        let i: u32 = std::mem::transmute_copy(&f);
        result = i;
    }
    result
}

/// Convert the bits of a 32-bit unsigned interger value and return the
/// representation as a 32-bit floating point value.
///
/// * `i` - The 32-bit unsigned interger.
pub fn bits_to_float(i: u32) -> f32 {
    let result: f32;
    unsafe {
        let f: f32 = std::mem::transmute_copy(&i);
        result = f;
    }
    result
}

/// Bump a floating point value up to the next greater representable floating
/// point value.
///
/// * `v` - Floating point value.
pub fn next_float_up(v: Float) -> Float {
    // Handle infinity and negative zero for next_float_up
    if v.is_infinite() && v > 0.0 {
        return v;
    }

    let nv = if v == -0.0 { 0.0 } else { v };

    // Advance v to next higher float
    let mut ui = float_to_bits(nv);
    if nv >= 0.0 {
        ui += 1;
    } else {
        ui -= 1;
    }

    bits_to_float(ui)
}

/// Bump a floating point value up to the next lower representable floating
/// point value.
///
/// * `v` - Floating point value.
pub fn next_float_down(v: Float) -> Float {
    // Handle infinity and positive zero for next_float_down
    if v.is_infinite() && v < 0.0 {
        return v;
    }

    // Advance v to next lower float
    let nv = if v == 0.0 { -0.0 } else { v };
    let mut ui = float_to_bits(v);
    if nv > 0.0 {
        ui -= 1;
    } else {
        ui += 1;
    }

    bits_to_float(ui)
}

/// Emulates the behavior of `upper_bound` but uses a function object to get
/// values at various indices instead of requiring access to an actual array.
/// It is used to bisect arrays that are procedurally generated such as those
/// interpolated from point samples.
///
/// * `size` - Size of array.
/// * `pred` - Function that returns a value at a given index.
pub fn find_interval<Predicate>(size: usize, pred: Predicate) -> usize
where
    Predicate: Fn(usize) -> bool,
{
    let (mut first, mut len) = (0, size);

    while len > 0 {
        let half = len >> 1;
        let middle = first + half;

        // Bisect range based on value of `pred` at `middle`.
        if pred(middle) {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }

    clamp(first - 1, 0, size - 2)
}

/// Return the cosine of an angle.
///
/// * `theta` - The angle in radians.
#[inline]
pub fn cos(theta: Float) -> Float {
    theta.cos()
}

/// Return the sine of an angle.
///
/// * `theta` - The angle in radians.
#[inline]
pub fn sin(theta: Float) -> Float {
    theta.sin()
}

/// Return the tangent of an angle.
///
/// * `theta` - The angle in radians.
#[inline]
pub fn tan(theta: Float) -> Float {
    theta.tan()
}

/// Return the arccosine of an angle.
///
/// * `theta` - The angle in radians.
#[inline]
pub fn acos(theta: Float) -> Float {
    theta.acos()
}

/// Return the arcsine of an angle.
///
/// * `theta` - The angle in radians.
#[inline]
pub fn asin(theta: Float) -> Float {
    theta.asin()
}

/// Computes the arctangent of a number. Return value is in radians in the range
/// [-π/2, π/2];
///
/// * `theta` - The angle in radians.
#[inline]
pub fn atan(theta: Float) -> Float {
    theta.atan()
}

/// Computes the four quadrant arctangent of `y/x`.
///
/// Return values are in the following ranges based on `y` and `x`:
/// * x = 0, y = 0 => 0
/// * x >= 0       => arctan(y/x) -> [-π/2, π/2]
/// * y >= 0       => arctan(y/x) + π -> (π/2, π]
/// * y < 0        =>  arctan(y/x) - π -> (-π, -π/2)
///
/// * `y` - Proportion of y-coordinate.
/// * `x` - Proportion of x-coordinate.
#[inline]
pub fn atan2(y: Float, x: Float) -> Float {
    y.atan2(x)
}

/// Trait to support base-2 logarithm
pub trait Log2<T: Num> {
    /// Returns log base 2 of a value in given type `T`.
    fn log2(self) -> T;
}

impl Log2<u32> for Float {
    /// Returns log base 2 of a value.
    fn log2(self) -> u32 {
        if self < 1.0 {
            0
        } else {
            let bits = float_to_bits(self);
            let r = (bits >> 23) - 127;
            let t = if bits & (1 << 22) == 0 { 0 } else { 1 };
            r + t
        }
    }
}

impl Log2<i32> for u32 {
    /// Returns log base 2 of a value.
    fn log2(self) -> i32 {
        31_i32 - self.leading_zeros() as i32
    }
}

impl Log2<i32> for i32 {
    /// Returns log base 2 of a value.
    fn log2(self) -> i32 {
        (self as u32).log2()
    }
}

impl Log2<i64> for u64 {
    /// Returns log base 2 of a value.
    fn log2(self) -> i64 {
        63_i64 - self.leading_zeros() as i64
    }
}

impl Log2<i64> for i64 {
    /// Returns log base 2 of a value.
    fn log2(self) -> i64 {
        (self as u64).log2()
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
use proptest::prelude::*;

#[cfg(test)]
pub fn axis_2d_strategy() -> impl Strategy<Value = Axis> {
    prop_oneof![Just(Axis::X), Just(Axis::Y)]
}

#[cfg(test)]
pub fn axis_3d_strategy() -> impl Strategy<Value = Axis> {
    prop_oneof![Just(Axis::X), Just(Axis::Y), Just(Axis::Z)]
}
