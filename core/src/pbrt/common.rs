//! Common

#![allow(dead_code)]

use super::clamp::*;
use num_traits::{Num, Zero};
use std::ops::{Add, Mul, Neg};

/// Use 32-bit precision for floating point numbers.
pub type Float = f32;

/// Default signed integer to 32-bit.
pub type Int = i32;

/// Infinty (∞)
pub const INFINITY: Float = Float::INFINITY;

/// PI (π)
pub const PI: Float = std::f32::consts::PI;

/// 1/PI (1/π)
pub const INV_PI: Float = 1.0 / PI;

/// PI/2 (π/2)
pub const PI_OVER_TWO: Float = PI * 0.5;

/// PI/4 (π/4)
pub const PI_OVER_FOUR: Float = PI * 0.25;

/// 2*PI (2π)
pub const TWO_PI: Float = PI * 2.0;

/// 1/2*PI (1/2π)
pub const INV_TWO_PI: Float = 1.0 / TWO_PI;

/// 4*PI (4π)
pub const FOUR_PI: Float = PI * 4.0;

/// 1/4*PI (1/4π)
pub const INV_FOUR_PI: Float = 1.0 / FOUR_PI;

/// Machine Epsilon
pub const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;

/// Shadow Epsilon
pub const SHADOW_EPSILON: Float = 0.0001;

/// Returns the absolute value of a number.
///
/// * `n` - The number.
#[inline(always)]
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
#[inline(always)]
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
#[inline(always)]
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

/// Computes a mod b (the remainder of a divided by b). This version
/// ensures that modulus of a negative number is zero or positive.
///
/// * `a` - Dividend.
/// * `b` - Divisor.
#[inline(always)]
pub fn rem<T>(a: T, b: T) -> T
where
    T: Num + Zero + PartialOrd + Copy,
{
    let result = a - (a / b) * b;
    if result < T::zero() {
        result + b
    } else {
        result
    }
}

/// Returns the error bound for adding n terms.
///
/// * `n` - Number of terms
#[inline(always)]
pub fn gamma(n: Int) -> Float {
    (n as Float * MACHINE_EPSILON) / (1.0 - n as Float * MACHINE_EPSILON)
}

/// Returns gamma corrected values for use in 8-bit images.
///
/// * `value` - Value to correct.
#[inline(always)]
pub fn gamma_correct(value: Float) -> Float {
    if value <= 0.0031308 {
        12.92 * value
    } else {
        1.055 * value.powf(1.0 / 2.4) - 0.055
    }
}

/// Returns inverse of a gamma corrected value.
///
/// * `value` - The value.
#[inline(always)]
pub fn inv_gamma_correct(value: Float) -> Float {
    if value <= 0.04045 {
        value * 1.0 / 12.92
    } else {
        ((value + 0.055) * 1.0 / 1.055).powf(2.4)
    }
}

/// Linearly interpolate between two points for parameters in [0, 1] and
/// extrapolate for parameters outside that interval.
///
/// * `t` - Parameter.
/// * `p0` - Point at t=0.
/// * `p1` - Point at t=1.
#[inline(always)]
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
    // SAFETY: f32 and u32 have same size (32-bits).
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
    // SAFETY: f32 and u32 have same size (32-bits).
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
    let mut ui = float_to_bits(nv);
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
#[inline(always)]
pub fn cos(theta: Float) -> Float {
    theta.cos()
}

/// Return the sine of an angle.
///
/// * `theta` - The angle in radians.
#[inline(always)]
pub fn sin(theta: Float) -> Float {
    theta.sin()
}

/// Return the tangent of an angle.
///
/// * `theta` - The angle in radians.
#[inline(always)]
pub fn tan(theta: Float) -> Float {
    theta.tan()
}

/// Return the arccosine of an angle.
///
/// * `theta` - The angle in radians.
#[inline(always)]
pub fn acos(theta: Float) -> Float {
    theta.acos()
}

/// Return the arcsine of an angle.
///
/// * `theta` - The angle in radians.
#[inline(always)]
pub fn asin(theta: Float) -> Float {
    theta.asin()
}

/// Computes the arctangent of a number. Return value is in radians in the range
/// [-π/2, π/2];
///
/// * `theta` - The angle in radians.
#[inline(always)]
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
#[inline(always)]
pub fn atan2(y: Float, x: Float) -> Float {
    y.atan2(x)
}

/// Returns `v^5`.
///
/// * `v` - The value.
#[inline(always)]
pub fn pow5<T: Mul<T, Output = T> + Copy>(v: T) -> T {
    (v * v) * (v * v) * v
}

/// Returns the error function for a given floating point value.
///
/// * `x` - The floating point value.
#[inline(always)]
pub fn erf(x: Float) -> Float {
    // constants
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    // Save the sign of x
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = abs(x);

    // A&S formula 7.1.26.
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

/// Returns the inverse of the error function for a given floating point value.
///
/// * `x` - The floating point value.
#[inline(always)]
pub fn erf_inv(x: Float) -> Float {
    let x = clamp(x, -0.99999, 0.99999);
    let mut w = -((1.0 - x) * (1.0 + x)).ln();
    if w < 5.0 {
        w -= 2.5;

        let mut p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
        p * x
    } else {
        w = w.sqrt() - 3.0;

        let mut p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
        p * x
    }
}

/// Deconstructs a floating point value as a (mantissa, exponent, sign) tuple
/// so it can be used in data structures that require hashing.
#[derive(Hash, Eq, PartialEq)]
pub struct HashableFloat(u64, i16, i8);
impl HashableFloat {
    /// Create a new `HashableFloat` from a floating point value.
    ///
    /// * `val` - The floating point value.
    pub fn new(val: Float) -> Self {
        let (mantissa, exponent, sign) = num_traits::Float::integer_decode(val);
        Self(mantissa, exponent, sign)
    }
}
