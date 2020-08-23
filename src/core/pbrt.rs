//! PBRT common stuff

#![allow(dead_code)]

use num_traits::Num;
use std::ops::{Add, Neg};

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
impl Add<usize> for Axis {
    type Output = Axis;
    fn add(self, i: usize) -> Self::Output {
        Axis::from((self as usize + i) % 3)
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

/// Returns log base 2 of a value.
///
/// * `v` - The floating point value.
pub fn log2(v: Float) -> u32 {
    if v < 1.0 {
        0
    } else {
        let bits = float_to_bits(v);
        let r = (bits >> 23) - 127;
        let t = if bits & (1 << 22) == 0 { 0 } else { 1 };
        r + t
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
