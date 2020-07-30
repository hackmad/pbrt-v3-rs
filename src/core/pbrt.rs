//! PBRT common stuff

#![allow(dead_code)]

use num_traits::Num;
use std::ops;

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

/// Axis enumeration
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Axis {
    X = 0,
    Y = 1,
    Z = 2,
}

/// Returns the absolute value of a number.
///
/// * `n` - The number.
pub fn abs<T>(n: T) -> T
where
    T: Num + ops::Neg<Output = T> + PartialOrd + Copy,
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
