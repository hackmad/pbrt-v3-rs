//! Clamp

use super::{Float, INFINITY};
use num_traits::Num;

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

/// Interface for clamping values.
pub trait Clamp<T: Copy> {
    /// Clamps the values to given [low, high] interval.
    ///
    /// * `low`  - Low value.
    /// * `high` - High value.
    fn clamp(&self, low: T, high: T) -> Self;

    /// Clamps the values to some default [low, high] interval determined by
    /// `T`.
    fn clamp_default(&self) -> Self;
}

impl Clamp<Float> for Float {
    /// Clamps the values to given [low, high] interval.
    ///
    /// * `low`  - Low value.
    /// * `high` - High value.
    fn clamp(&self, low: Float, high: Float) -> Self {
        clamp(*self, low, high)
    }

    /// Clamps the values to [0.0, INFINITY].
    fn clamp_default(&self) -> Self {
        clamp(*self, 0.0, INFINITY)
    }
}
