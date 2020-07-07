//! Common

use super::{abs, Float};
use num_traits::{Num, Zero};
use std::ops;

/// Dot product trait.
pub trait Dot<V> {
    type Output: Num + Zero + ops::Neg<Output = Self::Output> + PartialOrd + Copy;

    /// Returns the dot product.
    ///
    /// * `other` - The other vector/normal.
    fn dot(&self, other: &V) -> Self::Output;

    /// Returns the absolute value of dot product.
    ///
    /// * `other` - The other vector/normal.
    fn abs_dot(&self, other: &V) -> Self::Output {
        abs(self.dot(other))
    }
}

/// FaceForward trait allows pointing vectors in the same hemisphere as
/// another normal/vector.
///
/// NOTE:
/// This is only used by Normal3, so its probably overkill but a nice
/// example of type bounds and generics.
pub trait FaceForward<T, V>
where
    T: Num + Zero + ops::Neg<Output = T> + PartialOrd + Copy,
    Self: Dot<V, Output = T> + ops::Neg<Output = Self> + Sized + Copy,
{
    /// If the vector/normal is not in the same hemisphere as another,
    /// return flipped vector/normal. Otherwise, return itself.
    ///
    /// * `other` - The other vector.
    fn face_forward(&self, other: &V) -> Self {
        if self.dot(other) < T::zero() {
            -*self
        } else {
            *self
        }
    }
}

/// Linearly interpolate between two points for parameters in [0, 1] and
/// extrapolate for parameters outside that interval.
///
/// * `t` - Parameter.
/// * `p0` - Point at t=0.
/// * `p1` - Point at t=1.
pub fn lerp<P>(t: Float, p0: P, p1: P) -> P
where
    Float: ops::Mul<P, Output = P>,
    P: ops::Add<P, Output = P>,
{
    (1.0 - t) * p0 + t * p1
}
