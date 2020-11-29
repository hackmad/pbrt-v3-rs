//! Common

use super::abs;
use num_traits::{Num, Zero};
use std::ops::Neg;

/// Dot product trait.
pub trait Dot<V> {
    type Output: Num + Zero + Neg<Output = Self::Output> + PartialOrd + Copy;

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
    T: Num + Zero + Neg<Output = T> + PartialOrd + Copy,
    Self: Dot<V, Output = T> + Neg<Output = Self> + Sized + Copy,
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

/// Union trait allows union between two objects.
pub trait Union<T> {
    /// Return the result of a union with an another object of type `T`.
    ///
    /// * `other` - The other object.
    fn union(&self, other: &T) -> Self;
}

/// Intersect trait allows intersection between objects.
pub trait Intersect<T> {
    /// Return the result of an intersection with an another object of type `T`.
    ///
    /// * `other` - The other object.
    fn intersect(&self, other: &T) -> Self;
}
