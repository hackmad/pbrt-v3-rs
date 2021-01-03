//! Transform Set

#![allow(dead_code)]
use crate::core::geometry::Transform;
use std::ops::Index;

/// Number of transformations to store.
pub const MAX_TRANSFORMS: usize = 2;

/// Transformation for starting time.
pub const START_TRANSFORM_BITS: usize = 1 << 0;

/// Transformation for ending time.
pub const END_TRANSFORM_BITS: usize = 1 << 1;

/// All transformations.
pub const ALL_TRANSFORM_BITS: usize = (1 << MAX_TRANSFORMS) - 1;

/// Stores an array of transformations.
#[derive(Copy, Clone, Debug, Default)]
pub struct TransformSet {
    /// The transformations.
    t: [Transform; MAX_TRANSFORMS],
}

impl TransformSet {
    /// Returns a new `TransformSet` containing the inverse transformations.
    pub fn inverse(&self) -> Self {
        let mut t_inv = Self::default();
        for i in 0..MAX_TRANSFORMS {
            t_inv.t[i] = self.t[i].inverse();
        }
        t_inv
    }

    /// Returns `true` if 2 successive transformations are not the same
    /// indicating that this is storing animated transforms.
    pub fn is_animated(&self) -> bool {
        for i in 0..MAX_TRANSFORMS - 1 {
            if self.t[i] != self.t[i + 1] {
                return true;
            }
        }
        false
    }
}

impl Index<usize> for TransformSet {
    type Output = Transform;

    /// Return the `Transform` at the given index.
    ///
    /// * `index` -  The index.
    fn index(&self, index: usize) -> &Self::Output {
        &self.t[index]
    }
}
