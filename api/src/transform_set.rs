//! Transform Set

use core::geometry::{ArcTransform, Transform};
use std::ops::{Index, IndexMut};
use std::sync::Arc;

/// Number of transformations to store.
/// NOTE: TransformCache assumes only two transforms. Don't change this.
pub const MAX_TRANSFORMS: usize = 2;

/// Transformation for starting time.
pub const START_TRANSFORM_BITS: usize = 1 << 0;

/// Transformation for ending time.
pub const END_TRANSFORM_BITS: usize = 1 << 1;

/// Transformation for both starting and ending time.
pub const ALL_TRANSFORM_BITS: usize = (1 << MAX_TRANSFORMS) - 1;

/// Stores an array of transformations.
#[derive(Clone, Debug, Default)]
pub struct TransformSet {
    /// The transformations.
    t: [ArcTransform; MAX_TRANSFORMS],
}

impl TransformSet {
    /// Returns a new `TransformSet` containing the inverse transformations.
    pub fn inverse(&self) -> Self {
        let mut t_inv = Self::default();
        for i in 0..self.t.len() {
            t_inv.t[i] = Arc::new(self.t[i].inverse());
        }
        t_inv
    }

    /// Returns `true` if 2 successive transformations are not the same
    /// indicating that this is storing animated transforms.
    pub fn is_animated(&self) -> bool {
        for i in 0..self.t.len() - 1 {
            if self.t[i] != self.t[i + 1] {
                return true;
            }
        }
        false
    }

    /// Reset transforms to identity.
    pub fn reset(&mut self) {
        for i in 0..self.t.len() {
            self.t[i] = Arc::new(Transform::IDENTITY);
        }
    }
}

impl Index<usize> for TransformSet {
    type Output = ArcTransform;

    /// Return the `Transform` at the given index.
    ///
    /// * `index` - The index.
    fn index(&self, index: usize) -> &Self::Output {
        &self.t[index]
    }
}

impl IndexMut<usize> for TransformSet {
    /// Return mutable `Transform` at the given index.
    ///
    /// * `index` - The index.
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.t[index]
    }
}
