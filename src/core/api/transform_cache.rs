//! Transform Cache

#![allow(dead_code)]
use crate::core::geometry::ArcTransform;
use std::collections::HashSet;

/// Allocates and stores a single `Transform` reference for each unique
/// transformation.
pub struct TransformCache {
    /// Caches the transformations.
    transforms: HashSet<ArcTransform>,
}

impl TransformCache {
    /// Lookup a reference to a `Transform`. If it is cached, return a
    /// reference to it. Otherwise, insert it and return the cloned reference.
    ///
    /// * `t` - Reference to a transform to lookup.
    pub fn lookup(&mut self, t: ArcTransform) -> ArcTransform {
        match self.transforms.get(&t) {
            Some(transform) => transform.clone(),
            None => {
                self.transforms.insert(t.clone());
                t.clone()
            }
        }
    }

    /// Clear the cached transformations.
    pub fn clear(&mut self) {
        self.transforms.clear();
    }
}

impl Default for TransformCache {
    /// Returns an empty `TransformCache`.
    fn default() -> Self {
        Self {
            transforms: HashSet::new(),
        }
    }
}
