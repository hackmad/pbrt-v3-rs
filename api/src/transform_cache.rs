//! Transform Cache

use core::geometry::{ArcTransform, Transform};
use shared_arena::SharedArena;
use std::collections::HashSet;
use std::mem::MaybeUninit;
use std::sync::Arc;

/// Allocates and stores a single `Transform` reference for each unique transformation.
pub struct TransformCache {
    /// Caches the transformations.
    transforms: HashSet<ArcTransform>,

    /// Arena for memory allocations.
    arena: SharedArena<Transform>,
}

impl TransformCache {
    /// Creates a new `TransformCache`.
    pub fn new() -> Self {
        Self {
            transforms: HashSet::new(),
            arena: SharedArena::new(),
        }
    }

    /// Lookup a reference to a `Transform`. If it is cached, return a reference to it. Otherwise, insert it and return
    /// the cloned reference.
    ///
    /// * `t` - Reference to a transform to lookup.
    pub fn lookup(&mut self, t: &Transform) -> ArcTransform {
        match self.transforms.get(t) {
            Some(transform) => Arc::clone(transform),
            None => {
                let boxed_t = self.arena.alloc_with(|uninit| initialize_data(uninit, t));
                let ret = Arc::new(boxed_t.to_owned());
                self.transforms.insert(Arc::clone(&ret));
                ret
            }
        }
    }

    /// Clear the cached transformations.
    pub fn clear(&mut self) {
        self.transforms.clear();
        self.arena = SharedArena::new(); // no reset. so drop and make new.
    }
}

/// Initialize a Transform for which space is allocated by memory arena.
///
/// * `uninit` - Uninitialized Transform.
/// * `source` - Transform to use for initialization.
fn initialize_data<'a>(uninit: &'a mut MaybeUninit<Transform>, source: &Transform) -> &'a Transform {
    // SAFETY: `uninit` contains `Transform` and source points to `Transform`; so same size.
    //         Memory arena allocated space for a `Transform`.
    unsafe {
        let ptr = uninit.as_mut_ptr();
        std::ptr::copy(source, ptr, 1);
        &*ptr
    }
}
