//! Transform Cache

use core::geometry::{ArcTransform, Transform};
use core::{stat_inc, stat_memory_counter, stat_percent, stat_register_fns, stats::*};
use shared_arena::SharedArena;
use std::collections::HashSet;
use std::mem::MaybeUninit;
use std::sync::Arc;

stat_memory_counter!(
    "Memory/TransformCache",
    TRANSFORM_CACHE_BYTES,
    transform_cache_stats_bytes,
);
stat_percent!(
    "Scene/TransformCache hits",
    N_TRANSFORM_CACHE_HITS,
    N_TRANSFORM_CACHE_LOOKUPS,
    transform_cache_stats_hits,
);

stat_register_fns!(transform_cache_stats_bytes, transform_cache_stats_hits);

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
        register_stats();
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
        stat_inc!(N_TRANSFORM_CACHE_LOOKUPS, 1);

        match self.transforms.get(t) {
            Some(transform) => {
                stat_inc!(N_TRANSFORM_CACHE_HITS, 1);
                Arc::clone(transform)
            }
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
        let (arena_used, arena_free) = self.arena.stats();
        let bytes = (arena_used + arena_free) * std::mem::size_of::<Transform>()
            + self.transforms.len() * (std::mem::size_of::<ArcTransform>() + std::mem::size_of::<Transform>());
        stat_inc!(TRANSFORM_CACHE_BYTES, bytes as u64);

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
