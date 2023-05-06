//! Profile Sample

use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug)]
pub struct ProfileSample {
    pub(super) profiler_state: AtomicU64,
    pub(super) count: AtomicU64,
}

impl Default for ProfileSample {
    fn default() -> Self {
        Self {
            profiler_state: AtomicU64::new(0),
            count: AtomicU64::new(0),
        }
    }
}

impl Clone for ProfileSample {
    fn clone(&self) -> Self {
        Self {
            profiler_state: AtomicU64::new(self.profiler_state.load(Ordering::SeqCst)),
            count: AtomicU64::new(self.count.load(Ordering::SeqCst)),
        }
    }
}
