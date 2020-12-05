//! Sampler

#![allow(dead_code)]
use std::sync::Arc;

/// Sampler interface.
pub trait Sampler {}

/// Atomic reference counted `Sampler`.
pub type ArcSampler = Arc<dyn Sampler + Send + Sync>;
