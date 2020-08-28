//! Light

#![allow(dead_code)]
use std::sync::Arc;

/// Light trait provides common behavior.
pub trait Light {}

/// Atomic reference counted `Light`.
pub type ArcLight = Arc<dyn Light + Send + Sync>;

/// AreaLight trait provides common behavior for area lights.
pub trait AreaLight: Light {}

/// Atomic reference counted `AreaLight`.
pub type ArcAreaLight = Arc<dyn AreaLight + Send + Sync>;
