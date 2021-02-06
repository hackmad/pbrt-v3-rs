//! Light

#![allow(dead_code)]
use crate::core::scene::*;
use std::sync::Arc;

mod light_type;
mod visibility_tester;

/// Light trait provides common behavior.
pub trait Light {
    /// Initialize the light source before rendering begins.
    #[allow(unused)]
    fn preprocess(&mut self, scene: &Scene) {}

    /// Returns the type of light.
    fn get_type(&self) -> LightType;
}

/// Atomic reference counted `Light`.
pub type ArcLight = Arc<dyn Light + Send + Sync>;

pub trait AsLight {
    fn as_light(&self) -> &'static dyn Light;
}

/// AreaLight trait provides common behavior for area lights.
pub trait AreaLight: Light {}

/// Atomic reference counted `AreaLight`.
pub type ArcAreaLight = Arc<dyn AreaLight + Send + Sync>;

// Re-export
pub use light_type::*;
pub use visibility_tester::*;
