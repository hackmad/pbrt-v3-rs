//! Camera

#![allow(dead_code)]
mod environment_camera;
mod orthographic_camera;
mod perspective_camera;
mod realistic_camera;

// Re-export
pub use environment_camera::*;
pub use orthographic_camera::*;
pub use perspective_camera::*;
pub use realistic_camera::*;

use crate::core::film::Film;
use crate::core::geometry::AnimatedTransform;
use crate::core::medium::ArcMedium;
use crate::core::paramset::ParamSet;
use std::sync::Arc;

/// Stores properties for camera creation.
#[derive(Clone)]
pub struct CameraProps {
    /// Parameter set.
    pub params: ParamSet,

    /// Transformation from camera space to world space.
    pub cam2world: AnimatedTransform,

    /// Film.
    pub film: Arc<Film>,

    /// Medium
    pub medium: ArcMedium,
}
