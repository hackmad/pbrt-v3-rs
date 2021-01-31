//! Render options

#![allow(dead_code)]
use super::transform_set::*;
use crate::core::camera::*;
use crate::core::integrator::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::scene::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Stores rendering options.
#[derive(Clone)]
pub struct RenderOptions {
    /// Starting time for animation.
    pub transform_start_time: Float,

    /// Ending time for animation.
    pub transform_end_time: Float,

    /// Filter name.
    pub filter_name: String,

    /// Filter parameters.
    pub filter_params: ParamSet,

    /// Film name.
    pub film_name: String,

    /// Film parameters.
    pub film_params: ParamSet,

    /// Sampler name.
    pub sampler_name: String,

    /// Sampler parameters.
    pub sampler_params: ParamSet,

    /// Accelerator name.
    pub accelerator_name: String,

    /// Accelerator parameters.
    pub accelerator_params: ParamSet,

    /// Integrator name.
    pub integrator_name: String,

    /// Integrator parameters.
    pub integrator_params: ParamSet,

    /// Camera name.
    pub camera_name: String,

    /// Camera parameters.
    pub camera_params: ParamSet,

    /// Camera to world transformation.
    pub camera_to_world: TransformSet,

    /// Named media.
    pub named_media: HashMap<String, ArcMedium>,

    /// Lights.
    pub lights: Vec<ArcLight>,

    /// Primitives.
    pub primitives: Vec<ArcPrimitive>,

    /// Instances.
    pub instances: HashMap<String, Vec<ArcPrimitive>>,

    /// Current instance.
    pub current_instance: Option<Vec<ArcPrimitive>>,

    /// Is there scattering media in the scene.
    pub have_scattering_media: bool,
}

impl Default for RenderOptions {
    /// Return a default `RenderOptions` instance.
    fn default() -> Self {
        Self {
            transform_start_time: 0.0,
            transform_end_time: 1.0,
            filter_name: String::from("box"),
            filter_params: ParamSet::new(),
            film_name: String::from("image"),
            film_params: ParamSet::new(),
            sampler_name: String::from("halton"),
            sampler_params: ParamSet::new(),
            accelerator_name: String::from("bvh"),
            accelerator_params: ParamSet::new(),
            integrator_name: String::from("path"),
            integrator_params: ParamSet::new(),
            camera_name: String::from("perspective"),
            camera_params: ParamSet::new(),
            camera_to_world: TransformSet::default(),
            named_media: HashMap::new(),
            lights: vec![],
            primitives: vec![],
            instances: HashMap::new(),
            current_instance: None,
            have_scattering_media: false,
        }
    }
}

impl RenderOptions {
    /// Returns an `Integrator` based on the render options.
    pub fn make_integrator(&self) -> ArcIntegrator {
        todo!();
    }

    /// Returns a `Scene` based on the render options.
    pub fn make_scene(&self) -> Arc<Scene> {
        todo!();
    }

    /// Returns a `Camera` based on the render options.
    pub fn make_camera(&self) -> ArcCamera {
        todo!();
    }
}
