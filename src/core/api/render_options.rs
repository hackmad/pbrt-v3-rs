//! Render options

#![allow(dead_code)]
use super::graphics_state::GraphicsState;
use super::transform_set::*;
use crate::core::camera::*;
use crate::core::integrator::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::paramset::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::scene::*;
use crate::integrators::*;
use crate::{accelerators::BVHAccel, accelerators::SplitMethod};
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

    /// Object instances (each is a collection of primitives).
    pub instances: HashMap<String, Arc<Vec<ArcPrimitive>>>,

    /// Current instance (a collection of primitives).
    pub current_instance: Option<Arc<Vec<ArcPrimitive>>>,

    /// Is there scattering media in the scene.
    pub have_scattering_media: bool,
}

impl RenderOptions {
    /// Creates a new `RenderOptions` instance.
    pub fn new() -> Self {
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

    /// Returns an `Integrator` based on the render options.
    ///
    /// * `gs` - The `GraphicsState`.
    pub fn make_integrator(&self, gs: &GraphicsState) -> Result<ArcIntegrator, String> {
        let camera = self.make_camera(gs);
        let sampler = GraphicsState::make_sampler(
            &self.sampler_name,
            &self.sampler_params,
            camera.get_film_sample_bounds(),
        )?;

        let integrator: Result<ArcIntegrator, String> = match self.integrator_name.as_str() {
            "whitted" => {
                let p = (&self.integrator_params, sampler, camera);
                Ok(Arc::new(WhittedIntegrator::from(p)))
            }
            _ => Err(format!("Integrator '{}' unknown.", self.integrator_name)),
        };

        if integrator.is_ok() {
            if self.have_scattering_media
                && self.integrator_name != "volpath"
                && self.integrator_name != "bdpt"
                && self.integrator_name != "mlt"
            {
                warn!(
                    "Scene has scattering media but '{}' integrator doesn't support 
                volume scattering. Consider using 'volpath', 'bdpt', or 'mlt'.",
                    self.integrator_name
                );
            }

            // Warn if no light sources are defined.
            if self.lights.is_empty() {
                warn!("No light sources defined in scene; rendering a black image.");
            }
        }

        integrator
    }

    /// Returns a `Scene` based on the render options.
    pub fn make_scene(&mut self) -> Arc<Scene> {
        let scene = match GraphicsState::make_accelerator(
            &self.accelerator_name,
            &self.primitives,
            &self.accelerator_params,
        ) {
            Ok(accelerator) => Arc::new(Scene::new(accelerator, self.lights.clone())),
            Err(err) => {
                warn!("Error: {}. Using BVH.", err);
                let accelerator = Arc::new(BVHAccel::new(&self.primitives, 1, SplitMethod::SAH));
                Arc::new(Scene::new(accelerator, self.lights.clone()))
            }
        };
        self.primitives.clear();
        self.lights.clear();
        scene
    }

    /// Returns a `Camera` based on the render options.
    ///
    /// * `gs` - The `GraphicsState`.
    pub fn make_camera(&self, gs: &GraphicsState) -> ArcCamera {
        let filter = match GraphicsState::make_filter(&self.filter_name, &self.filter_params) {
            Ok(f) => f,
            Err(err) => panic!("{}", err),
        };
        let film = match GraphicsState::make_film(&self.film_name, &self.film_params, filter) {
            Ok(f) => f,
            Err(err) => panic!("{}", err),
        };

        let inside_medium =
            gs.current_inside_medium
                .clone()
                .map_or(None, |m| match self.named_media.get(&m) {
                    Some(medium) => Some(medium.clone()),
                    None => None,
                });
        let outside_medium =
            gs.current_outside_medium
                .clone()
                .map_or(None, |m| match self.named_media.get(&m) {
                    Some(medium) => Some(medium.clone()),
                    None => None,
                });

        let medium_interface = MediumInterface::new(inside_medium, outside_medium);

        match gs.make_camera(
            &self.camera_name,
            &self.camera_params,
            &self.camera_to_world,
            self.transform_start_time,
            self.transform_end_time,
            film,
            &medium_interface,
        ) {
            Ok(camera) => camera,
            Err(err) => panic!("{}", err),
        }
    }
}
