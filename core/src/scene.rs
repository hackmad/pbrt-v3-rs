//! Scene

use crate::geometry::*;
use crate::light::*;
use crate::primitive::*;
use crate::sampler::*;
use crate::spectrum::*;
use std::sync::Arc;

/// Scene.
#[derive(Clone)]
pub struct Scene {
    /// An aggregate of all primitives in the scene.
    pub aggregate: ArcPrimitive,

    /// All light sources in the scene.
    pub lights: Vec<ArcLight>,

    /// Infinite light sources in the scene.
    pub infinite_lights: Vec<ArcLight>,

    /// The bounding box of the scene geometry.
    pub world_bound: Bounds3f,
}

impl Scene {
    /// Creates a new `Scene`.
    ///
    /// * `aggregate` - An aggregate of all primitives in the scene.
    /// * `lights`    - All light sources in the scene.
    pub fn new(aggregate: ArcPrimitive, lights: Vec<ArcLight>) -> Self {
        let scene = Self {
            aggregate: Arc::clone(&aggregate),
            world_bound: aggregate.world_bound(),
            lights: lights.iter().map(|l| Arc::clone(&l)).collect(),
            infinite_lights: lights
                .iter()
                .filter(|l| l.get_type().matches(INFINITE_LIGHT))
                .map(|l| Arc::clone(&l))
                .collect(),
        };

        for light in lights {
            light.preprocess(&scene);
        }

        scene
    }

    /// Traces the ray into the scene and returns the `SurfaceInteraction` if
    /// an intersection occurred.
    ///
    /// * `ray` - The ray to trace.
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        self.aggregate.intersect(ray)
    }

    /// Traces the ray into the scene and returns whether or not an intersection
    /// occurred.
    ///
    /// * `ray` - The ray to trace.
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        self.aggregate.intersect_p(ray)
    }

    /// Traces the ray into the scene and returns the first intersection with a
    /// light scattering surface along the given ray as the beam transmittance
    /// up to that point.
    ///
    /// * `ray`     - The ray to trace.
    /// * `sampler` - Sampler.
    pub fn intersect_tr(
        &self,
        ray: &mut Ray,
        sampler: ArcSampler,
    ) -> Option<(SurfaceInteraction, Spectrum)> {
        let mut tr = Spectrum::new(1.0);

        loop {
            let hit_surface = self.intersect(ray);

            // Accumulate beam transmittance for ray segment
            if let Some(medium) = &ray.medium {
                tr *= medium.tr(ray, Arc::clone(&sampler));
            }

            // Initialize next ray segment or terminate transmittance computation.
            if let Some(isect) = hit_surface {
                if isect.primitive.unwrap().get_material().is_some() {
                    return Some((isect, tr));
                }

                *ray = isect.hit.spawn_ray(&ray.d);
            } else {
                return None;
            }
        }
    }
}
