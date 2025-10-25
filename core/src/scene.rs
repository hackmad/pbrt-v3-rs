//! Scene

use crate::geometry::*;
use crate::interaction::*;
use crate::light::*;
use crate::primitive::*;
use crate::sampler::*;
use crate::spectrum::*;
use crate::{stat_counter, stat_inc, stat_register_fns, stats::*};
use std::collections::HashMap;
use std::sync::Arc;

stat_counter!(
    "Intersections/Regular ray intersection tests",
    N_INTERSECTION_TESTS,
    scene_stats_n_intersection_tests,
);

stat_counter!(
    "Intersections/Shadow ray intersection tests",
    N_SHADOW_TESTS,
    scene_stats_n_shadow_tests,
);

stat_register_fns!(scene_stats_n_intersection_tests, scene_stats_n_shadow_tests);

/// Scene.
pub struct Scene {
    /// An aggregate of all primitives in the scene.
    pub aggregate: ArcPrimitive,

    /// All light sources in the scene.
    pub lights: Vec<ArcLight>,

    /// Infinite light sources in the scene.
    pub infinite_lights: Vec<ArcLight>,

    /// The bounding box of the scene geometry.
    pub world_bound: Bounds3f,

    /// Maps light indices by the Light ID field so we can correctly index into the `lights` vector given a light ID.
    light_id_to_index: HashMap<usize, usize>,
}

impl Scene {
    /// Creates a new `Scene`.
    ///
    /// * `aggregate` - An aggregate of all primitives in the scene.
    /// * `lights`    - All light sources in the scene.
    pub fn new(aggregate: ArcPrimitive, lights: Vec<ArcLight>) -> Self {
        register_stats();

        let mut light_id_to_index = HashMap::new();
        for (i, light) in lights.iter().enumerate() {
            light_id_to_index.insert(light.get_id(), i);
        }

        let scene = Self {
            aggregate: Arc::clone(&aggregate),
            world_bound: aggregate.world_bound(),
            lights: lights.iter().map(Arc::clone).collect(),
            infinite_lights: lights
                .iter()
                .filter(|l| l.get_type().matches(LightType::INFINITE_LIGHT))
                .map(Arc::clone)
                .collect(),
            light_id_to_index,
        };

        for light in lights {
            light.preprocess(&scene);
        }

        scene
    }

    /// Retrieve a light by its ID.
    ///
    /// * `id` - The light ID.
    pub fn light_by_id(&self, id: usize) -> Option<&ArcLight> {
        let index = self.light_id_to_index.get(&id);
        index.map(|i| &self.lights[*i])
    }

    /// Traces the ray into the scene and returns the `SurfaceInteraction` if an intersection occurred.
    ///
    /// * `ray` - The ray to trace.
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction<'_>> {
        stat_inc!(N_INTERSECTION_TESTS, 1);
        self.aggregate.intersect(ray)
    }

    /// Traces the ray into the scene and returns whether or not an intersection occurred.
    ///
    /// * `ray` - The ray to trace.
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        stat_inc!(N_SHADOW_TESTS, 1);
        self.aggregate.intersect_p(ray)
    }

    /// Traces the ray into the scene and returns the first intersection with a light scattering surface along the given
    /// ray as the beam transmittance up to that point.
    ///
    /// * `ray`     - The ray to trace.
    /// * `sampler` - Sampler.
    pub fn intersect_tr(&self, ray: &mut Ray, sampler: &mut dyn Sampler) -> (Option<SurfaceInteraction<'_>>, Spectrum) {
        let mut tr = Spectrum::ONE;

        loop {
            let hit_surface = self.intersect(ray);

            // Accumulate beam transmittance for ray segment
            if let Some(medium) = &ray.medium {
                tr *= medium.tr(ray, sampler);
            }

            // Initialize next ray segment or terminate transmittance computation.
            if let Some(isect) = hit_surface {
                if isect.primitive.unwrap().get_material().is_some() {
                    return (Some(isect), tr);
                }

                *ray = isect.hit.spawn_ray(&ray.d);
            } else {
                return (None, tr);
            }
        }
    }
}
