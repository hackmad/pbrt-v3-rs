//! Whitted Integrator

#![allow(dead_code)]

use crate::core::camera::*;
use crate::core::geometry::*;
use crate::core::integrator::*;
use crate::core::light::*;
use crate::core::material::*;
use crate::core::paramset::*;
use crate::core::reflection::*;
use crate::core::sampler::*;
use crate::core::scene::*;
use crate::core::spectrum::*;
use std::sync::Arc;

/// Implements Whitted's ray tracing algorithm.
pub struct WhittedIntegrator {
    /// The `SamplerIntegratorData`.
    data: SamplerIntegratorData,

    /// Maximum recursion depth.
    max_depth: usize,
}

impl WhittedIntegrator {
    /// Create a new `WhittedIntegrator`.
    ///
    /// * `max_depth`    - Maximum recursion depth.
    /// * `camera`       - The camera.
    /// * `sampler`      - The sampler.
    /// * `pixel_bounds` - Pixel bounds for the image.
    pub fn new(
        max_depth: usize,
        camera: ArcCamera,
        sampler: ArcSampler,
        pixel_bounds: Bounds2i,
    ) -> Self {
        Self {
            data: SamplerIntegratorData::new(camera, sampler, pixel_bounds),
            max_depth,
        }
    }
}

impl SamplerIntegrator for WhittedIntegrator {
    /// Returns the common data.
    fn get_data(&self) -> &SamplerIntegratorData {
        &self.data
    }
}

impl Integrator for WhittedIntegrator {
    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&mut self, _scene: Arc<Scene>) {}

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(
        &self,
        ray: &mut Ray,
        scene: Arc<Scene>,
        sampler: &mut ArcSampler,
        depth: usize,
    ) -> Spectrum {
        let mut l = Spectrum::new(0.0);

        // Find closest ray intersection or return background radiance.
        if let Some(mut isect) = scene.intersect(ray) {
            // Compute emitted and reflected light at ray intersection point.

            // Initialize common variables for Whitted integrator.
            let n = isect.shading.n;
            let wo = isect.hit.wo;

            // Compute scattering functions for surface interaction.
            isect.compute_scattering_functions(ray, false, TransportMode::Radiance);
            if isect.bsdf.is_none() {
                let mut new_ray = isect.hit.spawn_ray(&ray.d);
                return self.li(&mut new_ray, scene.clone(), sampler, depth);
            }

            // Compute emitted light if ray hit an area light source.
            l += isect.le(&wo);

            // Add contribution of each light source.
            for light in scene.lights.iter() {
                let sample = Arc::get_mut(sampler).unwrap().get_2d();
                let Li {
                    wi,
                    pdf,
                    visibility,
                    value: li,
                } = light.sample_li(&isect, &sample);

                if li.is_black() || pdf == 0.0 {
                    continue;
                }

                let f = isect
                    .bsdf
                    .clone()
                    .unwrap()
                    .f(&wo, &wi, BxDFType::from(BSDF_ALL));

                if !f.is_black() && visibility.unoccluded(scene.clone()) {
                    l += f * li * wi.abs_dot(&n) / pdf;
                }
            }
            if depth + 1 < self.max_depth {
                // Trace rays for specular reflection and refraction.
                l += SamplerIntegrator::specular_reflect(
                    self,
                    ray,
                    &isect,
                    scene.clone(),
                    sampler.clone(),
                    depth,
                );
                l += SamplerIntegrator::specular_transmit(
                    self,
                    ray,
                    &isect,
                    scene.clone(),
                    sampler.clone(),
                    depth,
                );
            }
        } else {
            if let Some(rd) = ray.differentials {
                for light in scene.lights.iter() {
                    l += light.le(&rd);
                }
            }
        }

        l
    }
}

impl From<(&ParamSet, ArcSampler, ArcCamera)> for WhittedIntegrator {
    /// Create a `WhittedIntegrator` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcSampler, ArcCamera)) -> Self {
        let (params, sampler, camera) = p;

        let max_depth = params.find_one_int("max_depth", 5) as usize;

        let pb = params.find_int("pixel_bounds");
        let np = pb.len();

        let mut pixel_bounds = camera.get_data().film.get_sample_bounds();
        if np > 0 {
            if np != 4 {
                error!("Expected 4 values for 'pixel_bounds' parameter. Got {}", np);
            } else {
                pixel_bounds = pixel_bounds.intersect(&Bounds2i::new(
                    Point2i::new(pb[0], pb[1]),
                    Point2i::new(pb[2], pb[3]),
                ));
                if pixel_bounds.area() == 0 {
                    error!("Degenerate 'pixel_bounds' specified.");
                }
            }
        }

        Self::new(max_depth, camera.clone(), sampler.clone(), pixel_bounds)
    }
}
