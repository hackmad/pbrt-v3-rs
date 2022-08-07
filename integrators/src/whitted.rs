//! Whitted Integrator

#![allow(dead_code)]

use core::camera::*;
use core::geometry::*;
use core::integrator::*;
use core::light::*;
use core::material::*;
use core::paramset::*;
use core::reflection::*;
use core::sampler::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements Whitted's ray tracing algorithm.
pub struct WhittedIntegrator {
    /// Common data for sampler integrators.
    pub data: SamplerIntegratorData,
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
            data: SamplerIntegratorData::new(max_depth, camera, sampler, pixel_bounds),
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
    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, _scene: &Scene) {}

    /// Render the scene.
    ///
    /// * `scene` - The scene.
    fn render(&self, scene: &Scene) {
        SamplerIntegrator::render(self, scene);
    }

    /// Returns the incident radiance at the origin of a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `scene`   - The scene.
    /// * `sampler` - The sampler.
    /// * `depth`   - The recursion depth.
    fn li(&self, ray: &mut Ray, scene: &Scene, sampler: &mut ArcSampler, depth: usize) -> Spectrum {
        let mut l = Spectrum::ZERO;

        // Find closest ray intersection or return background radiance.
        if let Some(mut isect) = scene.intersect(ray) {
            // Compute emitted and reflected light at ray intersection point.

            // Compute scattering functions for surface interaction.
            let mut bsdf: Option<BSDF> = None;
            let mut bssrdf: Option<BSDF> = None;
            isect.compute_scattering_functions(
                ray,
                false,
                TransportMode::Radiance,
                &mut bsdf,
                &mut bssrdf,
            );
            if bsdf.is_none() {
                let mut new_ray = isect.hit.spawn_ray(&ray.d);
                return self.li(&mut new_ray, scene, sampler, depth);
            }

            // Initialize common variables for `WhittedIntegrator`.
            let n = isect.shading.n;
            let wo = isect.hit.wo;

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
                } = light.sample_li(&isect.hit, &sample);

                if li.is_black() || pdf == 0.0 {
                    continue;
                }

                let f = bsdf
                    .as_ref()
                    .map_or(Spectrum::ZERO, |bsdf| bsdf.f(&wo, &wi, BxDFType::all()));
                if !f.is_black() {
                    // If no visiblity tester, then unoccluded = true.
                    let unoccluded = visibility.map_or(true, |vis| vis.unoccluded(scene));

                    if unoccluded {
                        l += f * li * wi.abs_dot(&n) / pdf;
                    }
                }
            }
            if depth + 1 < self.data.max_depth {
                // Trace rays for specular reflection and refraction.
                let refl = self.specular_reflect(ray, &isect, &bsdf, scene, sampler, depth);
                let trans = self.specular_transmit(ray, &isect, &bsdf, scene, sampler, depth);
                l += refl + trans;
            }
        } else {
            // Return background radiance.
            for light in scene.lights.iter() {
                l += light.le(&ray);
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

        let max_depth = params.find_one_int("maxdepth", 5) as usize;

        let pb = params.find_int("pixelbounds");
        let np = pb.len();

        let mut pixel_bounds = camera.get_data().film.get_sample_bounds();
        if np > 0 {
            if np != 4 {
                error!("Expected 4 values for 'pixel_bounds' parameter. Got {np}");
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

        Self::new(
            max_depth,
            Arc::clone(&camera),
            Arc::clone(&sampler),
            pixel_bounds,
        )
    }
}
