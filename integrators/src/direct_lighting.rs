//! Direct Lighting Integrator

use core::camera::*;
use core::geometry::*;
use core::integrator::*;
use core::interaction::*;
use core::material::*;
use core::paramset::*;
use core::reflection::*;
use core::sampler::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::Arc;

/// Direct light sampling strategy.
#[derive(Copy, Clone, Eq, PartialEq)]
pub enum DirectLightStrategy {
    /// Loops over all of the lights and takes a number of samples based on `n_samples` from each of them, summing the
    /// result.
    UniformSampleAll,

    /// Takes a single sample from just one of the lights, chosen at random.
    UniformSampleOne,
}

/// Implements the direct lighting integrator.
pub struct DirectLightingIntegrator {
    /// Common data for sampler integrators.
    pub data: SamplerIntegratorData,

    /// Direct light sampling strategy.
    strategy: DirectLightStrategy,

    /// Number of samples to use for each light source.
    n_light_samples: Vec<usize>,
}

impl DirectLightingIntegrator {
    /// Create a new `DirectLightingIntegrator`.
    ///
    /// * `strategy`     - Light sampling strategy.
    /// * `max_depth`    - Maximum recursion depth.
    /// * `camera`       - The camera.
    /// * `sampler`      - The sampler.
    /// * `pixel_bounds` - Pixel bounds for the image.
    pub fn new(
        strategy: DirectLightStrategy,
        max_depth: usize,
        camera: ArcCamera,
        sampler: ArcSampler,
        pixel_bounds: Bounds2i,
    ) -> Self {
        Self {
            data: SamplerIntegratorData::new(max_depth, camera, sampler, pixel_bounds),
            strategy,
            n_light_samples: vec![],
        }
    }
}

impl SamplerIntegrator for DirectLightingIntegrator {
    /// Returns the common data.
    fn get_data(&self) -> &SamplerIntegratorData {
        &self.data
    }
}

impl Integrator for DirectLightingIntegrator {
    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, scene: &Scene) {
        if self.strategy == DirectLightStrategy::UniformSampleAll {
            let sampler_mut = Arc::get_mut(&mut self.data.sampler).unwrap();

            // Compute number of samples to use for each light.
            for light in scene.lights.iter() {
                self.n_light_samples
                    .push(sampler_mut.round_count(light.get_num_samples()));
            }

            // Request samples for sampling all lights.
            for _i in 0..self.data.max_depth {
                for j in 0..scene.lights.len() {
                    sampler_mut.request_2d_array(self.n_light_samples[j]);
                    sampler_mut.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }

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
            // Compute scattering functions for surface interaction.
            let mut bsdf: Option<BSDF> = None;
            let mut bssrdf: Option<BSDF> = None;
            isect.compute_scattering_functions(ray, false, TransportMode::Radiance, &mut bsdf, &mut bssrdf);
            if bsdf.is_none() {
                let mut new_ray = isect.hit.spawn_ray(&ray.d);
                return self.li(&mut new_ray, scene, sampler, depth);
            }

            // Initialize common variables for `DirectLightingIntegrator`.
            let wo = isect.hit.wo;

            // Compute emitted light if ray hit an area light source.
            l += isect.le(&wo);

            // Create an `Interaction` for sampling lights.
            let it = Interaction::Surface { si: isect };

            if !scene.lights.is_empty() {
                // Compute direct lighting for `DirectLightingIntegrator`.
                let light_sample = match self.strategy {
                    DirectLightStrategy::UniformSampleAll => {
                        uniform_sample_all_lights(&it, bsdf.as_ref(), scene, sampler, &self.n_light_samples, false)
                    }
                    DirectLightStrategy::UniformSampleOne => {
                        uniform_sample_one_light(&it, bsdf.as_ref(), scene, sampler, false, None)
                    }
                };
                l += light_sample;
            }

            // Need to move the `isect` back out of `it`.
            let isect = match it {
                Interaction::Surface { si } => si,
                _ => unreachable!(), // isect was SurfaceInteraction. So this should not be possible.
            };

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

impl From<(&ParamSet, ArcSampler, ArcCamera)> for DirectLightingIntegrator {
    /// Create a `DirectLightingIntegrator` from given parameter set and camera.
    ///
    /// * `p` - A tuple containing parameter set and camera.
    fn from(p: (&ParamSet, ArcSampler, ArcCamera)) -> Self {
        let (params, sampler, camera) = p;

        let max_depth = params.find_one_int("maxdepth", 5) as usize;

        let st = params.find_one_string("strategy", "all".to_string());
        let strategy = match st.as_ref() {
            "one" => DirectLightStrategy::UniformSampleOne,
            "all" => DirectLightStrategy::UniformSampleAll,
            _ => {
                warn!("Strategy '{}' for direct lighting unknown. Using 'all'.", st);
                DirectLightStrategy::UniformSampleAll
            }
        };

        let pb = params.find_int("pixelbounds");
        let np = pb.len();

        let mut pixel_bounds = camera.get_data().film.get_sample_bounds();
        if np > 0 {
            if np != 4 {
                error!("Expected 4 values for 'pixel_bounds' parameter. Got {np}");
            } else {
                pixel_bounds =
                    pixel_bounds.intersect(&Bounds2i::new(Point2i::new(pb[0], pb[1]), Point2i::new(pb[2], pb[3])));
                if pixel_bounds.area() == 0 {
                    error!("Degenerate 'pixel_bounds' specified.");
                }
            }
        }

        Self::new(
            strategy,
            max_depth,
            Arc::clone(&camera),
            Arc::clone(&sampler),
            pixel_bounds,
        )
    }
}
