//! Volumetric Path Integrator

#![allow(dead_code)]

use core::camera::*;
use core::geometry::*;
use core::integrator::*;
use core::interaction::*;
use core::light_distrib::*;
use core::material::*;
use core::medium::MediumInterface;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::sampler::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements volumetric light transport path tracing algorithm.
pub struct VolPathIntegrator {
    /// Common data for sampler integrators.
    pub data: SamplerIntegratorData,

    /// Russian roulette threshold used to terminate path sampling.
    rr_threshold: Float,

    /// Light sampling strategy.
    light_sample_strategy: LightSampleStategy,

    /// Light distribution.
    light_distribution: Option<ArcLightDistribution>,
}

impl VolPathIntegrator {
    /// Create a new `VolPathIntegrator`.
    ///
    /// * `max_depth`              - Maximum recursion depth.
    /// * `camera`                 - The camera.
    /// * `sampler`                - The sampler.
    /// * `pixel_bounds`           - Pixel bounds for the image.
    /// * `rr_threshold`           - Russian roulette threshold used to terminate
    ///                              path sampling. (default to 1.0)
    /// * `light_sample_strategoy` - Light sampling strategy (default to Spatial)
    pub fn new(
        max_depth: usize,
        camera: ArcCamera,
        sampler: ArcSampler,
        pixel_bounds: Bounds2i,
        rr_threshold: Float,
        light_sample_strategy: LightSampleStategy,
    ) -> Self {
        Self {
            data: SamplerIntegratorData::new(max_depth, camera, sampler, pixel_bounds),
            rr_threshold,
            light_sample_strategy,
            light_distribution: None, // Will be set in SamplerIntegrator::preprocess()
        }
    }
}

impl SamplerIntegrator for VolPathIntegrator {
    /// Returns the common data.
    fn get_data(&self) -> &SamplerIntegratorData {
        &self.data
    }
}

impl Integrator for VolPathIntegrator {
    /// Preprocess the scene.
    ///
    /// * `scene` - The scene
    fn preprocess(&mut self, scene: &Scene) {
        self.light_distribution = Some(create_light_sample_distribution(
            self.light_sample_strategy,
            scene,
        ));
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
    fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut ArcSampler,
        _depth: usize,
    ) -> Spectrum {
        let mut l = Spectrum::ZERO;
        let mut beta = Spectrum::ONE;
        let mut specular_bounce = false;

        // Added after book publication: etaScale tracks the accumulated effect
        // of radiance scaling due to rays passing through refractive boundaries
        // (see the derivation on p. 527 of the third edition). We track this
        // value in order to remove it from beta when we apply Russian roulette;
        // this is worthwhile, since it lets us sometimes avoid terminating
        // refracted rays that are about to be refracted back out of a medium and
        // thus have their beta value increased.
        let mut eta_scale: Float = 1.0;

        let mut bounces = 0_usize;
        loop {
            debug!("Path tracer bounce {bounces}, current L = {l}, beta = {beta}");

            // Find closest ray intersection or return background radiance.
            let isect = scene.intersect(ray);

            // Sample the participating medium, if present.
            let mi = if let Some(medium) = ray.medium.as_ref() {
                let (s, mut mi) = medium.sample(&ray, sampler);
                if let Some(m) = mi.as_mut() {
                    let inside = Some(Arc::clone(medium));
                    let outside = Some(Arc::clone(medium));
                    m.hit.medium_interface = Some(MediumInterface::new(inside, outside));
                }
                beta *= s;
                mi
            } else {
                None
            };
            if beta.is_black() {
                break;
            }

            // Handle an interaction with a medium or a surface.
            if let Some(mi) = mi {
                // Terminate path if ray escaped or `max_depth` was reached.
                if bounces >= self.data.max_depth {
                    break;
                }

                // Handle scattering at point in medium for volumetric path tracer.
                let light_distribution = self.light_distribution.as_ref().unwrap();
                let distrib = light_distribution.lookup(&mi.hit.p);

                let it = Interaction::Medium { mi };
                let l_sample = uniform_sample_one_light(&it, None, scene, sampler, true, distrib);
                l += beta * l_sample;

                let sample_2d = {
                    let sampler = Arc::get_mut(sampler).unwrap();
                    sampler.get_2d()
                };
                let wo = -ray.d;

                match it {
                    Interaction::Medium { mi } => {
                        let (_, wi) = mi.phase.sample_p(&wo, &sample_2d);
                        *ray = mi.spawn_ray(&wi);
                    }
                    _ => unreachable!(),
                }

                specular_bounce = false;
            } else {
                // Handle scattering at point on surface for volumetric path tracer.

                // Possibly add emitted light at intersection.
                if bounces == 0 || specular_bounce {
                    // Add emitted light at path vertex or from the environment.
                    if let Some(isect) = isect.as_ref() {
                        l += beta * isect.le(&-ray.d);
                        debug!("Added Le -> L = {l}");
                    } else {
                        for light in scene.infinite_lights.iter() {
                            l += beta * light.le(&ray);
                        }
                        debug!("Added infinite area lights -> L = {l}");
                    }
                }

                // Terminate path if ray escaped or `max_depth` was reached.
                if isect.is_none() || bounces >= self.data.max_depth {
                    break;
                }
                let mut isect = isect.unwrap();

                // Compute scattering functions and skip over medium boundaries.
                let mut bsdf: Option<BSDF> = None;
                let mut bssrdf: Option<BSDF> = None;
                isect.compute_scattering_functions(
                    ray,
                    true,
                    TransportMode::Radiance,
                    &mut bsdf,
                    &mut bssrdf,
                );
                if bsdf.is_none() {
                    debug!("Skipping intersection due to null bsdf");
                    *ray = isect.spawn_ray(&ray.d);
                    continue; // No need to update bounces for skips.
                }
                let bsdf = bsdf.unwrap();

                // Initialize common variables for `VolPathIntegrator`.
                let shading_n = isect.shading.n.clone();

                let light_distribution = self.light_distribution.as_ref().unwrap();
                let distrib = light_distribution.lookup(&isect.hit.p);

                let it = Interaction::Surface { si: isect };

                // Sample illumination from lights to find path contribution. (But
                // skip this for perfectly specular BSDFs).
                let ld = beta
                    * uniform_sample_one_light(&it, Some(&bsdf), scene, sampler, true, distrib);
                debug!("Sampled direct lighting Ld = {ld}");
                assert!(ld.y() >= 0.0);
                l += ld;

                // Sample BSDF to get new path direction.
                let sample_2d = {
                    let sampler = Arc::get_mut(sampler).unwrap();
                    sampler.get_2d()
                };
                let wo = -ray.d;
                let BxDFSample {
                    f,
                    pdf,
                    wi,
                    bxdf_type: flags,
                } = bsdf.sample_f(&wo, &sample_2d, BxDFType::all());
                debug!("Sampled BSDF, f = {f}, pdf = {pdf}");

                if f.is_black() || pdf == 0.0 {
                    break;
                }
                beta *= f * wi.abs_dot(&shading_n) / pdf;
                debug!("Updated beta = {beta}");
                assert!(beta.y() >= 0.0);
                debug_assert!(!beta.y().is_infinite());
                specular_bounce = (flags & BxDFType::BSDF_SPECULAR) > BxDFType::BSDF_NONE;

                if (flags & BxDFType::BSDF_SPECULAR > BxDFType::BSDF_NONE)
                    && (flags & BxDFType::BSDF_TRANSMISSION) > BxDFType::BSDF_NONE
                {
                    let eta = bsdf.eta;
                    // Update the term that tracks radiance scaling for refraction
                    // depending on whether the ray is entering or leaving the medium.
                    eta_scale *= if wo.dot(&it.get_hit().n) > 0.0 {
                        eta * eta
                    } else {
                        1.0 / (eta * eta)
                    };
                }
                *ray = it.spawn_ray(&wi);

                // Account for subsurface scattering, if applicable.
                let bssrdf_bxdf = bssrdf
                    .map(|b| b.bxdfs)
                    .map(|b| {
                        if b.is_empty() {
                            None
                        } else {
                            Some(b[0].clone())
                        }
                    })
                    .flatten();
                if bssrdf_bxdf.is_some()
                    && (flags & BxDFType::BSDF_TRANSMISSION) > BxDFType::BSDF_NONE
                {
                    match bssrdf_bxdf.unwrap() {
                        BxDF::TabulatedBSSRDF(bxdf) => {
                            // Importance sample the BSSRDF.
                            let (sample_1d, sample_2d) = {
                                let sampler = Arc::get_mut(sampler).unwrap();
                                let sample_1d = sampler.get_1d();
                                let sample_2d = sampler.get_2d();
                                (sample_1d, sample_2d)
                            };
                            let SampleS {
                                si: pi,
                                bsdf,
                                value: s,
                                pdf,
                            } = bxdf.sample_s(scene, sample_1d, &sample_2d);
                            debug_assert!(!beta.y().is_infinite());
                            if s.is_black() || pdf == 0.0 {
                                break;
                            }
                            beta *= s / pdf;

                            // Account for the direct subsurface scattering component.
                            let pi = pi.unwrap();
                            let pi_hit_wo = pi.hit.wo.clone();
                            let pi_shading_n = pi.shading.n.clone();
                            let pi_light_pdf = light_distribution.lookup(&pi.hit.p);
                            let pi = Interaction::Surface { si: pi };
                            l += beta
                                * uniform_sample_one_light(
                                    &pi,
                                    bsdf.as_ref(),
                                    scene,
                                    sampler,
                                    true,
                                    pi_light_pdf,
                                );

                            // Account for the indirect subsurface scattering component.
                            let sample_2d = {
                                let sampler = Arc::get_mut(sampler).unwrap();
                                sampler.get_2d()
                            };
                            let BxDFSample {
                                f,
                                pdf,
                                wi,
                                bxdf_type: flags,
                            } = bsdf.map_or(BxDFSample::default(), |b| {
                                b.sample_f(&pi_hit_wo, &sample_2d, BxDFType::all())
                            });
                            if f.is_black() || pdf == 0.0 {
                                break;
                            }
                            beta *= f * wi.abs_dot(&pi_shading_n) / pdf;
                            debug_assert!(!beta.y().is_infinite());
                            specular_bounce =
                                (flags & BxDFType::BSDF_SPECULAR) > BxDFType::BSDF_NONE;
                            *ray = pi.spawn_ray(&wi);
                        }
                        _ => warn!("bssrdf reflection model should be BxDF::*BSSRDF."),
                    }
                }
            }
            // Possibly terminate the path with Russian roulette. Factor out
            // radiance scaling due to refraction in `rr_beta`.
            let rr_beta = beta * eta_scale;
            if rr_beta.max_component_value() < self.rr_threshold && bounces > 3 {
                let q = max(0.05, 1.0 - rr_beta.max_component_value());

                let sampler = Arc::get_mut(sampler).unwrap();
                if sampler.get_1d() < q {
                    break;
                };

                beta /= 1.0 - q;
                debug_assert!(!beta.y().is_infinite());
            }

            bounces += 1;
        }

        l
    }
}

impl From<(&ParamSet, ArcSampler, ArcCamera)> for VolPathIntegrator {
    /// Create a `VolPathIntegrator` from given parameter set and camera.
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

        let rr_threshold = params.find_one_float("rrthreshold", 1.0);

        let lss = params.find_one_string("lightsamplestrategy", "spatial".to_owned());
        let light_sample_strategy = LightSampleStategy::from(lss.as_ref());

        Self::new(
            max_depth,
            Arc::clone(&camera),
            Arc::clone(&sampler),
            pixel_bounds,
            rr_threshold,
            light_sample_strategy,
        )
    }
}
