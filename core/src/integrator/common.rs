//! Common

use crate::geometry::*;
use crate::interaction::*;
use crate::light::*;
use crate::pbrt::*;
use crate::reflection::*;
use crate::sampler::*;
use crate::sampling::*;
use crate::scene::*;
use crate::spectrum::*;
use std::sync::Arc;

/// Uniformly sample all lights in the scene for direct lighting.
///
/// * `it`              - The intersection information.
/// * `bsdf`            - The BSDF at the intersection.
/// * `scene`           - The scene.
/// * `sampler`         - The sampler.
/// * `n_light_samples` - The number of samples to take for each light.
/// * `handle_media`    - Indicates whether effects of volumetric attenuation should be considered (default to false).
pub fn uniform_sample_all_lights(
    it: &Interaction,
    bsdf: Option<&BSDF>,
    scene: &Scene,
    sampler: &mut dyn Sampler,
    n_light_samples: &[usize],
    handle_media: bool,
) -> Spectrum {
    let mut l = Spectrum::ZERO;

    for (j, light) in scene.lights.iter().enumerate() {
        // Accumulate contribution of j^th light to `l`.
        let n_samples = n_light_samples[j];

        let u_light_array = sampler.get_2d_array(n_samples);
        let u_scattering_array = sampler.get_2d_array(n_samples);

        let nl = u_light_array.len();
        let sl = u_scattering_array.len();

        if nl == 0 || sl == 0 {
            // Use a single sample for illumination from `light`.
            let u_light = sampler.get_2d();
            let u_scattering = sampler.get_2d();

            l += estimate_direct(
                it,
                bsdf,
                &u_scattering,
                Arc::clone(light),
                &u_light,
                scene,
                sampler,
                handle_media,
                false,
            );
        } else {
            // Estimate direct lighting using sample arrays.
            let mut ld = Spectrum::ZERO;
            for k in 0..n_samples {
                ld += estimate_direct(
                    it,
                    bsdf,
                    &u_scattering_array[k],
                    Arc::clone(light),
                    &u_light_array[k],
                    scene,
                    sampler,
                    handle_media,
                    false,
                );
            }
            l += ld / (n_samples as Float);
        }
    }
    l
}

/// Uniformly sample from one random light in the scene for direct lighting and multiply result by number of lights to
/// compensate.
///
/// * `it`            - The intersection information.
/// * `bsdf`          - The BSDF at the intersection.
/// * `scene`         - The scene.
/// * `sampler`       - The sampler.
/// * `handle_media`  - Indicates whether effects of volumetric attenuation should be considered (default to false).
/// * `light_distrib` - PDF for the light's distribution. (default to None).
pub fn uniform_sample_one_light(
    it: &Interaction,
    bsdf: Option<&BSDF>,
    scene: &Scene,
    sampler: &mut dyn Sampler,
    handle_media: bool,
    light_distrib: Option<Arc<Distribution1D>>,
) -> Spectrum {
    // Randomly choose a single light to sample, `light`.
    let n_lights = scene.lights.len();
    if n_lights == 0 {
        return Spectrum::ZERO;
    }

    let (light_num, light_pdf) = if let Some(ld) = light_distrib {
        let sample = sampler.get_1d();
        let (ln, pdf, _) = ld.sample_discrete(sample);
        if pdf == 0.0 {
            return Spectrum::ZERO;
        }
        (ln, pdf)
    } else {
        let sample = sampler.get_1d();
        let ln = min(sample * n_lights as Float, n_lights as Float - 1.0) as usize;
        let pdf = 1.0 / n_lights as Float;
        (ln, pdf)
    };

    let light = Arc::clone(&scene.lights[light_num]);
    let u_light = sampler.get_2d();
    let u_scattering = sampler.get_2d();

    let estimate = estimate_direct(
        it,
        bsdf,
        &u_scattering,
        light,
        &u_light,
        scene,
        sampler,
        handle_media,
        false,
    );
    estimate / light_pdf
}

/// Compute a direct lighting estimate for a light source sample by applying multiple importance sampling.
///
/// * `it`           - The intersection information.
/// * `bsdf`         - The BSDF at the intersection.
/// * `u_scattering` - Scattering sample.
/// * `light`        - The light.
/// * `u_light`      - Light sample.
/// * `scene`        - The scene.
/// * `sampler`      - The sampler.
/// * `handle_media` - Indicates whether effects of volumetric attenuation should be considered (default to false).
/// * `specular`     - Indicates whether perfectly specular lobes should be considered (default to false).
pub fn estimate_direct(
    it: &Interaction,
    bsdf: Option<&BSDF>,
    u_scattering: &Point2f,
    light: ArcLight,
    u_light: &Point2f,
    scene: &Scene,
    sampler: &mut dyn Sampler,
    handle_media: bool,
    specular: bool,
) -> Spectrum {
    let bsdf_flags = if specular {
        BxDFType::all()
    } else {
        BxDFType::all() & !BxDFType::BSDF_SPECULAR
    };
    let mut ld = Spectrum::ZERO;
    let hit = it.get_hit();
    let mut scattering_pdf = 0.0;

    // Sample light source with multiple importance sampling.
    let (mut wi, light_pdf, visibility, mut li) = if let Some(li) = light.sample_li(hit, u_light) {
        (li.wi, li.pdf, Some(li.visibility), li.value)
    } else {
        (Vector3f::ZERO, 0.0, None, Spectrum::ZERO)
    };
    debug!("EstimateDirect uLight: {u_light} -> Li: {li}, wi: {wi}, pdf: {light_pdf}");

    if light_pdf > 0.0 && !li.is_black() {
        let visibility = visibility.unwrap();

        // Compute BSDF or phase function's value for light sample.
        let mut f = Spectrum::ZERO;
        match it {
            Interaction::Surface { si } => {
                // Evaluate BSDF for light sampling strategy.
                if let Some(bsdf) = bsdf.as_ref() {
                    f = bsdf.f(&hit.wo, &wi, bsdf_flags) * wi.abs_dot(&si.shading.n);
                    scattering_pdf = bsdf.pdf(&hit.wo, &wi, bsdf_flags);
                    debug!("  surf f*dot : {f}, scatteringPdf: {scattering_pdf}");
                }
            }
            Interaction::Medium { mi } => {
                // Evaluate phase function for light sampling strategy.
                let p = mi.phase.p(&mi.hit.wo, &wi);
                f = Spectrum::new(p);
                scattering_pdf = p;
                debug!("  medium p: {p}");
            }
            Interaction::Endpoint { ei: _ } => {
                panic!("Endpoint interaction not supported");
            }
        }

        if !f.is_black() {
            // Compute effect of visibility for light source sample.
            if handle_media {
                li *= visibility.tr(scene, sampler);
                debug!("  after Tr, Li: {li}");
            } else if !visibility.unoccluded(scene) {
                debug!("  shadow ray blocked");
                li = Spectrum::ZERO;
            } else {
                debug!("  shadow ray unoccluded");
            }

            // Add light's contribution to reflected radiance.
            if !li.is_black() {
                if light.is_delta_light() {
                    ld += f * li / light_pdf;
                } else {
                    let weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
                    ld += f * li * weight / light_pdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling.
    if !light.is_delta_light() {
        let mut f = Spectrum::ZERO;
        let mut sampled_specular = false;
        match it {
            Interaction::Surface { si } => {
                // Sample scattered direction for surface interactions.
                if let Some(bsdf) = bsdf.as_ref() {
                    let BxDFSample {
                        f: f1,
                        pdf: scatter_pdf,
                        wi: wi2,
                        bxdf_type: sampled_type,
                    } = bsdf.sample_f(&hit.wo, u_scattering, bsdf_flags);
                    scattering_pdf = scatter_pdf;
                    wi = wi2;
                    f = f1 * wi.abs_dot(&si.shading.n);
                    sampled_specular = (sampled_type & BxDFType::BSDF_SPECULAR) != BxDFType::BSDF_NONE;
                }
            }
            Interaction::Medium { mi } => {
                // Sample scattered direction for medium interactions.
                let (p, wi2) = mi.phase.sample_p(&mi.hit.wo, u_scattering);
                f = Spectrum::new(p);
                scattering_pdf = p;
                wi = wi2;
            }
            Interaction::Endpoint { ei: _ } => {
                panic!("Endpoint interaction not supported");
            }
        }
        debug!("  BSDF / phase sampling f: {f}, scattering_pdf: {scattering_pdf}");

        if !f.is_black() && scattering_pdf > 0.0 {
            // Account for light contributions along sampled direction `wi`.
            let mut weight = 1.0;
            if !sampled_specular {
                let light_pdf = light.pdf_li(hit, &wi);
                if light_pdf == 0.0 {
                    return ld;
                }
                weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
            }

            // Find intersection and compute transmittance.
            let mut ray = hit.spawn_ray(&wi);
            let (light_isect, tr) = if handle_media {
                scene.intersect_tr(&mut ray, sampler)
            } else {
                (scene.intersect(&mut ray), Spectrum::ONE)
            };

            // Add light contribution from material sampling.
            let mut li = Spectrum::ZERO;
            if let Some(light_isect) = light_isect {
                if let Some(primitive) = light_isect.primitive {
                    if let Some(area_light) = primitive.get_area_light() {
                        let alt = Arc::as_ptr(&area_light) as *const usize;
                        let lt = Arc::as_ptr(&light) as *const usize;
                        if alt == lt {
                            li = light_isect.le(&(-wi));
                        }
                    }
                }
            } else {
                li = light.le(&ray);
            }

            if !li.is_black() {
                ld += f * li * tr * weight / scattering_pdf;
            }
        }
    }

    ld
}

/// Returns the light power distribution in a scene.
///
/// * `scene` - The scene.
pub fn compute_light_power_distribution(scene: &Scene) -> Option<Distribution1D> {
    if scene.lights.is_empty() {
        None
    } else {
        let light_power: Vec<Float> = scene.lights.iter().map(|light| light.power().y()).collect();
        Some(Distribution1D::new(light_power))
    }
}
