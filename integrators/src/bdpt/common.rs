//! Common

#![allow(dead_code)]
use core::geometry::*;
use core::interaction::*;
use core::material::*;
use core::pbrt::*;
use core::sampling::Distribution1D;
use core::scene::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Returns the correction term for adjoint BSDF with shading normals.
///
/// * `isect` - The surface interaction.
/// * `wo`    - Outgoing direction.
/// * `wi`    - Incident direction.
/// * `mode`  - Light transport mode.
pub(crate) fn correct_shading_normal(
    isect: &SurfaceInteraction,
    wo: &Vector3f,
    wi: &Vector3f,
    mode: TransportMode,
) -> Float {
    match mode {
        TransportMode::Importance => {
            let num = wo.abs_dot(&isect.shading.n) * wi.abs_dot(&isect.hit.n);
            let denom = wo.abs_dot(&isect.hit.n) * wi.abs_dot(&isect.shading.n);

            // wi is occasionally perpendicular to isect.shading.n; this is fine, but we don't want to
            // return an infinite or NaN value in that case.
            if denom == 0.0 {
                0.0
            } else {
                num / denom
            }
        }
        _ => 1.0,
    }
}

/// Calculates the spatial density of infinite area light endpoints expressed as a probability per
/// unit solid angle while accounting for presence of other infinite area lights. It performs a
/// weighted sum of directional densities of all infinite area lights.
///
/// * `scene`                - The scene.
/// * `light_distr`          - Light probabilities.
/// * `light_to_distr_index` - Map of light distribution indices.
/// * `w`                    - The ray direction.
pub(crate) fn infinite_light_density(
    scene: &Scene,
    light_distr: Arc<Distribution1D>,
    light_to_distr_index: Arc<HashMap<usize, usize>>,
    w: &Vector3f,
) -> Float {
    let mut pdf = 0.0;
    for light in scene.infinite_lights.iter() {
        let light_id = light.get_id();

        let index = light_to_distr_index
            .get(&light_id)
            .expect("Light not found in light_to_distr_index map");

        pdf += scene.lights[light_id].pdf_li(&Hit::default(), &-w) * light_distr.func[*index];
    }
    pdf / (light_distr.func_int * light_distr.count() as Float)
}
