//! Material

use crate::geometry::*;
use crate::interaction::*;
use crate::pbrt::*;
use crate::texture::*;
use bumpalo::Bump;
use std::sync::Arc;

// Light transport mode enumeration.
#[derive(Copy, Clone, PartialEq)]
pub enum TransportMode {
    /// Indicates incident ray that intersected a point started at the camera.
    Radiance,

    /// Indicates incident ray that intersected a point started at the light
    /// source.
    Importance,
}

/// Material trait provides common behavior.
pub trait Material {
    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `arena`                - The arena for memory allocations.
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available.
    fn compute_scattering_functions<'scene, 'arena>(
        &self,
        arena: &'arena Bump,
        si: &mut SurfaceInteraction<'scene, 'arena>,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    );

    /// Update the normal at the surface interaction using a bump map.
    ///
    /// * `d`  - Bump map.
    /// * `si` - Surface interaction.
    fn bump(&self, d: &ArcTexture<Float>, si: &mut SurfaceInteraction) {
        // Compute offset positions and evaluate displacement texture.
        let mut si_eval_hit = si.hit.clone();

        // Shift `si_eval` `du` in the `u` direction.
        let mut du = 0.5 * (abs(si.der.dudx) + abs(si.der.dudy));

        // The most common reason for du to be zero is for ray that start from
        // light sources, where no differentials are available. In this case,
        // we try to choose a small enough du so that we still get a decently
        // accurate bump value.
        if du == 0.0 {
            du = 0.0005;
        }

        si_eval_hit.p += du * si.shading.dpdu;
        let mut si_eval_uv = si.uv + Vector2f::new(du, 0.0);
        si_eval_hit.n = (Normal3f::from(si.shading.dpdu.cross(&si.shading.dpdv))
            + du * si.der.dndu)
            .normalize();
        let u_displace = d.evaluate(&si_eval_hit, &si_eval_uv, &si.der);

        // Shift `si_eval` `dv` in the `v` direction.
        let mut dv = 0.5 * (abs(si.der.dvdx) + abs(si.der.dvdy));
        if dv == 0.0 {
            dv = 0.0005;
        }
        si_eval_hit.p = si.hit.p + dv * si.shading.dpdv;
        si_eval_uv = si.uv + Vector2f::new(0.0, dv);
        si_eval_hit.n = (Normal3f::from(si.shading.dpdu.cross(&si.shading.dpdv))
            + dv * si.der.dndv)
            .normalize();
        let v_displace = d.evaluate(&si_eval_hit, &si_eval_uv, &si.der);
        let displace = d.evaluate(&si.hit, &si.uv, &si.der);

        // Compute bump-mapped differential geometry.
        let dpdu = si.shading.dpdu
            + (u_displace - displace) / du * Vector3f::from(si.shading.n)
            + displace * Vector3f::from(si.shading.dndu);
        let dpdv = si.shading.dpdv
            + (v_displace - displace) / dv * Vector3f::from(si.shading.n)
            + displace * Vector3f::from(si.shading.dndv);

        si.set_shading_geometry(dpdu, dpdv, si.shading.dndu, si.shading.dndv, false);
    }
}

/// Atomic reference counted `Material`.
pub type ArcMaterial = Arc<dyn Material + Send + Sync>;
