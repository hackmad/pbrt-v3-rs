//! Bidirectional scattering surface reflectance distribution function.

use crate::geometry::*;
use crate::interaction::*;
use crate::material::*;
use crate::pbrt::*;
use bumpalo::Bump;
use std::sync::Arc;

/// Implements a simple BSSRDF implementation that can support general `Shapes`.
pub struct SeparableBSSRDF {
    /// Current outgoing surface interaction hit point.
    pub po_hit: Hit,

    /// Current outgoing surface interaction shading information.
    pub po_shading: Shading,

    /// Index of refraction of the scattering medium.
    pub eta: Float,

    /// Local coordinate frame: normal at hit point.
    pub ns: Normal3f,

    /// Local coordinate frame: shading parametric partial derivative of the
    /// point ∂p/∂u (normalized).
    pub ss: Vector3f,

    /// Local coordinate frame: cross product of `ns` x `ss`.
    pub ts: Vector3f,

    /// The underlying material.
    pub material: ArcMaterial,

    /// Light transport mode.
    pub mode: TransportMode,
}

impl SeparableBSSRDF {
    /// Allocate a new instance of `SeparableBSSRDF`.
    ///
    /// * `po`       - Current outgoing surface interaction.
    /// * `eta`      - Index of refraction of the scattering medium.
    /// * `material` - The material.
    /// * `mode`     - Light transport mode.
    pub fn alloc<'arena>(
        arena: &'arena Bump,
        po: &SurfaceInteraction,
        eta: Float,
        material: ArcMaterial,
        mode: TransportMode,
    ) -> &'arena mut Self {
        let ns = po.shading.n;
        let ss = po.shading.dpdu.normalize();
        let ts = Vector3f::from(&ns.cross(&ss));

        arena.alloc(Self {
            po_hit: po.hit.clone(),
            po_shading: po.shading.clone(),
            eta,
            ns,
            ss,
            ts,
            material,
            mode,
        })
    }

    /// Clone into a newly allocated a new instance of `TabulatedBSSRDF`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut Self {
        arena.alloc(Self {
            po_hit: self.po_hit.clone(),
            po_shading: self.po_shading.clone(),
            eta: self.eta,
            ns: self.ns.clone(),
            ss: self.ss.clone(),
            ts: self.ts.clone(),
            material: Arc::clone(&self.material),
            mode: self.mode,
        })
    }
}

/// Evaluate first moment of the Fresnel reflectance function.
///
/// * `eta` - Index of refraction of the scattering medium.
pub fn fresnel_moment_1(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        0.45966 - 1.73965 * eta + 3.37668 * eta2 - 3.904945 * eta3 + 2.49277 * eta4 - 0.68441 * eta5
    } else {
        -4.61686 + 11.1136 * eta - 10.4646 * eta2 + 5.11455 * eta3 - 1.27198 * eta4 + 0.12746 * eta5
    }
}

/// Evaluate second moment of the Fresnel reflectance function.
///
/// * `eta` - Index of refraction of the scattering medium.
pub fn fresnel_moment_2(eta: Float) -> Float {
    let eta2 = eta * eta;
    let eta3 = eta2 * eta;
    let eta4 = eta3 * eta;
    let eta5 = eta4 * eta;
    if eta < 1.0 {
        0.27614 - 0.87350 * eta + 1.12077 * eta2 - 0.65095 * eta3 + 0.07883 * eta4 + 0.04860 * eta5
    } else {
        let r_eta = 1.0 / eta;
        let r_eta2 = r_eta * r_eta;
        let r_eta3 = r_eta2 * r_eta;
        -547.033 + 45.3087 * r_eta3 - 218.725 * r_eta2 + 458.843 * r_eta + 404.557 * eta
            - 189.519 * eta2
            + 54.9327 * eta3
            - 9.00603 * eta4
            + 0.63942 * eta5
    }
}
