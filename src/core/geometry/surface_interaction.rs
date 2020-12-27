//! Surface Interactions

#![allow(dead_code)]
use crate::core::bssrdf::*;
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::primitive::*;
use crate::core::reflection::*;

/// SurfaceInteraction represents geometry of a particular point on a surface.
#[derive(Clone)]
pub struct SurfaceInteraction<'a> {
    /// The common interaction data.
    pub hit: Hit,

    /// The uv coordinates from surface parametrization.
    pub uv: Point2f,

    /// Parametric partial derivative of the point ∂p/∂u.
    pub dpdu: Vector3f,

    /// Parametric partial derivative of the point ∂p/∂v.
    pub dpdv: Vector3f,

    /// Differential change ∂n/∂u in surface normal as we move along u.
    pub dndu: Normal3f,

    /// Differential change ∂n/∂v in surface normal as we move along v.
    pub dndv: Normal3f,

    /// Shading geometry used for perturbed values.
    pub shading: Shading,

    /// The shape.
    pub shape: Option<ArcShape>,

    /// The BSDF.
    pub bsdf: Option<ArcBSDF>,

    /// The BSSRDF.
    pub bssrdf: Option<ArcBSSRDF>,

    /// The primitive.
    pub primitive: Option<&'a dyn Primitive>,
}

impl<'a> SurfaceInteraction<'a> {
    /// Create a new surface interaction.
    ///
    /// `p`                - Point of interaction.
    /// `p_error`          - Floating point error for ray intersection points.
    /// `uv`               - The uv coordinates from surface parametrization.
    /// `wo`               - The negative ray direction (outgoing direction used
    ///                      when computing lighting at points).
    /// `dpdu`             - Parametric partial derivative of the point ∂p/∂u.
    /// `dpdv`             - Parametric partial derivative of the point ∂p/∂v.
    /// `dndu`             - Differential change ∂n/∂v in surface normal as we move along u.
    /// `dndv`             - Differential change ∂n/∂v in surface normal as we move along v.
    /// `time`             - Time when interaction occurred.
    /// `shape`            - The shape.
    pub fn new(
        p: Point3f,
        p_error: Vector3f,
        uv: Point2f,
        wo: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        shape: Option<ArcShape>,
    ) -> Self {
        // Calculate normal n from the partial derivatives.
        let mut n = Normal3f::from(dpdu.cross(&dpdv).normalize());

        // Adjust normal based on orientation and handedness
        if let Some(s) = shape.clone() {
            if s.get_data().reverse_orientation ^ s.get_data().transform_swaps_handedness {
                n *= -1.0;
            }
        }

        Self {
            hit: Hit::new(p, time, p_error, wo, n, None),
            uv,
            dpdu,
            dpdv,
            dndu,
            dndv,
            shape: shape.clone(),
            shading: Shading::new(n, dpdu, dpdv, dndu, dndv),
            bsdf: None,
            bssrdf: None,
            primitive: None,
        }
    }

    /// Returns updated shading geometry.
    ///
    /// * `dpdu` - Parametric partial derivative of the point ∂p/∂u.
    /// * `dpdv` - Parametric partial derivative of the point ∂p/∂v.
    /// * `dndu` - Differential change ∂n/∂v in surface normal as we move along u.
    /// * `dndv` - Differential change ∂n/∂v in surface normal as we move along v.
    pub fn update_shading_geometry(
        &self,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        orientation_is_authoritative: bool,
    ) -> (Normal3f, Shading) {
        // Compute normal.
        let mut hit_n = self.hit.n;
        let mut shading_n = Normal3::from(dpdu.cross(&dpdv)).normalize();

        if let Some(s) = self.shape.clone() {
            if s.get_data().reverse_orientation ^ s.get_data().transform_swaps_handedness {
                shading_n = -self.shading.n;
                if orientation_is_authoritative {
                    hit_n = hit_n.face_forward(&shading_n.into());
                } else {
                    shading_n = shading_n.face_forward(&hit_n.into());
                }
            }
        }

        // Initialize shading partial derivative values
        (hit_n, Shading::new(shading_n, dpdu, dpdv, dndu, dndv))
    }
}

impl<'a> Interaction for SurfaceInteraction<'a> {}

/// Shading geometry used for perturbed values for bump mapping.
#[derive(Clone)]
pub struct Shading {
    /// Surface normal.
    pub n: Normal3f,

    /// Parametric partial derivative of the point ∂p/∂u.
    pub dpdu: Vector3f,

    /// Parametric partial derivative of the point ∂p/∂v.
    pub dpdv: Vector3f,

    /// Differential change ∂n/∂v in surface normal as we move along u.
    pub dndu: Normal3f,

    /// Differential change ∂n/∂v in surface normal as we move along v.
    pub dndv: Normal3f,
}

impl Shading {
    /// Create a new shading struct.
    /// * `n`    - Surface normal.
    /// * `dpdu` - Parametric partial derivative of the point ∂p/∂u.
    /// * `dpdv` - Parametric partial derivative of the point ∂p/∂v.
    /// * `dndu` - Differential change ∂n/∂v in surface normal as we move along u.
    /// * `dndv` - Differential change ∂n/∂v in surface normal as we move along v.
    pub fn new(
        n: Normal3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
    ) -> Self {
        Self {
            n,
            dpdu,
            dpdv,
            dndu,
            dndv,
        }
    }
}
