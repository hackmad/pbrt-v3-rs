//! Surface Interactions

#![allow(dead_code)]
use super::{hit, ArcShape, Float, Hit, Interaction, Normal3f, Point2f, Point3f, Vector3f};

/// SurfaceInteraction represents geometry of a particular point on a surface.
#[derive(Clone)]
pub struct SurfaceInteraction {
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
}

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
pub fn surface_interaction(
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
) -> SurfaceInteraction {
    // Calculate normal n from the partial derivatives.
    let mut n = Normal3f::from(dpdu.cross(&dpdv).normalize());

    // Adjust normal based on orientation and handedness
    if let Some(s) = shape.clone() {
        if s.get_data().reverse_orientation ^ s.get_data().transform_swaps_handedness {
            n *= -1.0;
        }
    }

    SurfaceInteraction {
        hit: hit(p, time, p_error, wo, n, None),
        uv,
        dpdu,
        dpdv,
        dndu,
        dndv,
        shape: shape.clone(),
        shading: shading(n, dpdu, dpdv, dndu, dndv),
    }
}

impl Interaction for SurfaceInteraction {}

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

/// Create a new shading struct.
/// * `n`    - Surface normal.
/// * `dpdu` - Parametric partial derivative of the point ∂p/∂u.
/// * `dpdv` - Parametric partial derivative of the point ∂p/∂v.
/// * `dndu` - Differential change ∂n/∂v in surface normal as we move along u.
/// * `dndv` - Differential change ∂n/∂v in surface normal as we move along v.
pub fn shading(
    n: Normal3f,
    dpdu: Vector3f,
    dpdv: Vector3f,
    dndu: Normal3f,
    dndv: Normal3f,
) -> Shading {
    Shading {
        n,
        dpdu,
        dpdv,
        dndu,
        dndv,
    }
}
