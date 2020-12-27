//! Primitives

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::light::*;
use crate::core::material::*;
use std::sync::Arc;

/// Primitive trait provide common behavior.
pub trait Primitive {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f;

    /// Returns geometric details if a ray intersects the primitive and updates
    /// the t_max parameter of the ray. If there is no intersection, `None` is
    /// returned.
    ///
    /// * `r`                  - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction>;

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    fn intersect_p(&self, r: &Ray) -> bool;

    /// Returns a reference to the AreaLight that describes the primitive’s
    /// emission distribution, if the primitive is itself a light source.
    /// If the primitive is not emissive, this method should return `None`.  
    fn get_area_light(&self) -> Option<ArcAreaLight>;

    /// Returns a reference to the material instance assigned to the primitive.
    /// If `None` is returned, ray intersections with the primitive should be
    /// ignored; the primitive only serves to delineate a volume of space for
    /// participating media. This method is also used to check if two rays have
    /// intersected the same object by comparing their Material pointers.
    fn get_material(&self) -> Option<ArcMaterial>;

    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Indicates whether the material should use
    ///                            BxDFs that aggregate multiple types of
    ///                            scattering into a single BxDF when such BxDFs
    ///                            are available.
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    );
}

/// Atomic referenced counted `Primitive`.
pub type ArcPrimitive = Arc<dyn Primitive + Send + Sync>;

/// Aggregate trait defines common behaviours for ray intersection accelerators.
pub trait Aggregate: Primitive {}

/// Atomic referenced counted `Aggregate`.
pub type ArcAggregate = Arc<dyn Aggregate + Send + Sync>;
