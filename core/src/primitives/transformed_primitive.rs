//! Object Instancing and Animated Primitives.

use crate::geometry::*;
use crate::light::*;
use crate::material::*;
use crate::primitive::*;
use std::sync::Arc;

/// TransformedPrimitive stores an underlying primitive and animated transform
/// and is used for object instancing and animated transformations.
#[derive(Clone)]
pub struct TransformedPrimitive {
    /// The primitive.
    pub primitive: ArcPrimitive,

    /// The animated transform.
    pub primitive_to_world: AnimatedTransform,
}

impl TransformedPrimitive {
    /// Create a new transformed primitive.
    ///
    /// * `primitive`          - The primitive.
    /// * `primitive_to_world` - The animated transform.
    pub fn new(primitive: ArcPrimitive, primitive_to_world: AnimatedTransform) -> Self {
        Self {
            primitive: Arc::clone(&primitive),
            primitive_to_world,
        }
    }
}

impl Primitive for TransformedPrimitive {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f {
        self.primitive_to_world
            .motion_bounds(&self.primitive.world_bound())
    }

    /// Returns geometric details if a ray intersects the primitive and updates
    /// the t_max parameter of the ray. If there is no intersection, `None` is
    /// returned.
    ///
    /// * `r`                  - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction> {
        let interpolated_prim_to_world = self.primitive_to_world.interpolate(r.time);
        let mut ray = interpolated_prim_to_world.inverse().transform_ray(r);

        let mut it = self.primitive.intersect(&mut ray)?;
        r.t_max = ray.t_max;
        if !interpolated_prim_to_world.is_identity() {
            interpolated_prim_to_world.transform_surface_interaction(&mut it);
        }

        debug_assert!(it.hit.n.dot(&it.shading.n) >= 0.0);

        Some(it)
    }

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    fn intersect_p(&self, r: &Ray) -> bool {
        let interpolated_prim_to_world = self.primitive_to_world.interpolate(r.time);
        let mut ray = interpolated_prim_to_world.inverse().transform_ray(r);
        self.primitive.intersect_p(&mut ray)
    }

    /// Returns a reference to the AreaLight that describes the primitive’s
    /// emission distribution, if the primitive is itself a light source.
    /// If the primitive is not emissive, this method should return `None`.  
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_area_light() on the primitive from the ray-primitive intersection.
    fn get_area_light(&self) -> Option<ArcAreaLight> {
        error!(
            "TransformedPrimitive::get_area_light() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Returns a reference to the material instance assigned to the primitive.
    /// If `None` is returned, ray intersections with the primitive should be
    /// ignored; the primitive only serves to delineate a volume of space for
    /// participating media. This method is also used to check if two rays have
    /// intersected the same object by comparing their Material pointers.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// get_material() on the primitive from the ray-primitive intersection.
    fn get_material(&self) -> Option<ArcMaterial> {
        error!(
            "TransformedPrimitive::get_material() shouldn't be called; \
            should've gone to GeometricPrimitive."
        );
        None
    }

    /// Initializes representations of the light-scattering properties of the
    /// material at the intersection point on the surface.
    ///
    /// *NOTE*: This should never be called. Calling code should directly call
    /// compute_scattering_functions() on the primitive from the ray-primitive
    /// intersection.
    ///
    /// * `_si`                   - The surface interaction at the intersection.
    /// * `_mode`                 - Transport mode.
    /// * `_allow_multiple_lobes` - Allow multiple lobes.
    fn compute_scattering_functions(
        &self,
        _si: &mut SurfaceInteraction,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        error!(
            "TransformedPrimitive::compute_scattering_functions() shouldn't be \
            called; should've gone to GeometricPrimitive."
        );
    }
}
