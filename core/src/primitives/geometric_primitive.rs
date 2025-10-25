//! Geometric Primitives

use super::common::*;
use crate::bssrdf::*;
use crate::geometry::*;
use crate::interaction::*;
use crate::light::*;
use crate::material::*;
use crate::medium::*;
use crate::primitive::*;
use crate::reflection::*;
use crate::stat_inc;
use std::sync::Arc;

/// GeometricPrimitive represents a single shape in a scene.
#[derive(Clone)]
pub struct GeometricPrimitive {
    /// The shape.
    pub shape: ArcShape,

    /// The material.
    pub material: Option<ArcMaterial>,

    /// Optional area light that describes emmission characterisitics if it emits light.
    pub area_light: Option<ArcLight>,

    /// Information about the participating media on the inside and outside the primitive.
    pub medium_interface: MediumInterface,
}

impl GeometricPrimitive {
    /// Create a new geometric primitive.
    ///
    /// * `shape`            - The shape.
    /// * `material`         - The material.
    /// * `area_light`       - Optional area light that describes emmission characterisitics if it emits light.
    /// * `medium_interface` - Information about the participating media on the inside and outside the primitive.
    pub fn new(
        shape: ArcShape,
        material: Option<ArcMaterial>,
        area_light: Option<ArcLight>,
        medium_interface: MediumInterface,
    ) -> Self {
        register_stats();

        stat_inc!(PRIMITIVE_MEMORY, std::mem::size_of::<Self>() as u64);

        Self {
            shape,
            material,
            area_light,
            medium_interface,
        }
    }
}

impl Primitive for GeometricPrimitive {
    /// Returns a bounding box in the world space.
    fn world_bound(&self) -> Bounds3f {
        self.shape.world_bound()
    }

    /// Returns geometric details if a ray intersects the primitive and updates the t_max parameter of the ray. If there
    /// is no intersection, `None` is returned.
    ///
    /// * `r` - The ray.
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction<'_>> {
        let shape = Arc::clone(&self.shape);
        //debug!("GeometricPrimitive::intersect(): Shape: {}", shape.get_type());

        if let Some(mut it) = shape.intersect(r, true) {
            r.t_max = it.t;
            it.isect.primitive = Some(self);

            assert!(it.isect.hit.n.dot(&it.isect.shading.n) >= 0.0);

            // Initialize SurfaceInteraction::mediumInterface after Shape intersection.
            let is_medium_transition = self.medium_interface.is_medium_transition();
            it.isect.hit.medium_interface = if is_medium_transition {
                Some(self.medium_interface.clone())
            } else {
                r.medium.as_ref().map(Arc::clone).map(MediumInterface::from)
            };
            Some(it.isect)
        } else {
            None
        }
    }

    /// Returns `true` if a ray-primitive intersection succeeds; otherwise `false`.
    ///
    /// * `r` - The ray.
    fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r, true)
    }

    /// Returns a reference to the AreaLight that describes the primitive’s emission distribution, if the primitive is
    /// itself a light source. If the primitive is not emissive, this method should return `None`.  
    fn get_area_light(&self) -> Option<ArcLight> {
        self.area_light.as_ref().map(Arc::clone)
    }

    /// Returns a reference to the material instance assigned to the primitive. If `None` is returned, ray intersections
    /// with the primitive should be ignored; the primitive only serves to delineate a volume of space for participating
    /// media. This method is also used to check if two rays have intersected the same object by comparing their
    /// `Material` pointers.
    fn get_material(&self) -> Option<ArcMaterial> {
        self.material.as_ref().map(Arc::clone)
    }

    /// Initializes representations of the light-scattering properties of the material at the intersection point on the
    /// surface.
    ///
    /// * `si`                   - The surface interaction at the intersection.
    /// * `mode`                 - Transport mode.
    /// * `allow_multiple_lobes` - Allow multiple lobes.
    /// * `bsdf`                 - The computed BSDF.
    /// * `bssrdf`               - The computed BSSSRDF.
    fn compute_scattering_functions<'scene>(
        &self,
        si: &mut SurfaceInteraction<'scene>,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        bsdf: &mut Option<BSDF>,
        bssrdf: &mut Option<BSDF>,
    ) {
        if let Some(material) = self.material.as_ref() {
            let mut bssrdf_type: Option<BSSRDF> = None;
            material.compute_scattering_functions(si, mode, allow_multiple_lobes, bsdf, &mut bssrdf_type);

            *bssrdf = bssrdf_type.map(|b| b.new(si, Arc::clone(material), mode));
        }
    }
}
