//! Shapes

#![allow(dead_code)]
use super::{ArcTransform, Bounds3f, Float, Ray, SurfaceInteraction};
use std::sync::Arc;

/// Shape common functions
pub trait Shape {
    /// Returns the underlying shape data.
    fn get_data(&self) -> ShapeData;

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f;

    /// Returns a bounding box in the world space.
    ///
    /// Default is to transform the object bounds with the object-to0world
    /// transformation. Override for tighter bounds implementation.
    fn world_bound(&self) -> Bounds3f {
        self.get_data()
            .object_to_world
            .clone()
            .transform_bounds(&self.object_bound())
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests.
    fn intersect<'a>(&self, r: &Ray, test_alpha_texture: bool) -> Option<Intersection<'a>>;

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_textu re` - Perform alpha texture tests.
    fn intersect_p(&self, r: &Ray, test_alpha_texture: bool) -> bool {
        !self.intersect(r, test_alpha_texture).is_none()
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float;
}

/// Atomic reference counted `Shape`.
pub type ArcShape = Arc<dyn Shape + Send + Sync>;

/// Stores geometric information about a single ray-shape intersection.
#[derive(Clone)]
pub struct Intersection<'a> {
    /// The parameter along the ray where intersection occurred.
    pub t: Float,

    /// The surface interaction details.
    pub isect: SurfaceInteraction<'a>,
}

impl<'a> Intersection<'a> {
    /// Create a new intersection.
    ///
    /// * `t`     - The parameter along the ray where intersection occurred.
    /// * `isect` - The surface interaction details.
    pub fn new(t: Float, isect: SurfaceInteraction<'a>) -> Self {
        Self { t, isect }
    }
}

/// Store common shape data.
#[derive(Clone)]
pub struct ShapeData {
    /// The object to world transfomation.
    pub object_to_world: ArcTransform,

    /// The world to object transfomation.
    pub world_to_object: Option<ArcTransform>,

    /// Indicates whether their surface normal directions should be reversed
    /// from the default
    pub reverse_orientation: bool,

    /// Indicates if `object_to_world` transformation changes the handedness
    /// of the coordinate system.
    pub transform_swaps_handedness: bool,
}

impl ShapeData {
    /// Create a new instance of shape data.
    ///
    /// * `object_to_world`     - The object to world transfomation.
    /// * `world_to_object`     - The world to object transfomation.
    /// * `reverse_orientation` - Indicates whether their surface normal directions
    ///                           should be reversed from the default
    pub fn new(
        object_to_world: ArcTransform,
        world_to_object: Option<ArcTransform>,
        reverse_orientation: bool,
    ) -> Self {
        Self {
            object_to_world: object_to_world.clone(),
            world_to_object: world_to_object.clone(),
            reverse_orientation,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
        }
    }
}
