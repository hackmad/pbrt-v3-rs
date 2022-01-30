//! Shapes

#![allow(dead_code)]
use crate::geometry::*;
use crate::interaction::*;
use crate::low_discrepency::radical_inverse;
use crate::pbrt::*;
use std::sync::Arc;

/// Shape common functions
pub trait Shape {
    /// Returns the shape type. Usually these are behind ArcShape and harder to
    /// debug. So this will be helpful.
    fn get_type(&self) -> &'static str;

    /// Returns the underlying shape data.
    fn get_data(&self) -> Arc<ShapeData>;

    /// Returns a bounding box in the shapes object space.
    fn object_bound(&self) -> Bounds3f;

    /// Returns a bounding box in the world space.
    ///
    /// Default is to transform the object bounds with the object-to0world
    /// transformation. Override for tighter bounds implementation.
    fn world_bound(&self) -> Bounds3f {
        Arc::clone(&self.get_data().object_to_world).transform_bounds(&self.object_bound())
    }

    /// Returns geometric details if a ray intersects the shape intersection.
    /// If there is no intersection, `None` is returned.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests.
    fn intersect<'primitive, 'arena>(
        &self,
        r: &Ray,
        test_alpha_texture: bool,
    ) -> Option<Intersection<'primitive, 'arena>>;

    /// Returns `true` if a ray-shape intersection succeeds; otherwise `false`.
    ///
    /// * `r`                  - The ray.
    /// * `test_alpha_texture` - Perform alpha texture tests; default to true.
    fn intersect_p(&self, r: &Ray, test_alpha_texture: bool) -> bool {
        self.intersect(r, test_alpha_texture).is_some()
    }

    /// Returns the surface area of the shape in object space.
    fn area(&self) -> Float;

    /// Sample a point on the surface and return the PDF with respect to area on
    /// the surface.
    ///
    /// NOTE: The returned `Hit` value will have `wo` = Vector3f::ZERO.
    ///
    /// * `u` - Sample value to use.
    fn sample_area(&self, u: &Point2f) -> (Hit, Float);

    /// Sample a point on the shape given a reference point and return the PDF
    /// with respect to the solid angle from ref.
    ///
    /// * `hit` - Reference point on shape.
    /// * `u`   - Sample value to use.
    fn sample_solid_angle(&self, hit: &Hit, u: &Point2f) -> (Hit, Float) {
        let (intr, mut pdf) = self.sample_area(u);
        let mut wi = intr.p - hit.p;

        if wi.length_squared() == 0.0 {
            pdf = 0.0;
        } else {
            wi = wi.normalize();
            // Convert from area measure, as returned by the sample_area() call
            // above, to solid angle measure.
            pdf *= hit.p.distance_squared(intr.p) / intr.n.abs_dot(&(-wi));
            if pdf.is_infinite() {
                pdf = 0.0;
            }
        }

        (intr, pdf)
    }

    /// Return the PDF for the shape. By default it is 1/area.
    ///
    /// * `hit` - The interaction hit point.
    fn pdf(&self, _hit: &Hit) -> Float {
        1.0 / self.area()
    }

    /// Returns the PDF with respect to solid angle.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_solid_angle(&self, hit: &Hit, wi: &Vector3f) -> Float {
        // Intersect sample ray with area light geometry.
        let ray = hit.spawn_ray(wi);

        // Ignore any alpha textures used for trimming the shape when performing
        // this intersection. Hack for the "San Miguel" scene, where this is used
        // to make an invisible area light.
        if let Some(Intersection {
            t: _t_hit,
            isect: isect_light,
        }) = self.intersect(&ray, false)
        {
            // Convert light sample weight to solid angle measure.
            let pdf = hit.p.distance_squared(isect_light.hit.p)
                / (isect_light.hit.n.abs_dot(&(-*wi)) * self.area());
            if pdf.is_infinite() {
                0.0
            } else {
                pdf
            }
        } else {
            0.0
        }
    }

    /// Returns the solid angle subtended by the shape w.r.t. the reference
    /// point p, given in world space. Some shapes compute this value in
    /// closed-form, while the default implementation uses Monte Carlo
    /// integration.
    ///
    /// * `p`         - The reference point.
    /// * `n_samples` - The number of samples to use for Monte-Carlo integration.
    ///                 Default to 512.
    fn solid_angle(&self, p: &Point3f, n_samples: usize) -> Float {
        let hit = Hit::new(
            *p,
            0.0,
            Vector3f::ZERO,
            Vector3f::new(0.0, 0.0, 1.0),
            Normal3f::ZERO,
            None,
        );

        let mut solid_angle: f64 = 0.0;

        for i in 0..n_samples {
            let u = Point2f::new(radical_inverse(0, i as u64), radical_inverse(1, i as u64));
            let (p_shape, pdf) = self.sample_solid_angle(&hit, &u);
            let ray = Ray::new(*p, p_shape.p - *p, 0.999, 0.0, None);
            if pdf > 0.0 && !self.intersect_p(&ray, true) {
                solid_angle += 1.0_f64 / pdf as f64;
            }
        }
        (solid_angle / n_samples as f64) as Float
    }
}

/// Atomic reference counted `Shape`.
pub type ArcShape = Arc<dyn Shape + Send + Sync>;

/// Stores geometric information about a single ray-shape intersection.
pub struct Intersection<'primitive, 'arena> {
    /// The parameter along the ray where intersection occurred.
    pub t: Float,

    /// The surface interaction details.
    pub isect: SurfaceInteraction<'primitive, 'arena>,
}

impl<'primitive, 'arena> Intersection<'primitive, 'arena> {
    /// Create a new intersection.
    ///
    /// * `t`     - The parameter along the ray where intersection occurred.
    /// * `isect` - The surface interaction details.
    pub fn new(t: Float, isect: SurfaceInteraction<'primitive, 'arena>) -> Self {
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
            object_to_world: Arc::clone(&object_to_world),
            world_to_object,
            reverse_orientation,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
        }
    }
}
