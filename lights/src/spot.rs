//! Spot Light Source

use core::geometry::*;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::cos_theta;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements an spot light source that emits light in a cone of directions from
/// its position.
///
/// The spotlight light coordinate system to always be at position (0, 0, 0) and
/// pointing down the +z axis.
#[derive(Clone)]
pub struct SpotLight {
    /// Light source type.
    pub light_type: LightType,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// Transformation from world coordinate system to light coordinate system.
    pub world_to_light: ArcTransform,

    /// Position.
    pub p_light: Point3f,

    /// Intensity.
    pub intensity: Spectrum,

    /// Cosine of overall angular width of the cone.
    cos_total_width: Float,

    /// Cosine of angle at which falloff starts.
    cos_falloff_start: Float,
}

impl SpotLight {
    /// Returns a new `SpotLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `intensity`        - Intensity.
    /// * `total_width`      - Overall angular width of the cone in degrees.
    /// * `falloff_start`    - Angle at which fallof starts in degrees.
    pub fn new(
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        intensity: Spectrum,
        total_width: Float,
        falloff_start: Float,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());
        let p_light = light_to_world.transform_point(&Point3f::default());

        Self {
            light_type: LightType::DELTA_POSITION_LIGHT,
            medium_interface: medium_interface.clone(),
            light_to_world,
            world_to_light,
            p_light,
            intensity,
            cos_total_width: cos(total_width.to_radians()),
            cos_falloff_start: cos(falloff_start.to_radians()),
        }
    }

    /// Returns the distribution of light accounting for the spotlight cone.
    ///
    /// * `w` - Vector from hit point to light.
    fn falloff(&self, w: &Vector3f) -> Float {
        let wl = self.world_to_light.transform_vector(w).normalize();
        let cos_theta = wl.z;
        if cos_theta < self.cos_total_width {
            return 0.0;
        }
        if cos_theta >= self.cos_falloff_start {
            return 1.0;
        }

        // Compute falloff inside spotlight cone.
        let delta =
            (cos_theta - self.cos_total_width) / (self.cos_falloff_start - self.cos_total_width);
        (delta * delta) * (delta * delta)
    }
}

impl Light for SpotLight {
    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, _u: &Point2f) -> Li {
        let wi = (self.p_light - hit.p).normalize();
        let pdf = 1.0;

        let p0 = hit.clone();
        let p1 = Hit::new_minimal(self.p_light, hit.time, hit.medium_interface.clone());
        let vis = VisibilityTester::new(p0, p1);

        let value = self.intensity * self.falloff(&-wi) / self.p_light.distance_squared(hit.p);
        Li::new(wi, pdf, Some(vis), value)
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        self.intensity * TWO_PI * (1.0 - 0.5 * (self.cos_falloff_start + self.cos_total_width))
    }

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, _hit: &Hit, _wi: &Vector3f) -> Float {
        0.0
    }

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, _u2: &Point2f, time: Float) -> Le {
        let w = uniform_sample_cone(&u1, self.cos_total_width);

        let ray = Ray::new(
            self.p_light,
            self.light_to_world.transform_vector(&w),
            INFINITY,
            time,
            self.medium_interface.inside.clone(),
        );

        Le::new(
            ray,
            Normal3f::from(w),
            1.0,
            uniform_cone_pdf(self.cos_total_width),
            self.intensity * self.falloff(&w),
        )
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f) -> Pdf {
        let pdf_dir =
            if cos_theta(&self.world_to_light.transform_vector(&ray.d)) >= self.cos_total_width {
                uniform_cone_pdf(self.cos_total_width)
            } else {
                0.0
            };
        Pdf::new(0.0, pdf_dir)
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>)> for SpotLight {
    /// Create a `SpotLight` from given parameter set, light to world transform
    /// and medium.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform
    ///         and medium.
    #[rustfmt::skip]
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>)) -> Self {
        let (params, light_to_world, medium) = p;

        let intensity = params.find_one_spectrum("I", Spectrum::new(1.0));
        let sc = params.find_one_spectrum("scale", Spectrum::new(1.0));
        let cone_angle = params.find_one_float("coneangle", 30.0);
        let cone_delta = params.find_one_float("conedeltaangle", 5.0);

        // Compute spotlight world to light transformation.
        let from = params.find_one_point3f("from", Point3f::new(0.0, 0.0, 0.0));
        let to = params.find_one_point3f("to", Point3f::new(0.0, 0.0, 1.0));
        let dir = (to - from).normalize();
        let (du, dv) = coordinate_system(&dir);

        let dir_to_z = Transform::from(
            Matrix4x4::new(
                du.x,  du.y,  du.z,  0.0, 
                dv.x,  dv.y,  dv.z,  0.0,
                dir.x, dir.y, dir.z, 0.0,
                0.0,   0.0,   0.0,   1.0,
            )
        );

        let t = Vector3f::new(from.x, from.y, from.z);
        let l2w = (*light_to_world).clone() * &Transform::translate(&t) * &dir_to_z.inverse();
        Self::new(
            Arc::new(l2w), 
            MediumInterface::from(medium),
            intensity * sc,
            cone_angle,
            cone_angle - cone_delta,
        )
    }
}
