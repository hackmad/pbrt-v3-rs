//! Diffuse Area Light Source

use core::app::OPTIONS;
use core::geometry::*;
use core::light::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::rng::ONE_MINUS_EPSILON;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a basic area light source with uniform spatial and directional
/// radiance distribution.
#[derive(Clone)]
pub struct DiffuseAreaLight {
    /// Light source type.
    pub light_type: LightType,

    /// Used to trace multiple shadow rays to the light to compute soft shadows.
    pub n_samples: usize,

    /// Participating medium.
    pub medium_interface: MediumInterface,

    /// Transformation from light coordinate system to world coordinate system.
    pub light_to_world: ArcTransform,

    /// Transformation from world coordinate system to light coordinate system.
    pub world_to_light: ArcTransform,

    /// Emitted radiance.
    pub l_emit: Spectrum,

    /// Shape describing surface of the light source.
    pub shape: ArcShape,

    /// Surface area of the shape.
    pub area: Float,

    /// Indicates whether light source 2-sided.
    pub two_sided: bool,
}

impl DiffuseAreaLight {
    /// Returns a new `DiffuseAreaLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `l_emit`           - Emitted radiance.
    /// * `n_samples`        - Used to trace multiple shadow rays to the light
    ///                        to compute soft shadows. Default to 1.
    /// * `shape`            - Shape describing surface of the light source.
    /// * `two_sided`        - Indicates whether light source 2-sided.
    pub fn new(
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        l_emit: Spectrum,
        n_samples: usize,
        shape: ArcShape,
        two_sided: bool,
    ) -> Self {
        let world_to_light = Arc::clone(&light_to_world).inverse();
        let area = shape.area();
        Self {
            light_type: LightType(AREA_LIGHT),
            medium_interface: medium_interface.clone(),
            light_to_world: Arc::clone(&light_to_world),
            world_to_light: Arc::new(world_to_light),
            l_emit,
            n_samples,
            shape: Arc::clone(&shape),
            two_sided,
            area,
        }
    }

    /// Returns emitted radiance based on `two_sided` flag.
    ///
    /// * `intr` - The interaction point.
    /// * `w`    - Direction.
    fn l(&self, intr: &Hit, w: &Vector3f) -> Spectrum {
        if self.two_sided || intr.n.dot(w) > 0.0 {
            self.l_emit
        } else {
            Spectrum::new(0.0)
        }
    }
}

impl Light for DiffuseAreaLight {
    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, u: &Point2f) -> Li {
        let (mut p_shape_hit, pdf) = self.shape.sample_solid_angle(hit, u);
        p_shape_hit.medium_interface = Some(self.medium_interface.clone());

        let wi = p_shape_hit.p - hit.p;
        if pdf == 0.0 || wi.length_squared() == 0.0 {
            let pdf = 0.0;
            let visibility = None;
            let value = Spectrum::new(0.0);
            Li::new(wi, pdf, visibility, value)
        } else {
            let visibility = Some(VisibilityTester::new(hit.clone(), p_shape_hit.p));
            let value = self.l(&p_shape_hit, &(-wi));
            Li::new(wi, pdf, visibility, value)
        }
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        if self.two_sided {
            2.0 * self.l_emit * self.area * PI
        } else {
            self.l_emit * self.area * PI
        }
    }

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, hit: &Hit, wi: &Vector3f) -> Float {
        self.shape.pdf_solid_angle(hit, wi)
    }

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, u2: &Point2f, time: Float) -> Le {
        let (mut p_shape_hit, pdf_pos) = self.shape.sample_area(u1);
        p_shape_hit.medium_interface = Some(self.medium_interface.clone());
        let n_light = p_shape_hit.n;
        let pdf_dir: Float;
        let mut w: Vector3f;

        // Sample a cosine-weighted outgoing direction `w` for area light.
        if self.two_sided {
            let mut u = *u2;
            // Choose a side to sample and then remap u[0] to [0,1] before
            // applying cosine-weighted hemisphere sampling for the chosen side.
            if u[0] < 0.5 {
                u[0] = min(u[0] * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
            } else {
                u[0] = min((u[0] - 0.5) * 2.0, ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(&u);
                w.z *= -1.0;
            }
            pdf_dir = 0.5 * cosine_hemisphere_pdf(abs(w.z));
        } else {
            w = cosine_sample_hemisphere(&u2);
            pdf_dir = cosine_hemisphere_pdf(w.z);
        }

        let n = Vector3f::from(p_shape_hit.n);
        let (v1, v2) = coordinate_system(&n);
        w = w.x * v1 + w.y * v2 + w.z * n;

        let mut ray = p_shape_hit.spawn_ray(&w);
        ray.time = time;

        let value = self.l(&p_shape_hit, &w);
        Le::new(ray, n_light, pdf_pos, pdf_dir, value)
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        Pdf::new(0.0, uniform_sphere_pdf())
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>, ArcShape)> for DiffuseAreaLight {
    /// Create a `DiffuseAreaLight` from given parameter set, light to world transform
    /// medium, and shape.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform,
    ///         medium, and shape.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>, ArcShape)) -> Self {
        let (params, light_to_world, medium, shape) = p;

        let l = params.find_one_spectrum("L", Spectrum::new(1.0));
        let sc = params.find_one_spectrum("scale", Spectrum::new(1.0));
        let two_sided = params.find_one_bool("twosided", false);

        let mut n_samples = params.find_one_int("samples", params.find_one_int("nsamples", 1));
        if OPTIONS.quick_render {
            n_samples = max(1, n_samples / 4);
        }

        Self::new(
            light_to_world,
            MediumInterface::from(medium),
            l * sc,
            n_samples as usize,
            shape,
            two_sided,
        )
    }
}
