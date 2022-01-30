//! Projection Light Source

use core::geometry::*;
use core::image_io::read_image;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::mipmap::MIPMap;
use core::mipmap::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::cos_theta;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

// Near plane of projection.
const NEAR: Float = 1e-3;

// Far plane of projection.
const FAR: Float = 1e30;

/// Implements a light source that acts like a slide projector; it takes an image
/// map and projects its image out into the scene.
///
/// The light coordinate system to always be at position (0, 0, 0) and pointing
/// down the +z axis.
#[derive(Clone)]
pub struct ProjectionLight {
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

    /// The MipMAP images of the projection texture. This is a precise representation,
    /// compared to an image texture, of the projection function useful for being
    /// able to sample the projection distribution using Monte Carlo techniques.
    pub projection_map: Option<MIPMap<RGBSpectrum>>,

    /// Projective transformation.
    pub light_projection: Transform,

    /// Screen space extent of the projection.
    pub screen_bounds: Bounds2f,

    /// The cosine of the angle between the +z axis and the vector to a corner of
    /// the screen window. This value is used elsewhere to define the minimal
    /// cone of directions that encompasses the set of directions in which light
    /// is projected.
    pub cos_total_width: Float,
}

impl ProjectionLight {
    /// Returns a new `ProjectionLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `intensity`        - Intensity.
    /// * `texmap`           - Path to the image texture.
    /// * `fov`              - Field of view in degrees.
    pub fn new(
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        intensity: Spectrum,
        texmap: Option<&str>,
        fov: Float,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());
        let p_light = light_to_world.transform_point(&Point3f::ZERO);

        // Create `ProjectionLight` MIP map.
        let (projection_map, aspect) = match texmap {
            Some(texmap) => match read_image(texmap) {
                Ok(image) => {
                    let mipmap = MIPMap::new(
                        &image.resolution,
                        &image.pixels,
                        FilteringMethod::Ewa,
                        ImageWrap::Repeat,
                        8.0,
                    );
                    let aspect = image.resolution.x as Float / image.resolution.y as Float;
                    (Some(mipmap), aspect)
                }
                Err(e) => {
                    warn!(
                        "Error reading projection image texture '{}'. {}.",
                        texmap, e
                    );
                    (None, 1.0)
                }
            },
            None => {
                warn!("No projection image texture provided.");
                (None, 1.0)
            }
        };

        // Initialize `ProjectionLight` projection matrix.
        let screen_bounds = if aspect > 1.0 {
            Bounds2f::new(Point2f::new(-aspect, -1.0), Point2f::new(aspect, 1.0))
        } else {
            Bounds2f::new(
                Point2f::new(-1.0, -1.0 / aspect),
                Point2f::new(1.0, 1.0 / aspect),
            )
        };

        let light_projection = Transform::perspective(fov, NEAR, FAR);

        let screen_to_light = light_projection.inverse();
        let p_corner = Point3f::new(screen_bounds.p_max.x, screen_bounds.p_max.y, 0.0);
        let w_corner = Vector3f::from(screen_to_light.transform_point(&p_corner)).normalize();
        let cos_total_width = w_corner.z;

        Self {
            light_type: LightType::DELTA_POSITION_LIGHT,
            medium_interface: medium_interface.clone(),
            light_to_world,
            world_to_light,
            p_light,
            intensity,
            projection_map,
            light_projection,
            screen_bounds,
            cos_total_width,
        }
    }

    /// Returns how much light is projected in the given direction.
    ///
    /// * `w` - Vector from hit point to light.
    fn projection(&self, w: &Vector3f) -> Spectrum {
        let wl = self.world_to_light.transform_vector(w);

        // Discard directions behind projection light
        if wl.z < NEAR {
            return Spectrum::ZERO;
        }

        // Project point onto projection plane and compute light
        let p = self
            .light_projection
            .transform_point(&Point3f::new(wl.x, wl.y, wl.z));
        if !self.screen_bounds.contains(&Point2f::new(p.x, p.y)) {
            return Spectrum::ZERO;
        }

        match self.projection_map.as_ref() {
            None => Spectrum::ONE,
            Some(projection_map) => {
                let st = Point2f::from(self.screen_bounds.offset(&Point2f::new(p.x, p.y)));
                let rgb = projection_map.lookup_triangle(&st, 0.0).to_rgb();
                Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant))
            }
        }
    }
}

impl Light for ProjectionLight {
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

        let value = self.intensity * self.projection(&-wi) / self.p_light.distance_squared(hit.p);
        Li::new(wi, pdf, Some(vis), value)
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        let spectrum = match self.projection_map.as_ref() {
            None => Spectrum::ONE,
            Some(projection_map) => {
                let st = Point2f::new(0.5, 0.5);
                let rgb = projection_map.lookup_triangle(&st, 0.5).to_rgb();
                Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant))
            }
        };
        spectrum * self.intensity * TWO_PI * (1.0 - self.cos_total_width)
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
        let v = uniform_sample_cone(&u1, self.cos_total_width);
        let ray = Ray::new(
            self.p_light,
            self.light_to_world.transform_vector(&v),
            INFINITY,
            time,
            self.medium_interface.inside.as_ref().map(Arc::clone),
        );
        let n_light = Normal3f::from(&ray.d);
        let projection = self.projection(&ray.d);

        Le::new(
            ray,
            n_light,
            1.0,
            uniform_cone_pdf(self.cos_total_width),
            self.intensity * projection,
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

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>, &str)> for ProjectionLight {
    /// Create a `ProjectionLight` from given parameter set, light to world transform,
    /// medium and current working directory.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform,
    ///         medium and current working directory.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>, &str)) -> Self {
        let (params, light_to_world, medium, cwd) = p;

        let intensity = params.find_one_spectrum("I", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let fov = params.find_one_float("fov", 45.0);
        let texmap = params.find_one_filename("mapname", Some(cwd));

        Self::new(
            Arc::clone(&light_to_world),
            MediumInterface::from(medium),
            intensity * sc,
            texmap.as_ref().map(|t| t.as_ref()),
            fov,
        )
    }
}
