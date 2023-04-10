//! Goniophotometric Light Source

use core::geometry::*;
use core::image_io::read_image;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::mipmap::MIPMap;
use core::mipmap::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a light source that uses goniophotometric diagrams encoded in 2D image maps to describe the emission
/// distribution of the light.
///
/// A goniophotometric diagram describes the angular distribution of luminance from a point light source; it is widely
/// used in illumination engineering to characterize lights.
///
/// The light coordinate system to always be at position (0, 0, 0) and pointing down the +z axis.
#[derive(Clone)]
pub struct GonioPhotometricLight {
    /// Light ID. This is usually the index of the light in the scene's light sources. Usefull for adding lights into
    /// `std::collections::HashMap`.
    pub id: usize,

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

    /// The MipMAP images of the goniophotometric diagram that scales the intensity based on the angular distribution of
    /// light.
    pub mipmap: Option<MIPMap<RGBSpectrum>>,
}

impl GonioPhotometricLight {
    /// Returns a new `GonioPhotometricLight`.
    ///
    /// * `id`               - Light ID.
    /// * `light_to_world`   - Transformation from light coordinate system to world coordinate system.
    /// * `medium_interface` - Participating medium.
    /// * `intensity`        - Intensity.
    /// * `texmap`           - Path to the image texture.
    pub fn new(
        id: usize,
        light_to_world: ArcTransform,
        medium_interface: MediumInterface,
        intensity: Spectrum,
        texmap: Option<&str>,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());
        let p_light = light_to_world.transform_point(&Point3f::ZERO);

        // Create `GonioPhotometricLight` MIP map.
        let mipmap = match texmap {
            Some(texmap) => match read_image(texmap) {
                Ok(image) => {
                    let mipmap = MIPMap::new(
                        &image.resolution,
                        &image.pixels,
                        FilteringMethod::Ewa,
                        ImageWrap::Repeat,
                        8.0,
                    );
                    Some(mipmap)
                }
                Err(e) => {
                    warn!("Error reading goniophotometric image texture '{}'. {}.", texmap, e);
                    None
                }
            },
            None => {
                warn!("No goniophotometric image texture provided.");
                None
            }
        };

        Self {
            id,
            light_type: LightType::DELTA_POSITION_LIGHT,
            medium_interface: medium_interface.clone(),
            light_to_world,
            world_to_light,
            p_light,
            intensity,
            mipmap,
        }
    }

    /// Returns how much to scale the amount of radiance.
    ///
    /// * `w` - Vector from hit point to light.
    fn scale(&self, w: &Vector3f) -> Spectrum {
        let mut wp = self.world_to_light.transform_vector(w).normalize();
        std::mem::swap(&mut wp.y, &mut wp.z);
        let theta = spherical_theta(&wp);
        let phi = spherical_phi(&wp);
        match self.mipmap.as_ref() {
            None => Spectrum::ONE,
            Some(mipmap) => {
                let st = Point2f::new(phi * INV_TWO_PI, theta * INV_PI);
                let rgb = mipmap.lookup_triangle(&st, 0.0).to_rgb();
                Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant))
            }
        }
    }
}

impl Light for GonioPhotometricLight {
    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Returns the light unique id. Usually the index in the scene's light sources.
    fn get_id(&self) -> usize {
        self.id
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, _u: &Point2f) -> Option<Li> {
        let wi = (self.p_light - hit.p).normalize();
        let pdf = 1.0;

        let p0 = hit.clone();
        let p1 = Hit::new_minimal(self.p_light, hit.time, Some(self.medium_interface.clone()));
        let vis = VisibilityTester::new(p0, p1);

        let value = self.intensity * self.scale(&-wi) / self.p_light.distance_squared(hit.p);
        Some(Li::new(wi, pdf, vis, value))
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        let spectrum = match self.mipmap.as_ref() {
            None => Spectrum::ONE,
            Some(mipmap) => {
                let st = Point2f::new(0.5, 0.5);
                let rgb = mipmap.lookup_triangle(&st, 0.5).to_rgb();
                Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant))
            }
        };
        FOUR_PI * self.intensity * spectrum
    }

    /// Returns the probability density with respect to solid angle for the light’s `sample_li()`.
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
        let ray = Ray::new(
            self.p_light,
            uniform_sample_sphere(&u1),
            INFINITY,
            time,
            self.medium_interface.inside.as_ref().map(Arc::clone),
        );
        let n_light = Normal3f::from(&ray.d);
        let scale = self.scale(&ray.d);

        Le::new(ray, n_light, 1.0, uniform_sphere_pdf(), self.intensity * scale)
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f) -> Pdf {
        Pdf::new(0.0, uniform_sphere_pdf())
    }

    /// Returns the number of samples to use for the light source.
    fn get_num_samples(&self) -> usize {
        1
    }

    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `it` - Point on a surface to evaluate emitted radiance.
    /// * `w`  - Outgoing direction.
    fn l(&self, _hit: &Hit, _w: &Vector3f) -> Spectrum {
        panic!("Invalid call to Light::l() for non area lights.")
    }
}

impl From<(&ParamSet, ArcTransform, Option<ArcMedium>, &str, usize)> for GonioPhotometricLight {
    /// Create a `GonioPhotometricLight` from given parameter set, light to world transform, medium, current working
    /// directory and id.
    ///
    /// * `p` - A tuple containing the parameter set, light to world transform, medium, current working directory and id.
    fn from(p: (&ParamSet, ArcTransform, Option<ArcMedium>, &str, usize)) -> Self {
        let (params, light_to_world, medium, cwd, id) = p;

        let intensity = params.find_one_spectrum("I", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let texmap = params.find_one_filename("mapname", Some(cwd));

        Self::new(
            id,
            Arc::clone(&light_to_world),
            MediumInterface::from(medium),
            intensity * sc,
            texmap.as_ref().map(|t| t.as_ref()),
        )
    }
}
