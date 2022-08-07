//! Infinate Area Light Source

use core::app::OPTIONS;
use core::geometry::*;
use core::image_io::*;
use core::interaction::*;
use core::light::*;
use core::medium::*;
use core::mipmap::*;
use core::paramset::*;
use core::pbrt::*;
use core::sampling::*;
use core::scene::*;
use core::spectrum::*;
use std::sync::{Arc, Mutex, RwLock};

/// Implements an infinite area light source using a latitude-longitude radiance
/// map.
#[derive(Clone)]
pub struct InfiniteAreaLight {
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

    /// The radiance map.
    pub l_map: MIPMap<RGBSpectrum>,

    /// World center.
    pub world_center: Arc<RwLock<Point3f>>,

    /// World radius.
    pub world_radius: Arc<RwLock<Float>>,

    /// 2-d distribution
    pub distribution: Distribution2D,
}

impl InfiniteAreaLight {
    /// Returns a new `InfiniteAreaLight`.
    ///
    /// * `light_to_world`   - Transformation from light coordinate system to
    ///                        world coordinate system.
    /// * `l`                - Power `L`.
    /// * `n_samples`        - Used to trace multiple shadow rays to the light
    ///                        to compute soft shadows. Default to 1.
    /// * `texmap`           - Path to the image to use for the radiance map.
    pub fn new(
        light_to_world: ArcTransform,
        l: Spectrum,
        n_samples: usize,
        texmap: Option<&str>,
    ) -> Self {
        let light_to_world = Arc::clone(&light_to_world);
        let world_to_light = Arc::new(light_to_world.inverse());

        let lrgb = l.to_rgb_spectrum();

        // Read texel data from texmap and initialize `l_map`.
        let (texels, resolution) = match texmap {
            None => (vec![lrgb], Point2::new(1_usize, 1_usize)),
            Some(texmap) => match read_image(texmap) {
                Ok(RGBImage { pixels, resolution }) => {
                    let texels = pixels.iter().map(|texel| *texel * lrgb).collect();
                    (texels, resolution)
                }
                Err(err) => {
                    warn!("Problem reading file '{}'. {}", texmap, err);
                    (vec![lrgb], Point2::new(1_usize, 1_usize))
                }
            },
        };

        let l_map = MIPMap::new(
            &resolution,
            &texels,
            FilteringMethod::Trilinear,
            ImageWrap::Repeat,
            0.0,
        );

        // Initialize sampling PDFs for infinite area light.

        // Compute scalar-valued image `img` from environment map.
        let img = compute_scalar_image(&l_map);

        // Compute sampling distributions for rows and columns of image
        let distribution = Distribution2D::new(img.clone());

        Self {
            light_type: LightType::INFINITE_LIGHT,
            medium_interface: MediumInterface::vacuum(),
            light_to_world,
            world_to_light,
            n_samples,
            l_map,
            distribution,
            world_center: Arc::new(RwLock::new(Point3f::ZERO)), // Calculated in preprocess().
            world_radius: Arc::new(RwLock::new(1.0)),           // Calculated in preprocess().
        }
    }
}

impl Light for InfiniteAreaLight {
    /// Initialize the light source before rendering begins.
    ///
    /// * `scene` - The scene.
    fn preprocess(&self, scene: &Scene) {
        let (world_center, world_radius) = scene.world_bound.bounding_sphere();
        *self.world_center.write().unwrap() = world_center;
        *self.world_radius.write().unwrap() = world_radius;
    }

    /// Returns the type of light.
    fn get_type(&self) -> LightType {
        self.light_type
    }

    /// Return the radiance arriving at an interaction point.
    ///
    /// * `hit` - The interaction hit point.
    /// * `u`   - Sample value for Monte Carlo integration.
    fn sample_li(&self, hit: &Hit, u: &Point2f) -> Li {
        // Find `(u,v)` sample coordinates in infinite light texture.
        let (uv, map_pdf) = self.distribution.sample_continuous(u);
        if map_pdf == 0.0 {
            Li::new(Vector3f::ZERO, 0.0, None, Spectrum::ZERO)
        } else {
            // Convert infinite light sample point to direction.
            let theta = uv[1] * PI;
            let phi = uv[0] * TWO_PI;
            let cos_theta = cos(theta);
            let sin_theta = sin(theta);
            let sin_phi = sin(phi);
            let cos_phi = cos(phi);

            let wi = self.light_to_world.transform_vector(&Vector3f::new(
                sin_theta * cos_phi,
                sin_theta * sin_phi,
                cos_theta,
            ));

            // Compute PDF for sampled infinite light direction.
            let mut pdf = map_pdf / (TWO_PI * PI * sin_theta);
            if sin_theta == 0.0 {
                pdf = 0.0;
            }

            // Return radiance value for infinite light direction.
            let world_radius = *self.world_radius.read().unwrap();
            let p0 = hit.clone();
            let p1 = Hit::new_minimal(
                hit.p + wi * (2.0 * world_radius),
                hit.time,
                Some(self.medium_interface.clone()),
            );
            let vis = VisibilityTester::new(p0, p1);

            let rgb = self.l_map.lookup_triangle(&uv, 0.0).to_rgb();
            let spectrum = Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant));
            Li::new(wi, pdf, Some(vis), spectrum)
        }
    }

    /// Return the total emitted power.
    fn power(&self) -> Spectrum {
        let world_radius = *self.world_radius.read().unwrap();
        let spectrum = Spectrum::from_rgb(
            &self
                .l_map
                .lookup_triangle(&Point2f::new(0.5, 0.5), 0.5)
                .to_rgb(),
            Some(SpectrumType::Illuminant),
        );
        PI * world_radius * world_radius * spectrum
    }

    /// Returns emitted radiance due to that light along a ray that escapes the
    /// scene bounds.
    ///
    /// * `ray` - The ray with differentials.
    fn le(&self, ray: &Ray) -> Spectrum {
        let w = self.world_to_light.transform_vector(&ray.d).normalize();
        let st = Point2f::new(spherical_phi(&w) * INV_TWO_PI, spherical_theta(&w) * INV_PI);
        Spectrum::from_rgb(
            &self.l_map.lookup_triangle(&st, 0.0).to_rgb(),
            Some(SpectrumType::Illuminant),
        )
    }

    /// Returns the probability density with respect to solid angle for the light’s
    /// `sample_li()`.
    ///
    /// * `hit` - The interaction hit point.
    /// * `wi`  - The incident direction.
    fn pdf_li(&self, _hit: &Hit, wi: &Vector3f) -> Float {
        let wi = self.world_to_light.transform_vector(wi);
        let theta = spherical_theta(&wi);
        let phi = spherical_phi(&wi);
        let sin_theta = sin(theta);
        if sin_theta == 0.0 {
            0.0
        } else {
            self.distribution
                .pdf(&Point2f::new(phi * INV_TWO_PI, theta * INV_PI))
                / (TWO_PI * PI * sin_theta)
        }
    }

    /// Returns a sampled light-carrying ray leaving the light source.
    ///
    /// * `u1`   - Sample values for Monte Carlo.
    /// * `u2`   - Sample values for Monte Carlo.
    /// * `time` - Time to use for the ray.
    fn sample_le(&self, u1: &Point2f, u2: &Point2f, time: Float) -> Le {
        // Compute direction for infinite light sample ray.
        let u = *u1;

        // Find `(u,v)` sample coordinates in infinite light texture.
        let (uv, map_pdf) = self.distribution.sample_continuous(&u);
        if map_pdf == 0.0 {
            Le::new(Ray::default(), Normal3f::ZERO, 0.0, 0.0, Spectrum::ZERO)
        } else {
            let world_center = *self.world_center.read().unwrap();
            let world_radius = *self.world_radius.read().unwrap();

            let theta = uv[1] * PI;
            let phi = uv[0] * TWO_PI;
            let cos_theta = cos(theta);
            let sin_theta = sin(theta);
            let sin_phi = sin(phi);
            let cos_phi = cos(phi);
            let d = -self.light_to_world.transform_vector(&Vector3f::new(
                sin_theta * cos_phi,
                sin_theta * sin_phi,
                cos_theta,
            ));
            let n_light = Normal3f::from(d);

            // Compute origin for infinite light sample ray.
            let (v1, v2) = coordinate_system(&(-d));
            let cd = concentric_sample_disk(u2);
            let p_disk = world_center + world_radius * (cd.x * v1 + cd.y * v2);
            let ray = Ray::new(p_disk + world_radius * -d, d, INFINITY, time, None);

            // Compute `InfiniteAreaLight` ray PDFs.
            let pdf_dir = if sin_theta == 0.0 {
                0.0
            } else {
                map_pdf / (TWO_PI * PI * sin_theta)
            };
            let pdf_pos = 1.0 / (PI * world_radius * world_radius);

            let rgb = self.l_map.lookup_triangle(&uv, 0.0).to_rgb();
            let spectrum = Spectrum::from_rgb(&rgb, Some(SpectrumType::Illuminant));
            Le::new(ray, n_light, pdf_pos, pdf_dir, spectrum)
        }
    }

    /// Returns the probability density for the light’s `sample_le()`.
    ///
    /// * `ray`     - The ray.
    /// * `n_light` - The normal.
    fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f) -> Pdf {
        let world_radius = *self.world_radius.read().unwrap();
        let d = -self.world_to_light.transform_vector(&ray.d);
        let theta = spherical_theta(&d);
        let phi = spherical_phi(&d);
        let uv = Point2f::new(phi * INV_TWO_PI, theta * INV_PI);
        let map_pdf = self.distribution.pdf(&uv);
        let pdf_dir = map_pdf / (TWO_PI * PI * sin(theta));
        let pdf_pos = 1.0 / (PI * world_radius * world_radius);
        Pdf::new(pdf_pos, pdf_dir)
    }

    /// Returns the number of samples to use for the light source.
    fn get_num_samples(&self) -> usize {
        self.n_samples
    }

    /// Returns the area light's emitted radiance in a given outgoing direction.
    ///
    /// * `it` - Point on a surface to evaluate emitted radiance.
    /// * `w`  - Outgoing direction.
    fn l(&self, _hit: &Hit, _w: &Vector3f) -> Spectrum {
        panic!("Invalid call to Light::l() for non area lights.")
    }
}

impl From<(&ParamSet, ArcTransform, &str)> for InfiniteAreaLight {
    /// Create a `InfiniteAreaLight` from given parameter set and, light to world
    /// transform and current working directory.
    ///
    /// * `p` - A tuple containing the parameter set and light to world transform
    ///         and current working directory.
    fn from(p: (&ParamSet, ArcTransform, &str)) -> Self {
        let (params, light_to_world, cwd) = p;

        let l = params.find_one_spectrum("L", Spectrum::ONE);
        let sc = params.find_one_spectrum("scale", Spectrum::ONE);
        let texmap = params.find_one_filename("mapname", Some(cwd));

        let mut n_samples = params.find_one_int("samples", params.find_one_int("nsamples", 1));
        if OPTIONS.quick_render {
            n_samples = max(1, n_samples / 4);
        }

        Self::new(
            light_to_world,
            l * sc,
            n_samples as usize,
            texmap.as_ref().map(|t| t.as_ref()),
        )
    }
}

/// Compute scalar-valued image `img` from environment map.
///
/// * `l_map`  - The environment map's mipmap.
fn compute_scalar_image(l_map: &MIPMap<Spectrum>) -> Vec<Vec<Float>> {
    let width = 2 * l_map.width();
    let height = 2 * l_map.height();
    let fwidth = 0.5 / min(width, height) as Float;

    let img: Arc<Mutex<Vec<Vec<Float>>>> = Arc::new(Mutex::new(vec![vec![]; height]));

    let n_threads = OPTIONS.threads();

    crossbeam::scope(|scope| {
        let (tx, rx) = crossbeam_channel::bounded(n_threads);

        for _ in 0..n_threads {
            let rxc = rx.clone();
            let img = Arc::clone(&img);
            scope.spawn(move |_| {
                for v in rxc.iter() {
                    let vp = (v as Float + 0.5) / height as Float;
                    let sin_theta = sin(PI * (v as Float + 0.5) / height as Float);

                    let mut img = img.lock().unwrap();
                    (*img)[v] = (0..width)
                        .map(|u| {
                            let up = (u as Float + 0.5) / width as Float;
                            l_map.lookup_triangle(&Point2f::new(up, vp), fwidth).y() * sin_theta
                        })
                        .collect();
                }
            });
        }
        drop(rx); // Drop extra rx since we've cloned one for each worker.

        // Send work.
        for v in 0..height {
            tx.send(v).unwrap();
        }
    })
    .unwrap();

    let mut v = img.lock().unwrap();
    std::mem::replace(&mut v, vec![])
}
