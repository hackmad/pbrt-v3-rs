//! Realistic Camera

use core::app::OPTIONS;
use core::camera::*;
use core::efloat::*;
use core::film::*;
use core::float_file::parse_float_file;
use core::geometry::*;
use core::low_discrepency::*;
use core::medium::*;
use core::paramset::*;
use core::pbrt::*;
use core::reflection::*;
use core::report_stats;
use core::stats::*;
use core::{stat_inc, stat_percent, stat_register_fns};
use std::mem::swap;
use std::sync::{Arc, Mutex};
use std::thread;

stat_percent!(
    "Camera/Rays vignetted by lens system",
    VIGNETTED_RAYS,
    TOTAL_RAYS,
    realistic_camera_stats_vignetted_rays,
);

stat_register_fns!(realistic_camera_stats_vignetted_rays);

/// Number of samples for exit pupil bounds.
const N_SAMPLES: usize = 64_usize;

/// Realistic camera implements a camera consisting of multiple lens elements.
pub struct RealisticCamera {
    /// Common camera parameters.
    pub data: CameraData,

    /// If true, a modified version of weighting camera ray that only accounts for cos^4(θ) term is used.
    pub simple_weighting: bool,

    /// The lens element interfaces from left to right.
    pub element_interfaces: Vec<LensElementInterface>,

    /// Bounds of the exit pupil (the set of points on the rear element that do carry light through the lens system).
    pub exit_pupil_bounds: Vec<Bounds2f>,
}

impl RealisticCamera {
    /// Create a new realistic camera.
    ///
    /// * `camera_to_world`   - Animated transformation describing the camera's motion in the scene.
    /// * `aperture_diameter` - Diameter of the aperture stop in millimiters.
    /// * `focus_distance`    - Distance to the desired plane of focus.
    /// * `simple_weighting`  - If true, a modified version of weighting camera ray that only accounts for cos^4(θ) term
    ///                         is used.
    /// * `lens_data`         - A vector containing lens data as a flattened list of lens interfaces from left to right,
    ///                         each containing 4 properties: `curvature radius`, `thickness`, `index of refraction`,
    ///                         `aperture diameter`. These are related to the `LensElementInterface` properties (note
    ///                         that the data contains the aperture diameter as opposed to aperture radius).
    /// * `shutter_open`      - Time when shutter is open.
    /// * `shutter_close`     - Time when shutter is closed.
    /// * `lens_radius`       - Radius of camera lens.
    /// * `focal_distance`    - Focal distance.
    /// * `film`              - The film to capture the rendered image.
    /// * `medium`            - Scattering medium the camera lies in.
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        aperture_diameter: Float,
        focus_distance: Float,
        simple_weighting: bool,
        lens_data: &[Float],
        film: Film,
        medium: Option<ArcMedium>,
    ) -> Self {
        register_stats();

        let film_clone = film;
        let film_diagonal = film_clone.diagonal;

        let data = CameraData::new(camera_to_world, shutter_open, shutter_close, film_clone, medium);

        let n = lens_data.len();
        let mut element_interfaces = Vec::<LensElementInterface>::with_capacity(n);

        for i in (0..n).step_by(4) {
            let mut diameter: Float = lens_data[i + 3];

            if lens_data[i] == 0.0 {
                if aperture_diameter > lens_data[i + 3] {
                    warn!(
                        "Specified aperture diameter {} is greater than maximum \
                        possible {}. Clamping it.",
                        aperture_diameter,
                        lens_data[i + 3]
                    );
                } else {
                    diameter = aperture_diameter;
                }
            }

            element_interfaces.push(LensElementInterface::new(
                lens_data[i] * 0.001,
                lens_data[i + 1] * 0.001,
                lens_data[i + 2],
                diameter * 0.001 / 2.0,
            ));
        }

        // Compute exit pupil bounds at sampled points on the film.
        let mut camera = Self {
            data,
            simple_weighting,
            element_interfaces,
            exit_pupil_bounds: Vec::<Bounds2f>::with_capacity(N_SAMPLES),
        };

        // Compute lens-film distance for given focus distance.
        let fb = camera.focus_binary_search(focus_distance);
        info!("Binary search focus: {} -> {}", fb, camera.focus_distance(fb));

        let thickness = camera.focus_thick_lens(focus_distance);
        camera.element_interfaces.last_mut().unwrap().thickness = thickness;

        info!(
            "Thick lens focus: {} -> {}",
            thickness,
            camera.focus_distance(thickness)
        );

        // Compute exit pupil bounds at sampled points on the film.
        camera.exit_pupil_bounds = compute_exit_pupil_bounds(&camera, film_diagonal);

        if simple_weighting {
            error!(
                "'simple_weighting' option with RealisticCamera no longer \
                necessarily matches regular camera images. Further, pixel \
                values will vary a bit depending on the aperture size. See \
                this discussion for details: \
                https://github.com/mmp/pbrt-v3/issues/162#issuecomment-348625837"
            );
        }

        camera
    }

    /// Returns the z-depth value of the rear element.
    fn lens_rear_z(&self) -> Float {
        self.element_interfaces
            .last()
            .expect("element_interfaces is empty")
            .thickness
    }

    /// Returns the z-depth value of the front element.
    fn lens_front_z(&self) -> Float {
        self.element_interfaces.iter().map(|element| element.thickness).sum()
    }

    /// Returns the aperture radius of the rear element in meters.
    fn rear_element_radius(&self) -> Float {
        self.element_interfaces
            .last()
            .expect("element_interfaces is empty")
            .aperture_radius
    }

    /// Computes intersections of a ray fom the film with each element in turn, terminating it, and returning `None` if
    /// its path is blocked along the way through the lens system. Otherwise it returns a `Ray` initialized with the
    /// exiting ray in camera space.
    ///
    /// * `r_camera` - The camera ray to trace.
    fn trace_lenses_from_film(&self, r_camera: &Ray) -> Option<Ray> {
        // During traversal, element_z tracks the intercept of the current lens element. Because the ray is starting
        // from the film, the lenses are traversed in reverse order compared to how they are stored in
        // `element_interfaces`.
        let mut element_z = 0.0;

        // Transform `r_camera` from camera to lens system space.
        let camera_to_lens = Transform::scale(1.0, 1.0, -1.0);
        let mut r_lens = camera_to_lens.transform_ray(r_camera);

        for (i, element) in self.element_interfaces.iter().enumerate().rev() {
            // Update ray from film accounting for interaction with `element`.
            element_z -= element.thickness;

            // Compute intersection of ray with lens element.
            let t: Float;
            let mut n = Normal3f::ZERO;
            let is_stop = element.curvature_radius == 0.0;
            if is_stop {
                // The reflected ray computed in the previous lens element interface may be pointed towards film
                // plane(+z) in some extreme situations; in such cases; `t` becomes negative.
                if r_lens.d.z >= 0.0 {
                    return None;
                } else {
                    t = (element_z - r_lens.o.z) / r_lens.d.z;
                }
            } else {
                let z_center = element_z + element.curvature_radius;
                if let Some((t_hit, n_hit)) = intersect_spherical_element(element.curvature_radius, z_center, &r_lens) {
                    t = t_hit;
                    n = n_hit;
                } else {
                    return None;
                }
            }
            assert!(t >= 0.0);

            // Test intersection point against element aperture.
            let p_hit = r_lens.at(t);
            let r2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
            if r2 > element.aperture_radius * element.aperture_radius {
                return None;
            }
            r_lens.o = p_hit;

            // Update ray path for element interface interaction.
            if !is_stop {
                let eta_i = element.eta;

                let eta_t = if i > 0 && self.element_interfaces[i - 1].eta != 0.0 {
                    self.element_interfaces[i - 1].eta
                } else {
                    1.0
                };
                if let Some(w) = refract(&(-r_lens.d).normalize(), &n, eta_i / eta_t) {
                    r_lens.d = w;
                } else {
                    return None;
                }
            }
        }

        // Transform `r_len` from lens system space back to camera space.
        let lens_to_camera = Transform::scale(1.0, 1.0, -1.0);
        Some(lens_to_camera.transform_ray(&r_lens))
    }

    /// Computes intersections of a ray fom the scene with each element in turn, terminating it, and returning `None` if
    /// its path is blocked along the way through the lens system. Otherwise it returns a `Ray` initialized with the
    /// exiting ray in camera space.
    ///
    /// * `r_camera` - The camera ray to trace.
    fn trace_lenses_from_scene(&self, r_camera: &Ray) -> Option<Ray> {
        // During traversal, element_z tracks the intercept of the current lens element. Because the ray is starting
        // from the scene, the lenses are traversed in forward order as they are stored in `element_interfaces`.
        let mut element_z = -self.lens_front_z();

        // Transform `r_camera` from camera to lens system space.
        let camera_to_lens = Transform::scale(1.0, 1.0, -1.0);
        let mut r_lens = camera_to_lens.transform_ray(r_camera);

        for (i, element) in self.element_interfaces.iter().enumerate() {
            // Compute intersection of ray with lens element.
            let t: Float;
            let mut n = Normal3f::ZERO;
            let is_stop = element.curvature_radius == 0.0;
            if is_stop {
                t = (element_z - r_lens.o.z) / r_lens.d.z;
            } else {
                let z_center = element_z + element.curvature_radius;
                if let Some((t_hit, n_hit)) = intersect_spherical_element(element.curvature_radius, z_center, &r_lens) {
                    t = t_hit;
                    n = n_hit;
                } else {
                    return None;
                }
            }
            assert!(t >= 0.0);

            // Test intersection point against element aperture.
            let p_hit = r_lens.at(t);
            let r2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
            if r2 > element.aperture_radius * element.aperture_radius {
                return None;
            }
            r_lens.o = p_hit;

            // Update ray path for element interface interaction.
            if !is_stop {
                let eta_i = if i == 0 || self.element_interfaces[i - 1].eta == 0.0 {
                    1.0
                } else {
                    self.element_interfaces[i - 1].eta
                };

                let eta_t = if self.element_interfaces[i].eta != 0.0 {
                    self.element_interfaces[i].eta
                } else {
                    1.0
                };

                if let Some(wt) = refract(&(-r_lens.d).normalize(), &n, eta_i / eta_t) {
                    r_lens.d = wt;
                } else {
                    return None;
                }
            }
            element_z += element.thickness;
        }

        // Transform `r_len` from lens system space back to camera space.
        let lens_to_camera = Transform::scale(1.0, 1.0, -1.0);
        Some(lens_to_camera.transform_ray(&r_lens))
    }

    /// Computes both pairs of cardinal points ([pz0, pz1], [fz0, fz1]) for the lens system where `pz0`, `pz1` are
    /// z-depths of the focal point and `fz0`, `fz1` are the z-depths of the principal plane.
    fn compute_thick_lens_approximation(&self) -> ([Float; 2], [Float; 2]) {
        // Find height `x` from optical axis for parallel rays.
        //
        // Use a small fraction of the film's diagonal extent so that the rays we use don't get blocked by the aperture
        // stop. This works well unless aperture stop is extremely small.
        let x = 0.001 * self.data.film.diagonal;

        // Compute cardinal points for film side of lens system.
        let r_scene = Ray::new(
            Point3f::new(x, 0.0, self.lens_front_z() + 1.0),
            Vector3f::new(0.0, 0.0, -1.0),
            INFINITY,
            0.0,
            None,
        );

        let (pz0, fz0) = if let Some(r_film) = self.trace_lenses_from_scene(&r_scene) {
            compute_cardinal_points(&r_scene, &r_film)
        } else {
            panic!(
                "Unable to trace ray from scene to film for thick lens \
                approximation. Is aperture stop extremely small?"
            );
        };

        // Compute cardinal points for scene side of lens system.
        let r_film = Ray::new(
            Point3f::new(x, 0.0, self.lens_rear_z() - 1.0),
            Vector3f::new(0.0, 0.0, 1.0),
            INFINITY,
            0.0,
            None,
        );
        let (pz1, fz1) = if let Some(r_scene) = self.trace_lenses_from_film(&r_film) {
            compute_cardinal_points(&r_film, &r_scene)
        } else {
            panic!(
                "Unable to trace ray from film to scene for thick lens \
                approximation. Is aperture stop extremely small?"
            );
        };

        ([pz0, pz1], [fz0, fz1])
    }

    /// Focuses the lens system at a given depth and returns the offset along the z-axis from the film where the lens
    /// system should be placed.
    ///
    /// * `focus_distance` - Focus distance.
    fn focus_thick_lens(&self, focus_distance: Float) -> Float {
        // Get the cardinal points.
        let (pz, fz) = self.compute_thick_lens_approximation();
        info!(
            "Cardinal points: p' = {} f' = {}, p = {}, f = {}",
            pz[0], fz[0], pz[1], fz[1]
        );
        info!("Effective focal length {}", fz[0] - pz[0]);

        // Compute translation of lens `delta` to focus at `focus_distance`.
        let f = fz[0] - pz[0];
        let z = -focus_distance;
        let c = (pz[1] - z - pz[0]) * (pz[1] - z - 4.0 * f - pz[0]);
        assert!(
            c > 0.0,
            "Coefficient must be positive. It looks focus_distance: {} \
            is too short for a given lenses configuration",
            focus_distance
        );

        let delta = 0.5 * (pz[1] - z + pz[0] - c.sqrt());
        self.lens_rear_z() + delta
    }

    fn focus_binary_search(&self, focus_distance: Float) -> Float {
        // Find `film_distance_lower`, `film_distance_upper` that bound focus distance.
        let mut film_distance_lower = self.focus_thick_lens(focus_distance);
        let mut film_distance_upper = film_distance_lower;

        while self.focus_distance(film_distance_lower) > focus_distance {
            film_distance_lower *= 1.005;
        }

        while self.focus_distance(film_distance_upper) < focus_distance {
            film_distance_upper /= 1.005;
        }

        // Do binary search on film distances to focus.
        for _i in 0..20 {
            let fmid = 0.5 * (film_distance_lower + film_distance_upper);
            let mid_focus = self.focus_distance(fmid);
            if mid_focus < focus_distance {
                film_distance_lower = fmid;
            } else {
                film_distance_upper = fmid;
            }
        }
        0.5 * (film_distance_lower + film_distance_upper)
    }

    fn focus_distance(&self, film_distance: Float) -> Float {
        // Find offset ray from film center through lens.
        let bounds = self.bound_exit_pupil(0.0, 0.001 * self.data.film.diagonal);

        const SCALE_FACTORS: [Float; 3] = [0.1, 0.01, 0.001];
        let mut lu = 0.0;
        let mut ray: Option<Ray> = None;

        // Try some different and decreasing scaling factor to find focus ray more quickly when `aperturediameter` is
        // too small. (e.g. 2 [mm] for `aperturediameter` with wide.22mm.dat),
        for scale in SCALE_FACTORS {
            lu = scale * bounds.p_max[Axis::X];

            ray = self.trace_lenses_from_film(&Ray::new(
                Point3f::new(0.0, 0.0, self.lens_rear_z() - film_distance),
                Vector3f::new(lu, 0.0, film_distance),
                INFINITY,
                0.0,
                None,
            ));

            if ray.is_some() {
                break;
            }
        }

        // Compute distance `zFocus` where ray intersects the principal axis.
        if let Some(r) = ray {
            let t_focus = -r.o.x / r.d.x;
            let z_focus = r.at(t_focus).z;
            if z_focus < 0.0 {
                INFINITY
            } else {
                z_focus
            }
        } else {
            error!(
                "Focus ray at lens pos({}, 0) didn't make it through the lenses \
                with film distance {}!",
                lu, film_distance
            );
            INFINITY
        }
    }

    /// Compute a 2-d bounding box of the exit pupil as seen from a point along a segment on the film plane by tracing
    /// rays through the lens system at a set of points on a plane tangent to the rear lens element.
    ///
    /// Only the x-axis coordinates are used for the segment because the lens system is radially symmetric around the
    /// optical axis `z`. So the exit pupil bounds will be radially symmetric as well.
    ///
    /// * `p_film_x0` - First point along x-axis of film plane.
    /// * `p_film_x1` - Second point along x-axis of film plane.
    fn bound_exit_pupil(&self, p_film_x0: Float, p_film_x1: Float) -> Bounds2f {
        let mut pupil_bounds = Bounds2f::EMPTY;
        const N_SAMPLES: usize = 1024 * 1024;
        const N_SAMPLES_SQRT: usize = 1024;

        let mut n_exiting_rays = 0;

        // Compute bounding box of the projection of rear element on sampling plane.
        let rear_radius = self.rear_element_radius();
        let proj_rear_bounds = Bounds2f::new(
            Point2f::new(-1.5 * rear_radius, -1.5 * rear_radius),
            Point2f::new(1.5 * rear_radius, 1.5 * rear_radius),
        );
        for i in 0..N_SAMPLES {
            // Find location of sample points on `x` segment and rear lens element.
            let p_film = Point3f::new(
                lerp((i as Float + 0.5) / N_SAMPLES as Float, p_film_x0, p_film_x1),
                0.0,
                0.0,
            );

            let u = [radical_inverse(0, i as u64), radical_inverse(1, i as u64)];
            let p_rear = Point3f::new(
                lerp(u[0], proj_rear_bounds.p_min.x, proj_rear_bounds.p_max.x),
                lerp(u[1], proj_rear_bounds.p_min.y, proj_rear_bounds.p_max.y),
                self.lens_rear_z(),
            );

            // Expand pupil bounds if ray makes it through the lens system
            if pupil_bounds.contains(&Point2f::new(p_rear.x, p_rear.y))
                || self
                    .trace_lenses_from_film(&Ray::new(p_film, p_rear - p_film, INFINITY, 0.0, None))
                    .is_some()
            {
                pupil_bounds = pupil_bounds.union(&Point2f::new(p_rear.x, p_rear.y));
                n_exiting_rays += 1;
            }
        }

        // Return the entire element bounds if no rays made it through the lens system.
        if n_exiting_rays == 0 {
            // *NOTE*: This should happen a few times when called via focus_binary_search() -> focus_distance() ->
            // bound_exit_pupil(). Not when computing the actual exit pupil bounds.
            warn!(
                "Unable to find exit pupil in x = [{}, {}] on film.",
                p_film_x0, p_film_x1
            );
            proj_rear_bounds
        } else {
            pupil_bounds.expand(2.0 * proj_rear_bounds.diagonal().length() / N_SAMPLES_SQRT as Float)
        }
    }

    /// Returns the bounds and area on the exit pupil for a given point on the film plane.
    ///
    /// * `p_film`      - Point on the film plane.
    /// * `lens_sample` - Point on the lens the ray passes through for given sample.
    fn sample_exit_pupil(&self, p_film: &Point2f, lens_sample: &Point2f) -> (Point3f, Float) {
        // Find exit pupil bound for sample distance from film center.
        let r_film = (p_film.x * p_film.x + p_film.y * p_film.y).sqrt();
        let mut r_index = (r_film / (self.data.film.diagonal / 2.0) * self.exit_pupil_bounds.len() as Float) as usize;
        r_index = min(self.exit_pupil_bounds.len() - 1, r_index);
        let pupil_bounds = self.exit_pupil_bounds[r_index];
        let sample_bounds_area = pupil_bounds.area();

        // Generate sample point inside exit pupil bound.
        let p_lens = pupil_bounds.lerp(lens_sample);

        // Return sample point rotated by angle of `p_film` with `+x` axis.
        let sin_theta = if r_film != 0.0 { p_film.y / r_film } else { 0.0 };

        let cos_theta = if r_film != 0.0 { p_film.x / r_film } else { 1.0 };

        let p = Point3f::new(
            cos_theta * p_lens.x - sin_theta * p_lens.y,
            sin_theta * p_lens.x + cos_theta * p_lens.y,
            self.lens_rear_z(),
        );

        (p, sample_bounds_area)
    }
}

impl From<(&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>, &str)> for RealisticCamera {
    /// Create a `RealisticCamera` from given parameter set, animated transform, film, medium and current working
    /// directory.
    ///
    /// * `p` - A tuple containing  parameter set, animated transform, film, medium and current working directory.
    fn from(p: (&ParamSet, &AnimatedTransform, Film, Option<ArcMedium>, &str)) -> Self {
        let (params, cam2world, film, medium, cwd) = p;

        // Extract common camera parameters from `ParamSet`
        let mut shutter_open = params.find_one_float("shutteropen", 0.0);
        let mut shutter_close = params.find_one_float("shutterclose", 1.0);
        if shutter_close < shutter_open {
            warn!(
                "Shutter close time [{}] < shutter open [{}].  Swapping them.",
                shutter_close, shutter_open
            );
            swap(&mut shutter_close, &mut shutter_open);
        }

        // Realistic camera-specific parameters
        let lens_file = match params.find_one_filename("lensfile", Some(cwd)) {
            Ok(Some(f)) => f,
            Ok(None) => panic!("No lens description file supplied!"),
            Err(e) => panic!("{e}"),
        };
        let aperture_diameter = params.find_one_float("aperturediameter", 1.0);
        let focus_distance = params.find_one_float("focusdistance", 10.0);
        let simple_weighting = params.find_one_bool("simpleweighting", true);

        // Load element data from lens description file
        let lens_data = match parse_float_file(&lens_file) {
            Ok(data) => data,
            Err(err) => panic!("Error reading lens specification file '{}'. {}.", lens_file, err),
        };

        if lens_data.len() % 4 != 0 {
            panic!(
                "Excess values in lens specification file '{}'; must be \
                multiple-of-four values, read {}.",
                lens_file,
                lens_data.len()
            );
        }

        Self::new(
            cam2world.clone(),
            shutter_open,
            shutter_close,
            aperture_diameter,
            focus_distance,
            simple_weighting,
            &lens_data,
            film,
            medium,
        )
    }
}

impl Camera for RealisticCamera {
    /// Returns the camera data.
    fn get_data(&self) -> &CameraData {
        &self.data
    }

    /// Returns a ray corresponding to a given sample. It also returns, a floating point value that affects how much the
    /// radiance arriving at the film plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float) {
        stat_inc!(TOTAL_RAYS, 1);

        // Find point on film, `p_film`, corresponding to `sample.p_film`.
        let s = Point2f::new(
            sample.p_film.x / self.data.film.full_resolution.x as Float,
            sample.p_film.y / self.data.film.full_resolution.y as Float,
        );
        let p_film2 = self.data.film.get_physical_extent().lerp(&s);
        let p_film = Point3f::new(-p_film2.x, p_film2.y, 0.0);

        // Trace ray from `p_film` through lens system.
        let (p_rear, exit_pupil_bounds_area) =
            self.sample_exit_pupil(&Point2f::new(p_film.x, p_film.y), &sample.p_lens);

        let r_film = Ray::new(
            p_film,
            p_rear - p_film,
            INFINITY,
            lerp(sample.time, self.data.shutter_open, self.data.shutter_close),
            self.data.medium.as_ref().map(Arc::clone),
        );

        if let Some(ray) = self.trace_lenses_from_film(&r_film) {
            // Finish initialization of `RealisticCamera` ray.
            let mut ray = self.data.camera_to_world.transform_ray(&ray);
            ray.d = ray.d.normalize();

            // Return weighting for `RealisticCamera` ray.
            let cos_theta = r_film.d.normalize().z;
            let cos_4_theta = (cos_theta * cos_theta) * (cos_theta * cos_theta);
            let weight = if self.simple_weighting {
                cos_4_theta * exit_pupil_bounds_area / self.exit_pupil_bounds[0].area()
            } else {
                (self.data.shutter_close - self.data.shutter_open) * (cos_4_theta * exit_pupil_bounds_area)
                    / (self.lens_rear_z() * self.lens_rear_z())
            };
            (ray, weight)
        } else {
            stat_inc!(VIGNETTED_RAYS, 1);
            (Ray::default(), 0.0)
        }
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, _ray: &Ray) -> PDFResult {
        panic!("NOT IMPLEMENTED");
    }
}

/// Stores information about a single lens element interface. A lens interface intersects the optical axis at a position
/// `z`.
#[derive(Copy, Clone, Default)]
pub struct LensElementInterface {
    /// Radius of lens curvature:
    /// * Positive: convex curvature
    /// * Negative: concave curvature.
    /// * Zero:     aperture stop.
    pub curvature_radius: Float,

    /// Lens thickness is the distance to the next interface to the right or the distance to the film plane for the
    /// rearmost interface. This is in meters.
    pub thickness: Float,

    /// Index of refraction.
    pub eta: Float,

    /// Aperture radius describes the extent above and below the optical axis. This is in meters.
    pub aperture_radius: Float,
}

impl LensElementInterface {
    /// Create a new lens element interface.
    ///
    /// * `curvature_radius`  - Radius of lens curvature:
    ///                         * Positive: convex curvature
    ///                         * Negative: concave curvature.
    ///                         * Zero:     aperture stop.
    /// * `thicknesst`        - Lens thickness is the distance to the next interface to the right or the distance to the
    ///                         film plane for the rearmost interface.
    /// * `eta`               - Index of refraction.
    /// * `aperture_radius`   - Aperture radius describes the extent above and below the optical axis.
    fn new(curvature_radius: Float, thickness: Float, eta: Float, aperture_radius: Float) -> Self {
        Self {
            curvature_radius,
            thickness,
            eta,
            aperture_radius,
        }
    }

    /// Create a new lens element interface from 4 floating point values.
    ///
    /// * `aperture_diameter` - Diameter of the aperture stop in millimeters.
    /// * `lens_data`         - An array [curvature_radius, thickness, index of refraction, aperture_diameter]. The
    ///                         curvature_radius, thickness and aperture_diameter are in millimiters.
    #[allow(unused)]
    fn from_lens_data(aperture_diameter: Float, lens_data: &[Float]) -> Self {
        let ad = if lens_data[0] == 0.0 {
            if aperture_diameter > lens_data[3] {
                warn!(
                    "WARN: Specified aperture_diameter {} > max possible {}. \
                    Clamping it to max.",
                    aperture_diameter, lens_data[3]
                );
                lens_data[3]
            } else {
                aperture_diameter
            }
        } else {
            lens_data[3]
        };

        Self {
            curvature_radius: lens_data[0] * 0.001,
            thickness: lens_data[1] * 0.001,
            eta: lens_data[2],
            aperture_radius: ad * 0.001 / 2.0,
        }
    }
}

/// Calculate the parametric `t` value along a ray where it intersects a spherical element's interface.
///
/// * `radius`   - Radius of curvature.
/// * `z_center` - Z-axis intercept.
/// * `ray`      - Camera ray.
fn intersect_spherical_element(radius: Float, z_center: Float, ray: &Ray) -> Option<(Float, Normal3f)> {
    // Compute `t0` and `t1` for ray-element intersection.
    let o = ray.o - Vector3f::new(0.0, 0.0, z_center);
    let a = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
    let b = 2.0 * (ray.d.x * o.x + ray.d.y * o.y + ray.d.z * o.z);
    let c = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;

    if let Some((t0, t1)) = Quadratic::solve_float(a, b, c) {
        // Select intersection `t` based on ray direction and element curvature.
        let use_closer_t = (ray.d.z > 0.0) ^ (radius < 0.0);
        let t = if use_closer_t { min(t0, t1) } else { max(t0, t1) };

        if t < 0.0 {
            None
        } else {
            // Compute surface normal of element at ray intersection point
            let n = Normal3f::from(Vector3f::from(o + t * ray.d))
                .normalize()
                .face_forward(&-ray.d);
            Some((t, n))
        }
    } else {
        None
    }
}

/// Computes the z-depths of the focal point `pz` and of the principal plane `fz` for the given rays. Note that it
/// assumes that the rays are in camera space but returns values along the optical axis in lens space.
///
/// * `r_in`  - Ray in.
/// * `r_out` - Ray out.
fn compute_cardinal_points(r_in: &Ray, r_out: &Ray) -> (Float, Float) {
    let tf = -r_out.o.x / r_out.d.x;
    let fz = -r_out.at(tf).z;
    let tp = (r_in.o.x - r_out.o.x) / r_out.d.x;
    let pz = -r_out.at(tp).z;

    (pz, fz)
}

/// Compute the exit pupil bounds.
///
/// * `camera`        - The camera.
/// * `film_diagonal` - The diaogonal of the film's physical area in meters.
fn compute_exit_pupil_bounds(camera: &RealisticCamera, film_diagonal: Float) -> Vec<Bounds2f> {
    let exit_pupil_bounds = Arc::new(Mutex::new(vec![Bounds2f::EMPTY; N_SAMPLES]));

    let progress = create_progress_bar(N_SAMPLES as u64);
    progress.set_message("Calculate exit pupil");

    let fac = 1.0 / N_SAMPLES as Float * film_diagonal / 2.0;

    let n_threads = OPTIONS.threads();

    thread::scope(|scope| {
        let (tx, rx) = crossbeam_channel::bounded(n_threads);

        for _ in 0..n_threads {
            let rxc = rx.clone();
            let exit_pupil_bounds = Arc::clone(&exit_pupil_bounds);
            let progress = &progress;
            scope.spawn(move || {
                for i in rxc.iter() {
                    let r0 = i as Float * fac;
                    let r1 = (i + 1) as Float * fac;
                    let mut ep = exit_pupil_bounds.lock().unwrap();
                    (*ep)[i] = camera.bound_exit_pupil(r0, r1);
                    progress.inc(1);
                }

                // Report per thread statistics.
                report_stats!();
            });
        }
        drop(rx); // Drop extra rx since we've cloned one for each worker.

        // Send work.
        for i in 0..N_SAMPLES {
            tx.send(i).unwrap();
        }
    });

    progress.finish_with_message("Exit pupil calculated");

    let mut ep = exit_pupil_bounds.lock().unwrap();
    std::mem::replace(&mut ep, vec![])
}
