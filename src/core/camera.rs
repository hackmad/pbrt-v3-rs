//! Camera

#![allow(dead_code)]
use crate::core::film::*;
use crate::core::geometry::*;
use crate::core::light::*;
use crate::core::medium::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use std::fmt;
use std::sync::Arc;

/// Light trait provides common behavior.
pub trait Camera {
    /// Returns the common camera data.
    fn get_data(&self) -> &CameraData;

    /// Returns a ray corresponding to a given sample. It also returns, a floating
    /// point value that affects how much the radiance arriving at the film plane
    /// will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray(&self, sample: &CameraSample) -> (Ray, Float);

    /// Returns a main ray and rays shifted one pixel in x and y directions on
    /// the film plane for corresponding to a given sample. It also returns a
    /// floating point value that affects how much the radiance arriving at the
    /// film plane will contribute to final image.
    ///
    /// * `sample` - The sample.
    fn generate_ray_differential(&self, sample: &CameraSample) -> (Ray, Float) {
        let (mut ray, wt) = self.generate_ray(sample);
        let mut rd = RayDifferential::default();

        if wt == 0.0 {
            return (ray, 0.0);
        }

        // Find camera ray after shifting a fraction of a pixel in the x-direction.
        let mut wtx = 0.0;
        for eps in [0.05, -0.05].iter() {
            let mut sshift = *sample;
            sshift.p_film.x += eps;

            let (rx, w) = self.generate_ray(&sshift);
            rd.rx_origin = ray.o + (rx.o - ray.o) / *eps;
            rd.rx_direction = ray.d + (rx.d - ray.d) / *eps;

            wtx = w;
            if wtx != 0.0 {
                break;
            }
        }
        if wtx == 0.0 {
            return (ray, 0.0);
        }

        // Find camera ray after shifting a fraction of a pixel in the y-direction.
        let mut wty = 0.0;
        for eps in [0.05, -0.05].iter() {
            let mut sshift = *sample;
            sshift.p_film.y += eps;

            let (ry, w) = self.generate_ray(&sshift);
            rd.ry_origin = ray.o + (ry.o - ray.o) / *eps;
            rd.ry_direction = ray.d + (ry.d - ray.d) / *eps;

            wty = w;
            if wty != 0.0 {
                break;
            }
        }
        if wty == 0.0 {
            return (ray, 0.0);
        }

        ray.differentials = Some(rd);

        (ray, wt)
    }

    /// Evaluate the importance emitted from the point on the camera in a
    /// direction. The `include_raster_point` is true, then a raster position
    /// associated with the ray on the film is returned as well.
    ///
    /// * `ray`                  - The ray.
    /// * `include_raster_point` - Indicates whether or not to return the raster
    ///                            position.
    fn we(&self, _ray: &Ray, _include_raster: bool) -> (Spectrum, Option<Point2f>) {
        panic!("Camera::we() is not implemented");
    }

    /// Return the spatial and directional PDFs, as a tuple, for sampling a
    /// particular ray leaving the camera.
    ///
    /// * `ray` - The ray.
    fn pdf_we(&self, ray: &Ray) -> PDFResult;

    /// Returns a PDF value with respect to the solid angle at a reference point.
    ///
    /// * `interaction`          - The interaction point.
    /// * `u`                    - Used to sample point on the lens.
    /// * `include_raster_point` - Indicates whether or not to return the raster
    ///                            position.
    fn sample_wi(&self, _interaction: &dyn Interaction, _u: &Point2f) -> SampleResult {
        panic!("Camera::sample_wi() is not implemented");
    }
}

/// Atomic reference counted `Camera`.
pub type ArcCamera = Arc<dyn Camera + Send + Sync>;

/// Stores all of the sample values needed to specify a camera ray.
#[derive(Copy, Clone, Default)]
pub struct CameraSample {
    /// Point on the film to which the generated ray carries radiance.
    pub p_film: Point2f,

    /// The point on the lens the ray passes through (for cameras that
    /// support lenses).
    pub p_lens: Point2f,

    /// Time at which the ray should sample the scene. This should be linearly
    /// interpolated between shutter open and close time range.
    pub time: Float,
}

impl CameraSample {
    /// Create a new `CameraSample`.
    ///
    /// * `p_film` - Point on the film to which the generated ray carries radiance.
    /// * `p_lens` - The point on the lens the ray passes through (for cameras
    ///              that support lenses).
    /// * `time`   - Time at which the ray should sample the scene. This should
    ///              be linearly interpolated between shutter open and close
    ///              time range.
    pub fn new(p_film: Point2f, p_lens: Point2f, time: Float) -> Self {
        Self {
            p_film,
            p_lens,
            time,
        }
    }
}

impl fmt::Display for CameraSample {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "CameraSample<p_film: {:}, p_lens: {:}, time: {}>",
            self.p_film, self.p_lens, self.time
        )
    }
}

/// Stores the sampling result for a ray.
#[derive(Clone)]
pub struct SampleResult {
    /// The sample value.
    spectrum: Spectrum,

    /// Direction from lens to interaction point.
    wi: Vector3f,

    /// The PDF value.
    pdf: Float,

    /// Raster position.
    p_raster: Option<Point2f>,

    /// Visibility tester.
    vis: VisibilityTester,
}

impl SampleResult {
    /// Create a new `SampleResult`.
    ///
    /// * `spectrum` - The sample value.
    /// * `wi`       - Direction from lens to interaction point.
    /// * `pdf`      - The PDF value.
    /// * `p_raster` - Raster position.
    /// * `vis`      - Visibility tester.
    pub fn new(
        spectrum: Spectrum,
        wi: Vector3f,
        pdf: Float,
        p_raster: Option<Point2f>,
        vis: VisibilityTester,
    ) -> Self {
        Self {
            spectrum,
            wi,
            pdf,
            p_raster,
            vis,
        }
    }
}

/// Stores the spatial and directional PDFs for sampling a ray.
#[derive(Copy, Clone, Default)]
pub struct PDFResult {
    /// Spatial PDF.
    pos: Float,

    /// Directional PDF.
    dir: Float,
}

impl PDFResult {
    /// Create a new `PDFResult`.
    ///
    /// * `pos` - Spatial PDF.
    /// * `dir` - Directional PDF.
    pub fn new(pos: Float, dir: Float) -> Self {
        Self { pos, dir }
    }
}

/// Stores common camera parameters.
#[derive(Clone)]
pub struct CameraData {
    /// Animated transformation describing the camera's motion in the scene.
    pub camera_to_world: AnimatedTransform,

    /// Time when shutter is open.
    pub shutter_open: Float,

    /// Time when shutter is closed.
    pub shutter_close: Float,

    /// The film to capture the rendered image.
    pub film: Arc<Film>,

    /// Scattering medium the camera lies in.
    pub medium: Option<ArcMedium>,
}

impl CameraData {
    /// Creates a new instance of `CameraData`.
    ///
    /// * `camera_to_world` - Animated transformation describing the camera's
    ///                       motion in the scene.
    /// * `shutter_open`    - Time when shutter is open.
    /// * `shutter_close`   - Time when shutter is closed.
    /// * `film`            - The film to capture the rendered image.
    /// * `medium`          - Scattering medium the camera lies in.
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: Arc<Film>,
        medium: Option<ArcMedium>,
    ) -> Self {
        Self {
            camera_to_world,
            shutter_open,
            shutter_close,
            film,
            medium: medium.clone(),
        }
    }
}

/// Stores data for projective cameras.
#[derive(Clone)]
pub struct ProjectiveCameraData {
    /// Transformation from camera space to screen space.
    pub camera_to_screen: Transform,

    /// Transformation from raster space to camera space.
    pub raster_to_camera: Transform,

    /// Transformation from raster space to screen space.
    pub raster_to_screen: Transform,

    /// Radius of camera lens.
    pub lens_radius: Float,

    /// Focal distance.
    pub focal_distance: Float,
}

impl ProjectiveCameraData {
    /// Create a new instance of `ProjectiveCamera`.
    ///
    /// * `camera_data`      - Common camera parameters.
    /// * `camera_to_screen` - Transformation from camera space to screen space.
    /// * `screen_window`    - Bounds of screen space.
    /// * `lens_radius`      - Radius of camera lens.
    /// * `focal_distance`   - Focal distance.
    pub fn new(
        camera_data: &CameraData,
        camera_to_screen: Transform,
        screen_window: Bounds2f,
        lens_radius: Float,
        focal_distance: Float,
    ) -> Self {
        // Compute projective camera transformations.
        // Compute projective camera screen transformations.
        let screen_to_raster = Transform::scale(
            camera_data.film.full_resolution.x as Float,
            camera_data.film.full_resolution.y as Float,
            1.0,
        ) * Transform::scale(
            1.0 / (screen_window.p_max.x - screen_window.p_min.x),
            1.0 / (screen_window.p_min.y - screen_window.p_max.y),
            1.0,
        ) * Transform::translate(&Vector3f::new(
            -screen_window.p_min.x,
            -screen_window.p_max.y,
            0.0,
        ));

        let raster_to_screen = screen_to_raster.inverse();
        let raster_to_camera = camera_to_screen.inverse() * raster_to_screen;

        Self {
            camera_to_screen,
            raster_to_camera,
            raster_to_screen,
            lens_radius,
            focal_distance,
        }
    }
}
