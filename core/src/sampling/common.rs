//! Common sampling functions.

use crate::geometry::*;
use crate::pbrt::*;
use crate::rng::*;

/// Generate 1D samples.
///
/// * `rng`       - Random number generator.
/// * `n_samples` - Number of samples.
/// * `jitter`    - Jitter the samples.
pub fn stratified_sample_1d(rng: &mut RNG, n_samples: usize, jitter: bool) -> Vec<Float> {
    let inv_n_samples = 1.0 / n_samples as Float;

    (0..n_samples)
        .map(|i| {
            let delta = if jitter { rng.uniform_float() } else { 0.5 };
            min((i as Float + delta) * inv_n_samples, ONE_MINUS_EPSILON)
        })
        .collect::<Vec<Float>>()
}

/// Generate 2D samples.
///
/// * `rng`    - Random number generator.
/// * `nx`     - Number of samples in x-direction.
/// * `ny`     - Number of samples in y-direction.
/// * `jitter` - Jitter the samples.
pub fn stratified_sample_2d(rng: &mut RNG, nx: usize, ny: usize, jitter: bool) -> Vec<Point2f> {
    let dx = 1.0 / nx as Float;
    let dy = 1.0 / ny as Float;

    (0..ny)
        .zip(0..nx)
        .map(|(y, x)| {
            let jx = if jitter { rng.uniform_float() } else { 0.5 };
            let jy = if jitter { rng.uniform_float() } else { 0.5 };
            Point2f::new(
                min((x as Float + jx) * dx, ONE_MINUS_EPSILON),
                min((y as Float + jy) * dy, ONE_MINUS_EPSILON),
            )
        })
        .collect::<Vec<Point2f>>()
}

/// Generate Latin Hypercube samples.
///
/// * `rng`       - Random number generator.
/// * `n_samples` - Number of samples.
/// * `n_dim`     - Number of dimensions.
pub fn latin_hypercube(rng: &mut RNG, n_samples: usize, n_dim: usize) -> Vec<Float> {
    let mut samples = vec![0.0; n_samples];
    let inv_n_samples = 1.0 / n_samples as Float;

    // Generate LHS samples along diagonal.
    for i in 0..n_samples {
        for j in 0..n_dim {
            let r = rng.uniform_float();
            let sj = (i as Float + r) * inv_n_samples;
            samples[n_dim * i + j] = min(sj, ONE_MINUS_EPSILON);
        }
    }

    // Permute LHS samples in each dimension.
    for i in 0..n_dim {
        for j in 0..n_samples {
            let other = j + rng.bounded_uniform_u32(0, (n_samples - j) as u32) as usize;
            samples.swap(n_dim * j + i, n_dim * other + i);
        }
    }

    samples
}

/// Sample points inside a unit circle by sampling a square of length 2 and
/// rejecting points outside the circle.
///
/// * `rng`       - Random number generator.
pub fn rejection_sample_disk(rng: &mut RNG) -> Point2f {
    loop {
        let rx = rng.uniform_float();
        let ry = rng.uniform_float();

        let x = 1.0 - 2.0 * rx;
        let y = 1.0 - 2.0 * ry;

        if x * x + y * y > 1.0 {
            continue;
        }

        return Point2f::new(x, y);
    }
}

/// Uniformly sample a direction on a hemisphere.
///
/// * `u` - The random sample point.
pub fn uniform_sample_hemisphere(u: &Point2f) -> Vector3f {
    let z = u[0];
    let r = max(0.0, 1.0 - z * z).sqrt();
    let phi = TWO_PI * u[1];
    Vector3f::new(r * cos(phi), r * sin(phi), z)
}

/// Returns the PDF for uniformly sampling a direction from a hemisphere.
#[inline]
pub fn uniform_hemisphere_pdf() -> Float {
    INV_TWO_PI
}

/// Uniformly sample a direction from a sphere.
///
/// * `u` - The random sample point.
pub fn uniform_sample_sphere(u: &Point2f) -> Vector3f {
    let z = 1.0 - 2.0 * u[0];
    let r = max(0.0, 1.0 - z * z).sqrt();
    let phi = TWO_PI * u[1];
    Vector3f::new(r * cos(phi), r * sin(phi), z)
}

/// Returns the PDF for uniformly sampling a direction from a sphere.
#[inline]
pub fn uniform_sphere_pdf() -> Float {
    INV_FOUR_PI
}

/// Uniformly sample a point on the disk.
///
/// * `u` - The random sample point.
pub fn uniform_sample_disk(u: &Point2f) -> Point2f {
    let r = u[0].sqrt();
    let theta = TWO_PI * u[1];
    Point2f::new(r * cos(theta), r * sin(theta))
}

/// Sample a point on a unit disk by mapping from a unit square to the unit
/// circle. The concentric mapping takes points in [-1, 1]^2 to unit disk by
/// uniformly mapping concentric squares to concentric circles.
///
/// * `u` - The random sample point.
pub fn concentric_sample_disk(u: &Point2f) -> Point2f {
    // Map uniform random numbers to [-1,1]^2.
    let u_offset = 2.0 * u - Vector2f::new(1.0, 1.0);

    // Handle degeneracy at the origin.
    if u_offset.x == 0.0 && u_offset.y == 0.0 {
        return Point2f::zero();
    }

    // Apply concentric mapping to point
    let (r, theta) = if abs(u_offset.x) > abs(u_offset.y) {
        (u_offset.x, PI_OVER_FOUR * (u_offset.y / u_offset.x))
    } else {
        (
            u_offset.y,
            PI_OVER_TWO - PI_OVER_FOUR * (u_offset.x / u_offset.y),
        )
    };

    r * Point2f::new(cos(theta), sin(theta))
}

/// Uniformly sample a direction from a cone of directions about the `(0, 0, 1)`
/// axis.
///
/// * `u`             - The random sample point.
/// * `cos_theta_max` - Cosine of the maximum angle of the beam.
pub fn uniform_sample_cone(u: &Point2f, cos_theta_max: Float) -> Vector3f {
    let cos_theta = (1.0 - u[0]) + u[0] * cos_theta_max;
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let phi = u[1] * TWO_PI;
    Vector3f::new(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta)
}

/// Uniformly sample a direction from a cone of directions about the z-axis in a
/// given coordinate system.
///
/// * `u`             - The random sample point.
/// * `x`             - The x-axis basis vector.
/// * `y`             - The x-axis basis vector.
/// * `z`             - The x-axis basis vector.
/// * `cos_theta_max` - Cosine of the maximum angle of the beam.
pub fn uniform_sample_cone_coordinate_system(
    u: &Point2f,
    cos_theta_max: Float,
    x: &Vector3f,
    y: &Vector3f,
    z: &Vector3f,
) -> Vector3f {
    let cos_theta = lerp(u[0], cos_theta_max, 1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let phi = u[1] * TWO_PI;
    cos(phi) * sin_theta * x + sin(phi) * sin_theta * y + cos_theta * z
}

/// Returns the PDF for sampling a direction from a cone of directions.
#[inline]
pub fn uniform_cone_pdf(cos_theta_max: Float) -> Float {
    1.0 / (TWO_PI * (1.0 - cos_theta_max))
}

/// Uniformly sample a point on an isosceles triangle.
///
/// * `u` - The random sample point.
pub fn uniform_sample_triangle(u: &Point2f) -> Point2f {
    let su0 = u[0].sqrt();
    Point2f::new(1.0 - su0, u[1] * su0)
}

/// Sample a direction on a hemisphere using cosine-weighted sampling.
///
/// * `u` - The random sample point.
#[inline]
pub fn cosine_sample_hemisphere(u: &Point2f) -> Vector3f {
    let d = concentric_sample_disk(u);
    let z = max(0.0, 1.0 - d.x * d.x - d.y * d.y).sqrt();
    Vector3f::new(d.x, d.y, z)
}

/// Returns the PDF for cosine-weighted sampling a direction from a hemisphere.
///
/// * `cos_theta` - Cosine term of incident radiance.
#[inline]
pub fn cosine_hemisphere_pdf(cos_theta: Float) -> Float {
    cos_theta * INV_PI
}

/// Weight samples using the balance heuristic.
///
/// * `nf`    - Number of samples taken from `f_pdf`.
/// * `f_pdf` - First sampling distribution.
/// * `ng`    - Number of samples taken from `g_pdf`.
/// * `g_pdf` - Second sampling distribution.
#[inline]
pub fn balance_heuristic(nf: Int, f_pdf: Float, ng: Int, g_pdf: Float) -> Float {
    (nf as Float * f_pdf) / (nf as Float * f_pdf + ng as Float * g_pdf)
}

/// Weight samples using the power heuristic.
///
/// * `nf`    - Number of samples taken from `f_pdf`.
/// * `f_pdf` - First sampling distribution.
/// * `ng`    - Number of samples taken from `g_pdf`.
/// * `g_pdf` - Second sampling distribution.
#[inline]
pub fn power_heuristic(nf: Int, f_pdf: Float, ng: Int, g_pdf: Float) -> Float {
    let f = nf as Float * f_pdf;
    let g = ng as Float * g_pdf;
    (f * f) / (f * f + g * g)
}
