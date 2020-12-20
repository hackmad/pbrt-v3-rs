//! Sampling functions

#![allow(dead_code)]
use super::geometry::*;
use super::pbrt::*;
use super::rng::*;

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

/// Generate 1D samples.
///
/// * `rng`       - Random number generator.
/// * `n_samples` - Number of samples.
/// * `jitter`    - Jitter the samples.
pub fn stratified_sample_1d(rng: &mut RNG, n_samples: usize, jitter: bool) -> Vec<Float> {
    let inv_n_samples = 1.0 / n_samples as Float;

    (0..n_samples)
        .map(|i| {
            let delta: Float = if jitter { rng.uniform() } else { 0.5 };
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
            let jx = if jitter { rng.uniform() } else { 0.5 };
            let jy = if jitter { rng.uniform() } else { 0.5 };
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
            let r: Float = rng.uniform();
            let sj = (i as Float + r) * inv_n_samples;
            samples[n_dim * i + j] = min(sj, ONE_MINUS_EPSILON);
        }
    }

    // Permute LHS samples in each dimension.
    for i in 0..n_dim {
        for j in 0..n_samples {
            let other = j + rng.bounded_uniform(0, n_samples - j);
            let t = samples[n_dim * j + i];
            samples[n_dim * j + i] = samples[n_dim * other + i];
            samples[n_dim * other + i] = t;
        }
    }

    samples
}
