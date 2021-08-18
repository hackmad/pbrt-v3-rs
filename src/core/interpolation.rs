//! Interpolation Functions

#![allow(dead_code)]
use crate::core::pbrt::*;

/// Interpolate the Catmull-Rom spline.
///
/// * `nodes`  - Interpolations nodes.
/// * `values` - Value of the function.
/// * `x`      - Variable to interpolate.
pub fn catmull_rom(nodes: &[Float], values: &[Float], x: Float) -> Float {
    let size = nodes.len();
    if !(x >= nodes[0] && x <= nodes[size - 1]) {
        return 0.0;
    }

    let idx = find_interval(size, |i| nodes[i] <= x);
    let x0 = nodes[idx];
    let x1 = nodes[idx + 1];
    let f0 = values[idx];
    let f1 = values[idx + 1];
    let width = x1 - x0;

    let d0 = if idx > 0 {
        width * (f1 - values[idx - 1]) / (x1 - nodes[idx - 1])
    } else {
        f1 - f0
    };

    let d1 = if idx + 2 < size {
        width * (values[idx + 2] - f0) / (nodes[idx + 2] - x0)
    } else {
        f1 - f0
    };

    let t = (x - x0) / (x1 - x0);
    let t2 = t * t;
    let t3 = t2 * t;
    (2.0 * t3 - 3.0 * t2 + 1.0) * f0
        + (-2.0 * t3 + 3.0 * t2) * f1
        + (t3 - 2.0 * t2 + t) * d0
        + (t3 - t2) * d1
}

/// Returns the weights and the index offset for Catmull-Rom spline.
///
/// * `nodes` - Interpolations nodes.
/// * `x`     - Variable to interpolate.
pub fn catmull_rom_weights(nodes: &[Float], x: Float) -> Option<([Float; 4], usize)> {
    // Return None if `x` is out of bounds.
    let size = nodes.len();
    if !(x >= nodes[0] && x <= nodes[size - 1]) {
        return None;
    }

    // Search for the interval `idx` containing `x`.
    let idx = find_interval(size, |i| nodes[i] <= x);
    let offset = idx - 1;
    let x0 = nodes[idx];
    let x1 = nodes[idx + 1];

    // Compute the `t` parameter and powers.
    let t = (x - x0) / (x1 - x0);
    let t2 = t * t;
    let t3 = t2 * t;

    // Compute initial node weights `w_1` and `w_2`.
    let mut weights = [0.0; 4];
    weights[1] = 2.0 * t3 - 3.0 * t2 + 1.0;
    weights[2] = -2.0 * t3 + 3.0 * t2;

    // Compute first node weight `w_0`.
    if idx > 0 {
        let w0 = (t3 - 2.0 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        let w0 = t3 - 2.0 * t2 + t;
        weights[0] = 0.0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    // Compute last node weight `w_3`.
    if idx + 2 < size {
        let w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        let w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0.0;
    }

    Some((weights, offset))
}

/// Returns the sampled value, fval and pdf.
///
/// * `x` - Locations x0, ..., xn where the function `f` is evaluated.
/// * `f` - Value of the function at each point `xi`.
/// * `F` - Discrete CDF values computed via `integrate_catmull_rom()`.
/// * `u` - Uniform random variate ξ.
#[allow(non_snake_case)]
pub fn sample_catmull_rom(
    x: &[Float],
    f: &[Float],
    F: &[Float],
    u: Float,
) -> (Float, Float, Float) {
    // Get number of samples.
    let n = x.len();

    // Map `u` to a spline interval by inverting `F`.
    let u = u * F[n - 1];
    let i = find_interval(n, |j| F[j] <= u);

    // Look up `x_i` and function values of spline segment `i`.
    let x0 = x[i];
    let x1 = x[i + 1];
    let f0 = f[i];
    let f1 = f[i + 1];
    let width = x1 - x0;

    // Approximate derivatives using finite differences.
    let d0 = if i > 0 {
        width * (f1 - f[i - 1]) / (x1 - x[i - 1])
    } else {
        f1 - f0
    };
    let d1 = if i + 2 < n {
        width * (f[i + 2] - f0) / (x[i + 2] - x0)
    } else {
        f1 - f0
    };

    // Re-scale `u` for continous spline sampling step.
    let u = (u - F[i]) / width;

    // Invert definite integral over spline segment and return solution.

    // Set initial guess for `t` by importance sampling a linear interpolant.
    let mut t = if f0 != f1 {
        (f0 - (max(0.0, f0 * f0 + 2.0 * u * (f1 - f0))).sqrt()) / (f0 - f1)
    } else {
        u / f0
    };

    let mut a: Float = 0.0;
    let mut b: Float = 1.0;
    let mut Fhat: Float;
    let mut fhat: Float;
    loop {
        // Fall back to a bisection step when `t` is out of bounds.
        if !(t > a && t < b) {
            t = 0.5 * (a + b);
        }

        // Evaluate target function and its derivative in Horner form.
        Fhat = t
            * (f0
                + t * (0.5 * d0
                    + t * ((1.0 / 3.0) * (-2.0 * d0 - d1) + f1 - f0
                        + t * (0.25 * (d0 + d1) + 0.5 * (f0 - f1)))));
        fhat = f0
            + t * (d0 + t * (-2.0 * d0 - d1 + 3.0 * (f1 - f0) + t * (d0 + d1 + 2.0 * (f0 - f1))));

        // Stop the iteration if converged.
        if abs(Fhat - u) < 1e-6 || b - a < 1e-6 {
            break;
        }

        // Update bisection bounds using updated `t`.
        if Fhat - u < 0.0 {
            a = t;
        } else {
            b = t;
        }

        // Perform a Newton step.
        t -= (Fhat - u) / fhat;
    }

    let fval = fhat;
    let pdf = fhat / F[n - 1];
    let sample = x0 + width * t;
    (sample, fval, pdf)
}

/// Returns sampled Catmull-Rom spline value, fval and pdf in 2D.
///
/// * `nodes1` - Locations x0, ..., xn where the function `f` is evaluated.
/// * `nodes2` - Locations y0, ..., yn where the function `f` is evaluated.
/// * `values` - Matrix of values in row-major order of the function at each
///              point (`xi`, `yi`).
/// * `cdf`    - Matrix in row-major order of discrete CDFs where each row is
///              computed via `integrate_catmull_rom()` on corresponding row
///              of `values`.
/// * `u`      - Uniform random variate ξ.
#[allow(non_snake_case)]
pub fn sample_catmull_rom_2d(
    nodes1: &[Float],
    nodes2: &[Float],
    values: &[Float],
    cdf: &[Float],
    alpha: Float,
    u: Float,
) -> (Float, Float, Float) {
    // Get number of nodes.
    let size2 = nodes2.len();

    // Determine offset and coefficients for the `alpha` parameter.
    let (weights, offset) = if let Some((w, o)) = catmull_rom_weights(nodes1, alpha) {
        (w, o)
    } else {
        return (0.0, 0.0, 0.0);
    };

    // Define a lambda function to interpolate table entries.
    let interpolate = |array: &[Float], idx: usize| -> Float {
        (0..4).fold(0.0, |a, i| {
            if weights[i] != 0.0 {
                a + array[(offset + i) * size2 + idx] * weights[i]
            } else {
                a
            }
        })
    };

    // Map `u` to a spline interval by inverting the interpolated `cdf`.
    let maximum = interpolate(cdf, size2 - 1);
    let u = u * maximum;
    let idx = find_interval(size2, |i| interpolate(cdf, i) <= u);

    // Look up node positions and interpolated function values.
    let f0 = interpolate(values, idx);
    let f1 = interpolate(values, idx + 1);
    let x0 = nodes2[idx];
    let x1 = nodes2[idx + 1];
    let width = x1 - x0;

    // Re-scale `u` using the interpolated `cdf`.
    let u = (u - interpolate(cdf, idx)) / width;

    // Approximate derivatives using finite differences of the interpolant.
    let d0 = if idx > 0 {
        width * (f1 - interpolate(values, idx - 1)) / (x1 - nodes2[idx - 1])
    } else {
        f1 - f0
    };
    let d1 = if idx + 2 < size2 {
        width * (interpolate(values, idx + 2) - f0) / (nodes2[idx + 2] - x0)
    } else {
        f1 - f0
    };

    // Invert definite integral over spline segment and return solution.

    // Set initial guess for `t` by importance sampling a linear interpolant.
    let mut t = if f0 != f1 {
        (f0 - (max(0.0, f0 * f0 + 2.0 * u * (f1 - f0))).sqrt()) / (f0 - f1)
    } else {
        u / f0
    };
    let mut a: Float = 0.0;
    let mut b: Float = 1.0;
    let mut Fhat: Float;
    let mut fhat: Float;
    loop {
        // Fall back to a bisection step when `t` is out of bounds.
        if !(t >= a && t <= b) {
            t = 0.5 * (a + b);
        }

        // Evaluate target function and its derivative in Horner form.
        Fhat = t
            * (f0
                + t * (0.5 * d0
                    + t * ((1.0 / 3.0) * (-2.0 * d0 - d1) + f1 - f0
                        + t * (0.25 * (d0 + d1) + 0.5 * (f0 - f1)))));
        fhat = f0
            + t * (d0 + t * (-2.0 * d0 - d1 + 3.0 * (f1 - f0) + t * (d0 + d1 + 2.0 * (f0 - f1))));

        // Stop the iteration if converged.
        if abs(Fhat - u) < 1e-6 || b - a < 1e-6 {
            break;
        }

        // Update bisection bounds using updated `t`.
        if Fhat - u < 0.0 {
            a = t;
        } else {
            b = t;
        }

        // Perform a Newton step.
        t -= (Fhat - u) / fhat;
    }

    // Return the sample position and function value.
    let fval = fhat;
    let pdf = fhat / maximum;
    let sample = x0 + width * t;
    (sample, fval, pdf)
}

/// Computes the integral and the auxilliary CDF for importance sampling.
///
/// * `x`      - Samples values.
/// * `values` - Value of the function.
pub fn integrate_catmull_rom(x: &[Float], values: &[Float]) -> (Vec<Float>, Float) {
    let n = x.len();
    let mut sum = 0.0;
    let mut cdf = vec![0.0; n];

    for i in 0..n - 1 {
        // Look up `x_i` and function values of spline segment `i`.
        let x0 = x[i];
        let x1 = x[i + 1];
        let f0 = values[i];
        let f1 = values[i + 1];
        let width = x1 - x0;

        // Approximate derivatives using finite differences.
        let d0 = if i > 0 {
            width * (f1 - values[i - 1]) / (x1 - x[i - 1])
        } else {
            f1 - f0
        };
        let d1 = if i + 2 < n {
            width * (values[i + 2] - f0) / (x[i + 2] - x0)
        } else {
            f1 - f0
        };

        // Keep a running sum and build a cumulative distribution function.
        sum += ((d0 - d1) * (1.0 / 12.0) + (f0 + f1) * 0.5) * width;
        cdf[i + 1] = sum;
    }

    (cdf, sum)
}

/// Inverts the Catmull-Rom spline function (as opposed to the definite
/// integral.
///
/// * `x`      - Samples values.
/// * `values` - Value of the function.
/// * `u`      - Uniform random variate ξ.
#[allow(non_snake_case)]
fn invert_catmull_rom(x: &[Float], values: &[Float], u: Float) -> Float {
    let n = x.len();

    // Stop when `u` is out of bounds.
    if u <= values[0] {
        return x[0];
    } else if u >= values[n - 1] {
        return x[n - 1];
    }

    // Map `u` to a spline interval by inverting `values`.
    let i = find_interval(n, |i| values[i] <= u);

    // Look up `x_i` and function values of spline segment `i`.
    let x0 = x[i];
    let x1 = x[i + 1];
    let f0 = values[i];
    let f1 = values[i + 1];
    let width = x1 - x0;

    // Approximate derivatives using finite differences
    let d0 = if i > 0 {
        width * (f1 - values[i - 1]) / (x1 - x[i - 1])
    } else {
        f1 - f0
    };
    let d1 = if i + 2 < n {
        width * (values[i + 2] - f0) / (x[i + 2] - x0)
    } else {
        f1 - f0
    };

    // Invert the spline interpolant using Newton-Bisection.
    let mut a: Float = 0.0;
    let mut b: Float = 1.0;
    let mut t: Float = 0.5;
    let mut Fhat: Float;
    let mut fhat: Float;
    loop {
        // Fall back to a bisection step when `t` is out of bounds.
        if !(t > a && t < b) {
            t = 0.5 * (a + b);
        }

        // Compute powers of `t`.
        let t2 = t * t;
        let t3 = t2 * t;

        // Set `Fhat` using Equation (8.27).
        Fhat = (2.0 * t3 - 3.0 * t2 + 1.0) * f0
            + (-2.0 * t3 + 3.0 * t2) * f1
            + (t3 - 2.0 * t2 + t) * d0
            + (t3 - t2) * d1;

        // Set `fhat` using Equation (not present).
        fhat = (6.0 * t2 - 6.0 * t) * f0
            + (-6.0 * t2 + 6.0 * t) * f1
            + (3.0 * t2 - 4.0 * t + 1.0) * d0
            + (3.0 * t2 - 2.0 * t) * d1;

        // Stop the iteration if converged.
        if abs(Fhat - u) < 1e-6 || b - a < 1e-6 {
            break;
        }

        // Update bisection bounds using updated `t`.
        if Fhat - u < 0.0 {
            a = t;
        } else {
            b = t;
        }

        // Perform a Newton step.
        t -= (Fhat - u) / fhat;
    }

    x0 + t * width
}

/// Compute the BSDF value using Fourier interpolation.
///
/// * `a`       - The weighted coefficients from a Fourier BSDF for order `m`.
/// * `cos_phi` - Cosine of the angle ΔΦ between pair of directions.
pub fn fourier(a: &[Float], cos_phi: f64) -> Float {
    let mut value = 0.0_f64;

    // Initialize cosine iterates.
    let mut cos_k_minus_one_phi = cos_phi;
    let mut cos_k_phi = 1.0_f64;
    for ak in a {
        // Add the current summand and update the cosine iterates.
        value += (*ak as f64) * cos_k_phi;
        let cos_k_plus_one_phi = 2.0_f64 * cos_phi * cos_k_phi - cos_k_minus_one_phi;
        cos_k_minus_one_phi = cos_k_phi;
        cos_k_phi = cos_k_plus_one_phi;
    }

    value as Float
}

/// Sample fourier
///
/// * `ak`    - The weighted coefficients from a Fourier BSDF for order `m`.
/// * `recip` - Contains 1 / i for i in [0..`m_max`].
/// * `u`     - Uniform random variate ξ.
#[allow(non_snake_case)]
pub fn sample_fourier(ak: &[Float], recip: &[Float], u: Float) -> (Float, Float, Float) {
    let m = ak.len();

    // Pick a side and declare bisection variables.
    let flip = u >= 0.5;
    let u = if flip { 1.0 - 2.0 * (u - 0.5) } else { u * 2.0 };

    let mut a: f64 = 0.0;
    let mut b = PI as f64;
    let mut phi = 0.5 * PI as f64;
    let mut F: f64;
    let mut f: f64;
    loop {
        // Evaluate `F(ϕ)` and its derivative `f(ϕ)`.

        // Initialize sine and cosine iterates.
        let cos_phi = phi.cos(); // Use f64 variant.
        let sin_phi = (max(0.0, 1.0 - cos_phi * cos_phi)).sqrt();
        let mut cos_phi_prev = cos_phi;
        let mut cos_phi_cur: f64 = 1.0;
        let mut sin_phi_prev = -sin_phi;
        let mut sin_phi_cur: f64 = 0.0;

        // Initialize `F` and `f` with the first series term.
        F = ak[0] as f64 * phi;
        f = ak[0] as f64;
        for k in 1..m {
            // Compute next sine and cosine iterates.
            let sin_phi_next = 2.0 * cos_phi * sin_phi_cur - sin_phi_prev;
            let cos_phi_next = 2.0 * cos_phi * cos_phi_cur - cos_phi_prev;
            sin_phi_prev = sin_phi_cur;
            sin_phi_cur = sin_phi_next;
            cos_phi_prev = cos_phi_cur;
            cos_phi_cur = cos_phi_next;

            // Add the next series term to `F` and `f`.
            F += ak[k] as f64 * recip[k] as f64 * sin_phi_next;
            f += ak[k] as f64 * cos_phi_next;
        }
        F -= (u * ak[0] * PI) as f64;

        // Update bisection bounds using updated ϕ.
        if F > 0.0 {
            b = phi;
        } else {
            a = phi;
        }

        // Stop the Fourier bisection iteration if converged.
        if abs(F) < 1e-6 || b - a < 1e-6 {
            break;
        }

        // Perform a Newton step given `f(ϕ)` and `F(ϕ)`.
        phi -= F / f;

        // Fall back to a bisection step when ϕ is out of bounds.
        if !(phi > a && phi < b) {
            phi = 0.5 * (a + b);
        }
    }

    // Potentially flip `ϕ` and return the result.
    if flip {
        phi = TWO_PI as f64 - phi;
    }
    let pdf = (INV_TWO_PI as f64 * f / ak[0] as f64) as Float;
    (f as Float, pdf, phi as Float)
}
