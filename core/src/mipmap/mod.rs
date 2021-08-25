//! MIPMap

#![allow(dead_code)]

use crate::geometry::*;
use crate::memory::*;
use crate::pbrt::*;
use crate::texture::*;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign};
use std::sync::Arc;

mod cache;
mod convert_in;
mod tex_info;

// Re-export
pub use cache::*;
pub use convert_in::*;
pub use tex_info::*;

/// Size of the weights lookup table.
const WEIGHT_LUT_SIZE: usize = 128;

/// Enumeration for the image wrapping convention for out-of-bounds texels.
#[derive(Copy, Clone, Hash, PartialEq)]
pub enum ImageWrap {
    /// Repeat.
    Repeat,
    /// Black.
    Black,
    /// Clamp.
    Clamp,
}

/// Holds details for image reconstruction.
#[derive(Copy, Clone, Default)]
pub struct ResampleWeight {
    /// Offset to the first texel.
    pub first_texel: usize,

    /// The weight values for 4 texels.
    pub weight: [Float; 4],
}

/// MIPMap texture filtering methods.
#[derive(Copy, Clone, PartialEq, Hash)]
pub enum FilteringMethod {
    /// Trilinear interpolation.
    Trilinear,

    /// Elliptically weighted average.
    Ewa,
}

/// Implements methods for efficient texture filtering with spatially varying
/// filter widths.
#[derive(Clone)]
pub struct MIPMap<T> {
    /// MIP-Map method to use.
    filtering_method: FilteringMethod,

    /// Determines how to handle out-of-bounds texels.
    wrap_mode: ImageWrap,

    /// Image resolution.
    resolution: Point2<usize>,

    /// Stores the image pyramid of increasingly lower resolution prefiltered
    /// versions of the original image.
    pyramid: Vec<BlockedArray<T>>,

    /// Precomputed lookup table of Gaussian filter function values.
    weight_lut: [Float; WEIGHT_LUT_SIZE],

    /// Used to clamp the ellipse eccentricity (EWA).
    /// Set to 0 if EWA is not being used.
    max_anisotropy: Float,
}

/// Atomic reference counted `MIPMap`.
pub type ArcMIPMap<T> = Arc<MIPMap<T>>;

impl<T> MIPMap<T>
where
    T: Copy
        + Clone
        + Default
        + Mul<Float, Output = T>
        + MulAssign<Float>
        + Div<Float, Output = T>
        + DivAssign<Float>
        + Add<T, Output = T>
        + AddAssign
        + Clamp<Float>,
{
    /// * `resolution`       - Image resolution.
    /// * `img`              - Image data.
    /// * `filtering_method` - MIPMap filtering method to use.
    /// * `wrap_mode`        - Determines how to handle out-of-bounds texels.
    /// * `max_anisotropy`   - Used to clamp the ellipse eccentricity (EWA).
    ///                        Set to 0 if EWA is not being used.
    pub fn new(
        resolution: &Point2<usize>,
        img: &[T],
        filtering_method: FilteringMethod,
        wrap_mode: ImageWrap,
        max_anisotropy: Float,
    ) -> Self {
        let mut resampled_image: Vec<T> = vec![];

        let resolution = if !resolution[0].is_power_of_two() || !resolution[1].is_power_of_two() {
            // Resample image to power-of-two resolution.
            let res_pow2 = Point2::new(
                resolution[0].next_power_of_two(),
                resolution[1].next_power_of_two(),
            );
            info!(
                "Resampling MIPMap from {}x{} to {}x{}",
                resolution[0], resolution[1], res_pow2[0], res_pow2[1],
            );

            // Resample image in `s` direction.
            let s_weights = resample_weights(resolution[0], res_pow2[0]);
            resampled_image.resize(res_pow2[0] * res_pow2[1], T::default());

            // Apply `s_weights` in the `s` direction.
            for t in 0..resolution[1] {
                for s in 0..res_pow2[0] {
                    // Compute texel `(s, t)` in `s`-zoomed image.
                    resampled_image[t * res_pow2[0] + s] = T::default(); // Should be zero.
                    for j in 0..4 {
                        let orig_s = s_weights[s].first_texel + j;
                        let orig_s = match wrap_mode {
                            ImageWrap::Repeat => rem(orig_s, resolution[0]),
                            ImageWrap::Clamp => clamp(orig_s, 0, resolution[0] - 1),
                            _ => orig_s,
                        };

                        if orig_s < resolution[0] {
                            resampled_image[t * res_pow2[0] + s] +=
                                img[t * resolution[0] + orig_s] * s_weights[s].weight[j];
                        }
                    }
                }
            }

            // Resample image in `t` direction.
            let t_weights = resample_weights(resolution[1], res_pow2[1]);

            // Apply `t_weights` in the `t` direction.
            let mut work_data = vec![T::default(); res_pow2[1]];
            for s in 0..res_pow2[0] {
                for t in 0..res_pow2[1] {
                    work_data[t] = T::default(); // Should be zero.
                    for j in 0..4 {
                        let offset = t_weights[t].first_texel + j;
                        let offset = match wrap_mode {
                            ImageWrap::Repeat => rem(offset, resolution[1]),
                            ImageWrap::Clamp => clamp(offset, 0, resolution[1] - 1),
                            _ => offset,
                        };

                        if offset < resolution[1] {
                            work_data[t] +=
                                resampled_image[offset * res_pow2[0] + s] * t_weights[t].weight[j];
                        }
                    }
                }
                for t in 0..res_pow2[1] {
                    resampled_image[t * res_pow2[0] + s] = work_data[t].clamp_default();
                }
            }

            res_pow2
        } else {
            *resolution
        };

        // Initialize levels of MIPMap from image.
        let n_levels = 1 + Log2::log2(max(resolution[0], resolution[1])) as usize;
        let mut pyramid: Vec<BlockedArray<T>> = Vec::with_capacity(n_levels);

        // Initialize most detailed level of MIPMap
        pyramid.push(BlockedArray::from_slice(
            resolution[0],
            resolution[1],
            if resampled_image.len() > 0 {
                &resampled_image
            } else {
                img
            },
        ));

        for i in 1..n_levels {
            // Initialize i^th MIPMap level from `i-1` level.
            let s_res = max(1, pyramid[i - 1].u_size() / 2);
            let t_res = max(1, pyramid[i - 1].v_size() / 2);
            pyramid.push(BlockedArray::new(s_res, t_res));

            // Filter four texels from finer level of pyramid.
            for t in 0..t_res {
                for s in 0..s_res {
                    let tx0 = texel(&pyramid, wrap_mode, i - 1, 2 * s, 2 * t);
                    let tx1 = texel(&pyramid, wrap_mode, i - 1, 2 * s + 1, 2 * t);
                    let tx2 = texel(&pyramid, wrap_mode, i - 1, 2 * s, 2 * t + 1);
                    let tx3 = texel(&pyramid, wrap_mode, i - 1, 2 * s + 1, 2 * t + 1);
                    pyramid[i][(s, t)] = (tx0 + tx1 + tx2 + tx3) * 0.25;
                }
            }
        }

        // Initialize EWA filter weights.
        let mut weight_lut = [0.0; WEIGHT_LUT_SIZE];
        let alpha = 2.0;
        for i in 0..WEIGHT_LUT_SIZE {
            let r2 = i as Float / (WEIGHT_LUT_SIZE - 1) as Float;
            weight_lut[i] = (-alpha * r2).exp() - (-alpha).exp();
        }

        Self {
            filtering_method,
            wrap_mode,
            resolution,
            pyramid,
            weight_lut,
            max_anisotropy,
        }
    }

    /// Returns the width of the highest resolution level.
    pub fn width(&self) -> usize {
        self.resolution[0]
    }

    /// Returns the height of the highest resolution level.
    pub fn height(&self) -> usize {
        self.resolution[1]
    }

    /// Returns the number of MIPMap levels.
    pub fn levels(&self) -> usize {
        self.pyramid.len()
    }

    /// Applies the appropriate filter method based on `method` over the texture
    /// samples to remove high frequencies.
    ///
    /// * `st`   - The sample point coordinates (s, t).
    /// * `dst0` - Length of first elliptical axis.
    /// * `dst1` - Length of second elliptical axis.
    pub fn lookup(&self, st: &Point2f, dst0: &Vector2f, dst1: &Vector2f) -> T {
        match self.filtering_method {
            FilteringMethod::Trilinear => {
                let width = max(
                    max(abs(dst0[0]), abs(dst0[1])),
                    max(abs(dst1[0]), abs(dst1[1])),
                );
                self.lookup_triangle(st, width)
            }
            FilteringMethod::Ewa => self.lookup_ewa(st, &dst0, &dst1),
        }
    }

    /// Uses a triangle filter over the texture samples to remove high
    /// frequencies.
    ///
    /// * `st`    - The sample point coordinates (s, t).
    /// * `width` - Filter width (default to 0).
    pub fn lookup_triangle(&self, st: &Point2f, width: Float) -> T {
        // Compute MIPMap level for trilinear filtering.
        let levels = self.levels();
        let level = (levels - 1) as Float + max(width, 1e-8).log2();

        // Perform trilinear interpolation at appropriate MIPMap level.
        if level < 0.0 {
            self.triangle(0, st)
        } else if level >= (levels - 1) as Float {
            texel(&self.pyramid, self.wrap_mode, levels - 1, 0, 0)
        } else {
            // Do lerp() manually to avoid adding trait bound on T such that
            // `Float: Mul<T, Output=T>` and messing up Float multiplications.
            let i_level = level.floor() as usize;
            let delta = level - i_level as Float;
            self.triangle(i_level, st) * (1.0 - delta) + self.triangle(i_level + 1, st) * delta
        }
    }

    /// Uses the EWA filter over the texture samples to remove high frequencies.
    ///
    /// * `st`   - The sample point coordinates (s, t).
    /// * `dst0` - Length of first elliptical axis.
    /// * `dst1` - Length of second elliptical axis.
    fn lookup_ewa(&self, st: &Point2f, dst0: &Vector2f, dst1: &Vector2f) -> T {
        // Compute ellipse minor and major axes.
        let (dst0, mut dst1) = if dst0.length_squared() < dst1.length_squared() {
            (*dst1, *dst0)
        } else {
            (*dst0, *dst1)
        };

        let major_length = dst0.length();
        let mut minor_length = dst1.length();

        // Clamp ellipse eccentricity if too large.
        let adjusted_minor_length = minor_length * self.max_anisotropy;
        if adjusted_minor_length < major_length && minor_length > 0.0 {
            let scale = major_length / adjusted_minor_length;
            dst1 *= scale;
            minor_length *= scale;
        }
        if minor_length == 0.0 {
            return self.triangle(0, st);
        }

        // Choose level of detail for EWA lookup and perform EWA filtering
        let lod = max(0.0, self.levels() as Float - 1.0 + minor_length.log2());
        let i_lod = lod.floor() as usize;

        // NOTE: If we add a bound on T like this `Float: Mul<T, Output=T>`
        // in order to use `lerp()`, the Rust compiler gets confused and won't
        // allow the multiplication between Float values and code in the earlier
        // part of this function will fail to compile.
        //
        // So we do lerp manually :(
        let t = lod - i_lod as Float;
        self.ewa(i_lod, st, &dst0, &dst1) * (1.0 - t) + self.ewa(i_lod + 1, st, &dst0, &dst1) * t
    }

    /// Interpolates using a triangle filter between 4 texels that surround
    /// a given sample point.
    ///
    /// * `level` - The MIPMap level.
    /// * `st`    - The sample point coordinates (s, t).
    fn triangle(&self, level: usize, st: &Point2f) -> T {
        let level = clamp(level, 0, self.levels() - 1);

        let s = st[0] * self.pyramid[level].u_size() as Float - 0.5;
        let t = st[1] * self.pyramid[level].v_size() as Float - 0.5;

        let s0 = s.floor() as usize;
        let t0 = t.floor() as usize;

        let ds = s - s0 as Float;
        let dt = t - t0 as Float;

        let tx0 = texel(&self.pyramid, self.wrap_mode, level, s0, t0);
        let tx1 = texel(&self.pyramid, self.wrap_mode, level, s0, t0 + 1);
        let tx2 = texel(&self.pyramid, self.wrap_mode, level, s0 + 1, t0);
        let tx3 = texel(&self.pyramid, self.wrap_mode, level, s0 + 1, t0 + 1);

        tx0 * (1.0 - ds) * (1.0 - dt)
            + tx1 * (1.0 - ds) * dt
            + tx2 * ds * (1.0 - dt)
            + tx3 * ds * dt
    }

    /// Interpolates using EWA filter between 4 texels that surround a given
    /// sample point.
    ///
    ///
    /// * `level` - The MIPMap level.
    /// * `st`    - The sample point coordinates (s, t).
    /// * `dst0`  - Length of first elliptical axis.
    /// * `dst1`  - Length of second elliptical axis.
    fn ewa(&self, level: usize, st: &Point2f, dst0: &Vector2f, dst1: &Vector2f) -> T {
        let levels = self.levels();
        if level >= levels {
            return texel(&self.pyramid, self.wrap_mode, levels - 1, 0, 0);
        }

        let u_size = self.pyramid[level].u_size();
        let v_size = self.pyramid[level].v_size();

        // Convert EWA coordinates to appropriate scale for level.
        let st = [st[0] * u_size as Float - 0.5, st[1] * v_size as Float - 0.5];
        let dst0 = [dst0[0] * u_size as Float, dst0[1] * v_size as Float];
        let dst1 = [dst1[0] * u_size as Float, dst1[1] * v_size as Float];

        // Compute ellipse coefficients to bound EWA filter region.
        let mut a = dst0[1] * dst0[1] + dst1[1] * dst1[1] + 1.0;
        let mut b = -2.0 * (dst0[0] * dst0[1] + dst1[0] * dst1[1]);
        let mut c = dst0[0] * dst0[0] + dst1[0] * dst1[0] + 1.0;
        let inv_f = 1.0 / (a * c - b * b * 0.25);
        a *= inv_f;
        b *= inv_f;
        c *= inv_f;

        // Compute the ellipse's `(s,t)` bounding box in texture space.
        let det = -b * b + 4.0 * a * c;
        let inv_det = 1.0 / det;
        let u_sqrt = (det * c).sqrt();
        let v_sqrt = (a * det).sqrt();
        let s0 = (st[0] - 2.0 * inv_det * u_sqrt).ceil() as usize;
        let s1 = (st[0] + 2.0 * inv_det * u_sqrt).floor() as usize;
        let t0 = (st[1] - 2.0 * inv_det * v_sqrt).ceil() as usize;
        let t1 = (st[1] + 2.0 * inv_det * v_sqrt).floor() as usize;

        // Scan over ellipse bound and compute quadratic equation.
        let mut sum = T::default();
        let mut sum_wts = 0.0;
        for it in t0..t1 + 1 {
            let tt = it as Float - st[1];
            for is in s0..s1 + 1 {
                let ss = is as Float - st[0];
                // Compute squared radius and filter texel if inside ellipse.
                let r2 = a * ss * ss + b * ss * tt + c * tt * tt;
                if r2 < 1.0 {
                    let index = min(
                        (r2 * WEIGHT_LUT_SIZE as Float) as usize,
                        WEIGHT_LUT_SIZE - 1,
                    );
                    let weight = self.weight_lut[index];
                    sum += texel(&self.pyramid, self.wrap_mode, level, is, it) * weight;
                    sum_wts += weight;
                }
            }
        }
        sum / sum_wts
    }
}

/// Resample the weights for texels at a new resolution of the image size.
///
/// * `old_res` - The old resolution.
/// * `new_res` - The new resolution.
fn resample_weights(old_res: usize, new_res: usize) -> Vec<ResampleWeight> {
    assert!(new_res > old_res);

    let mut wt: Vec<ResampleWeight> = vec![ResampleWeight::default(); new_res];

    let filterwidth = 2.0;
    for i in 0..new_res {
        // Compute image resampling weights for i^th texel.
        let center = (i as Float + 0.5) * old_res as Float / new_res as Float;

        wt[i].first_texel = ((center - filterwidth) + 0.5).floor() as usize;
        for j in 0..4 {
            let pos = wt[i].first_texel as Float + j as Float + 0.5;
            wt[i].weight[j] = lanczos((pos - center) / filterwidth, 2.0);
        }

        // Normalize filter weights for texel resampling
        let inv_sum_wts =
            1.0 / (wt[i].weight[0] + wt[i].weight[1] + wt[i].weight[2] + wt[i].weight[3]);
        for j in 0..4 {
            wt[i].weight[j] *= inv_sum_wts;
        }
    }
    wt
}

/// Returns the texel from the MIPMap pyramid level.
///
/// * `pyramid`   - The MIPMap pyramid.
/// * `wrap_mode` - The image wrap mode.
/// * `level`     - MIPMap Level.
/// * `s`         - s-index.
/// * `t`         - t-index.
fn texel<T>(
    pyramid: &[BlockedArray<T>],
    wrap_mode: ImageWrap,
    level: usize,
    s: usize,
    t: usize,
) -> T
where
    T: Copy + Default,
{
    assert!(level < pyramid.len());

    let l = &pyramid[level];
    let u_size = l.u_size();
    let v_size = l.v_size();

    // Compute texel `(s, t)` accounting for boundary conditions.
    match wrap_mode {
        ImageWrap::Repeat => l[(rem(s, u_size), rem(t, v_size))],
        ImageWrap::Clamp => l[(clamp(s, 0, u_size - 1), clamp(t, 0, v_size - 1))],
        ImageWrap::Black => {
            if s >= u_size || t >= v_size {
                T::default()
            } else {
                l[(s, t)]
            }
        }
    }
}
