//! Fourier BSDF Table

use super::bsdf_reader::*;
use crate::interpolation::*;
use crate::pbrt::*;
use std::str;

/// Stores the measured Fourier BSDF data.
#[derive(Clone, Debug)]
pub struct FourierBSDFTable {
    /// Relative index of refraction over the surface boundary between two media.
    pub eta: Float,

    /// Maximum order `m` for any pair of μi, μj directions used to allocate
    /// ak coefficients.
    pub m_max: usize,

    /// Number of spectral channels available:
    /// 1: Monochromatic BSDF
    /// 3: BSDF with RGB colors (stored as luminence, red, blue values)
    pub n_channels: usize,

    /// Zenith angles stored in sorted order from low to high.
    pub mu: Vec<Float>,

    /// The order, m, of the Fourier representation bounded by `m_max` that
    /// varies with respect to incident and outgoing zenith angle cosine μi, μo
    /// Number of orders needed is determined by querying this mu.len() x mu.len()
    // integer matrix.
    pub m: Vec<usize>,

    /// The coefficients for all pairs of discretized direction `mu`.
    pub a: Vec<Float>,

    /// Stores offsets into `a`. The `m` coefficients starting at `a[offset] give
    // ak values for the corresponding pairs of directions.
    ///
    /// For the 3 color channel case, the first `m` coefficients after `a[offset]`
    /// encode coefficients for luminance, the next `m` correspond to the red
    /// channel, and then blue follows.
    pub a_offset: Vec<usize>,

    /// First coefficient a0.
    pub a0: Vec<Float>,

    /// CDF values.
    pub cdf: Vec<Float>,

    /// Contains 1 / i for i in [0..`m_max`].
    pub recip: Vec<Float>,
}

impl FourierBSDFTable {
    /// Loads a `FourierBSDF` from a binary file.
    ///
    /// * `path` - The path to the BSDF binary file.
    pub fn from_file(path: &str) -> Result<Self, String> {
        let mut file = open_file(path)?;
        file.check_header()?;

        let flags = file.read_i32()?;
        let n_mu = file.read_i32()? as usize;
        let n_coeffs = file.read_i32()? as usize;
        let m_max = file.read_i32()? as usize;
        let n_channels = file.read_i32()? as usize;
        let n_bases = file.read_i32()? as usize;
        let _unused = file.read_i32_vec(3)?;
        let eta = file.read_f32()?;
        let _unused = file.read_i32_vec(4)?;

        // Only a subset of BSDF files are supported for simplicity. In particular
        // monochromatic and RGB files with uniform (i.e. non-textured) material
        // properties.
        if flags != 1 || (n_channels != 1 && n_channels != 3) || n_bases != 1 {
            return Err(String::from("Unsupported BSDF file format"));
        }

        let mu = file.read_f32_vec(n_mu)?;
        let cdf = file.read_f32_vec(n_mu * n_mu)?;
        let offset_and_length = file.read_i32_vec(n_mu * n_mu * 2)?;
        let a = file.read_f32_vec(n_coeffs)?;

        let mut a0 = vec![0.0; n_mu * n_mu];
        let mut a_offset = vec![0_usize; n_mu * n_mu];
        let mut m = vec![0_usize; n_mu * n_mu];

        for i in 0..n_mu * n_mu {
            let offset = offset_and_length[2 * i] as usize;
            let length = offset_and_length[2 * i + 1] as usize;

            a_offset[i] = offset;
            m[i] = length;
            a0[i] = if length > 0 { a[offset] } else { 0.0 };
        }

        let recip: Vec<Float> = (0..m_max).map(|i| 1.0 / i as Float).collect();

        Ok(Self {
            eta,
            m_max,
            n_channels,
            mu,
            m,
            a,
            a_offset,
            a0,
            cdf,
            recip,
        })
    }

    /// For offsets into the `mu` array for incident and outgoing direction cosines,
    /// returns the order `m` of coefficients for them and their coefficients in `a`.
    ///
    /// * `offset_i` - Offset for incident direction.
    /// * `offset_o` - Offset for outgoing direction.
    pub fn get_ak(&self, offset_i: usize, offset_o: usize) -> (usize, &[Float]) {
        let offset = offset_o * self.mu.len() + offset_i;
        let mptr = self.m[offset];

        let a_offset = self.a_offset[offset];
        let ak = &self.a[a_offset..];

        (mptr, ak)
    }

    /// Returns Catmull-Rom weights and index offset for a given zenith angle.
    ///
    /// * `cos_theta` - The zenith angle to interpolate from `mu`.
    pub fn get_weights_and_offset(&self, cos_theta: Float) -> Option<([Float; 4], usize)> {
        catmull_rom_weights(&self.mu, cos_theta)
    }
}
