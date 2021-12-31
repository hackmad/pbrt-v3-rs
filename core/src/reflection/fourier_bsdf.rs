//! Fourier Basis BSDF Model

use super::*;
use crate::interpolation::*;
use crate::material::*;
use bumpalo::Bump;
use std::sync::Arc;

/// BSDF for modeling materials like metals with smooth or rough coatings and
/// fabrics which are often partially retro-reflective.
#[derive(Clone)]
pub struct FourierBSDF {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// The BSDF data.
    bsdf_table: Arc<FourierBSDFTable>,

    /// Indicates whether incident ray started from a light source or from camera.
    mode: TransportMode,
}

impl FourierBSDF {
    /// Allocate a new instance of `FourierBSDF`.
    ///
    /// * `arena`      - The arena for memory allocations.
    /// * `bsdf_table` - The BSDF data.
    /// * `mode`       - Indicates whether incident ray started from a light source
    ///                  or from camera.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc<'arena>(
        arena: &'arena Bump,
        bsdf_table: Arc<FourierBSDFTable>,
        mode: TransportMode,
    ) -> &'arena mut BxDF {
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_REFLECTION
                | BxDFType::BSDF_TRANSMISSION
                | BxDFType::BSDF_GLOSSY,
            bsdf_table: Arc::clone(&bsdf_table),
            mode,
        });
        arena.alloc(BxDF::FourierBSDF(model))
    }

    /// Clone into a newly allocated a new instance of `FourierBSDF`.
    ///
    /// * `arena` - The memory arena.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc<'arena>(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            bsdf_table: Arc::clone(&self.bsdf_table),
            mode: self.mode,
        });
        arena.alloc(BxDF::FourierBSDF(model))
    }

    /// Returns the BxDF type.
    pub fn get_type(&self) -> BxDFType {
        self.bxdf_type
    }

    /// Returns the value of the distribution function for the given pair of
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `wi` - Incident direction.
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        // Find the zenith angle cosines and azimuth difference angle.
        let mu_i = cos_theta(&(-wi));
        let mu_o = cos_theta(wo);
        let cos_phi = cos_d_phi(&(-wi), wo) as f64;

        // Compute Fourier coefficients `ak` for `(μi, μo)`.

        // Determine offsets and weights for `(μi, μo)`.
        let (weights_i, offset_i) = match self.bsdf_table.get_weights_and_offset(mu_i) {
            Some((weights, offset)) => (weights, offset),
            None => return Spectrum::new(0.0),
        };
        let (weights_o, offset_o) = match self.bsdf_table.get_weights_and_offset(mu_o) {
            Some((weights, offset)) => (weights, offset),
            None => return Spectrum::new(0.0),
        };

        // Allocate storage to accumulate `ak` coefficients.
        let ak_size = self.bsdf_table.m_max * self.bsdf_table.n_channels;
        let mut ak = vec![0.0; ak_size];

        // Accumulate weighted sums of nearby `ak` coefficients.
        let mut m_max = 0;
        for (b, wtob) in weights_o.iter().enumerate() {
            for (a, wtia) in weights_i.iter().enumerate() {
                // Add contribution of `(a, b)` to `ak` values.
                let weight = wtia * wtob;
                if weight != 0.0 {
                    let (m, ap) = self.bsdf_table.get_ak(offset_i + a, offset_o + b);
                    m_max = max(m_max, m);
                    for c in 0..self.bsdf_table.n_channels {
                        for k in 0..m {
                            ak[c * self.bsdf_table.m_max + k] += weight * ap[c * m + k];
                        }
                    }
                }
            }
        }

        // Evaluate Fourier expansion for angle ϕ.
        let y = max(0.0, fourier(&ak[0..m_max], cos_phi));
        let mut scale = if mu_i != 0.0 { 1.0 / abs(mu_i) } else { 0.0 };

        // Update `scale` to account for adjoint light transport.
        if self.mode == TransportMode::Radiance && mu_i * mu_o > 0.0 {
            let eta = if mu_i > 0.0 {
                1.0 / self.bsdf_table.eta
            } else {
                self.bsdf_table.eta
            };
            scale *= eta * eta;
        }
        if self.bsdf_table.n_channels == 1 {
            Spectrum::new(y * scale)
        } else {
            // Compute and return RGB colors for tabulated BSDF.
            let rs = self.bsdf_table.m_max;
            let re = rs + m_max;
            let r = fourier(&ak[rs..re], cos_phi);

            let bs = 2 * self.bsdf_table.m_max;
            let be = bs + m_max;
            let b = fourier(&ak[bs..be], cos_phi);

            let g = 1.39829 * y - 0.100913 * b - 0.297375 * r;
            let rgb = [r * scale, g * scale, b * scale];

            Spectrum::from_rgb(&rgb, None).clamp_default()
        }
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        // Sample zenith angle component for _FourierBSDF_
        let mu_o = cos_theta(wo);
        let (mu_i, _, pdf_mu) = sample_catmull_rom_2d(
            &self.bsdf_table.mu,
            &self.bsdf_table.mu,
            &self.bsdf_table.a0,
            &self.bsdf_table.cdf,
            mu_o,
            u[1],
        );

        // Compute Fourier coefficients `ak` for `(μi, μo)`.

        // Determine offsets and weights for `(μi, μo)`.
        let (weights_i, offset_i) = match self.bsdf_table.get_weights_and_offset(mu_i) {
            Some((weights, offset)) => (weights, offset),
            None => return BxDFSample::from(self.bxdf_type),
        };
        let (weights_o, offset_o) = match self.bsdf_table.get_weights_and_offset(mu_o) {
            Some((weights, offset)) => (weights, offset),
            None => return BxDFSample::from(self.bxdf_type),
        };

        // Allocate storage to accumulate `ak` coefficients.
        let ak_size = self.bsdf_table.m_max * self.bsdf_table.n_channels;
        let mut ak = vec![0.0; ak_size];

        // Accumulate weighted sums of nearby `a_k` coefficients.
        let mut m_max = 0;
        for (b, wtob) in weights_o.iter().enumerate() {
            for (a, wtia) in weights_i.iter().enumerate() {
                // Add contribution of `(a, b)` to `ak` values.
                let weight = wtia * wtob;
                if weight != 0.0 {
                    let (m, ap) = self.bsdf_table.get_ak(offset_i + a, offset_o + b);
                    m_max = max(m_max, m);
                    for c in 0..self.bsdf_table.n_channels {
                        for k in 0..m {
                            ak[c * self.bsdf_table.m_max + k] += weight * ap[c * m + k];
                        }
                    }
                }
            }
        }

        // Importance sample the luminance Fourier expansion.
        let (y, pdf_phi, phi) = sample_fourier(&ak[0..m_max], &self.bsdf_table.recip, u[0]);
        let pdf = max(0.0, pdf_phi * pdf_mu);

        // Compute the scattered direction for `FourierBSDF`.
        let sin2_theta_i = max(0.0, 1.0 - mu_i * mu_i);
        let mut norm = (sin2_theta_i / sin_2_theta(wo)).sqrt();
        if norm.is_infinite() {
            norm = 0.0;
        }

        let sin_phi = sin(phi);
        let cos_phi = cos(phi);
        let mut wi = -Vector3f::new(
            norm * (cos_phi * wo.x - sin_phi * wo.y),
            norm * (sin_phi * wo.x + cos_phi * wo.y),
            mu_i,
        );

        // Mathematically, `wi` will be normalized (if `wo` was). However, in
        // practice, floating-point rounding error can cause some error to
        // accumulate in the computed value of wi here. This can be catastrophic:
        // if the ray intersects an object with the FourierBSDF again and the `wo`
        // (based on such a `wi`) is nearly perpendicular to the surface, then the
        // `wi` computed at the next intersection can end up being substantially
        // (like 4x) longer than normalized, which leads to all sorts of errors,
        // including negative spectral values. Therefore, we normalize again here.
        wi = wi.normalize();

        // Evaluate remaining Fourier expansions for angle ϕ.
        let mut scale = if mu_i != 0.0 { 1.0 / abs(mu_i) } else { 0.0 };
        if self.mode == TransportMode::Radiance && mu_i * mu_o > 0.0 {
            let eta = if mu_i > 0.0 {
                1.0 / self.bsdf_table.eta
            } else {
                self.bsdf_table.eta
            };
            scale *= eta * eta;
        }

        if self.bsdf_table.n_channels == 1 {
            BxDFSample::new(Spectrum::new(y * scale), pdf, wi, self.bxdf_type)
        } else {
            let rs = self.bsdf_table.m_max;
            let re = rs + m_max;
            let r = fourier(&ak[rs..re], cos_phi as f64);

            let bs = 2 * self.bsdf_table.m_max;
            let be = bs + m_max;
            let b = fourier(&ak[bs..be], cos_phi as f64);

            let g = 1.39829 * y - 0.100913 * b - 0.297375 * r;
            let rgb = [r * scale, g * scale, b * scale];

            let spectrum = Spectrum::from_rgb(&rgb, None).clamp_default();
            BxDFSample::new(spectrum, pdf, wi, self.bxdf_type)
        }
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // Find the zenith angle cosines and azimuth difference angle.
        let mu_i = cos_theta(&(-(wi)));
        let mu_o = cos_theta(wo);
        let cos_phi = cos_d_phi(&(-(wi)), wo);

        // Compute luminance Fourier coefficients `a_k` for `(μi, μo)`.
        let (weights_i, offset_i) = match self.bsdf_table.get_weights_and_offset(mu_i) {
            Some((weights, offset)) => (weights, offset),
            None => return 0.0,
        };
        let (weights_o, offset_o) = match self.bsdf_table.get_weights_and_offset(mu_o) {
            Some((weights, offset)) => (weights, offset),
            None => return 0.0,
        };

        let ak_size = self.bsdf_table.m_max * self.bsdf_table.n_channels;
        let mut ak = vec![0.0; ak_size];

        let mut m_max = 0;
        for (o, wto) in weights_o.iter().enumerate() {
            for (i, wti) in weights_i.iter().enumerate() {
                // Add contribution of `(a, b)` to `ak` values.
                let weight = wti * wto;
                if weight == 0.0 {
                    continue;
                }

                let (order, coeffs) = self.bsdf_table.get_ak(offset_i + i, offset_o + o);
                m_max = max(m_max, order);

                for k in 0..order {
                    ak[k] += coeffs[k] * weight;
                }
            }
        }

        // Evaluate probability of sampling `wi`.
        let n_mu = self.bsdf_table.mu.len();
        let rho = (0..4).fold(0.0, |a, o| {
            if weights_o[o] == 0.0 {
                a
            } else {
                a + weights_o[o] * self.bsdf_table.cdf[(offset_o + o) * n_mu + n_mu - 1] * TWO_PI
            }
        });

        let y = fourier(&ak[0..m_max], cos_phi as f64);
        if rho > 0.0 && y > 0.0 {
            y / rho
        } else {
            0.0
        }
    }
}
