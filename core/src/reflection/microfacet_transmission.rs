//! Torrance-Sparrow Microfacet Transmission Model

#![allow(dead_code)]

use super::*;
use crate::material::*;
use crate::microfacet::*;
use bumpalo::Bump;
use std::fmt;

/// BTDF for modeling glossy transmissive surfaces using a microfacet distribution.
pub struct MicrofacetTransmission<'arena> {
    /// BxDF type.
    bxdf_type: BxDFType,

    /// Fresnel interface for dielectrics.
    fresnel: &'arena mut Fresnel<'arena>,

    /// Spectrum used to scale the transmitted colour.
    t: Spectrum,

    /// Index of refraction above the surface (same side as surface normal).
    eta_a: Float,

    /// Index of refraction below the surface (opposite side as surface normal).
    eta_b: Float,

    /// Indicates whether incident ray started from a light source or from camera.
    mode: TransportMode,

    /// The microfacet distribution model.
    distribution: &'arena mut MicrofacetDistribution<'arena>,
}

impl<'arena> MicrofacetTransmission<'arena> {
    /// Allocate a new instance of `MicrofacetTransmission`.
    ///
    /// * `arena`        - The arena for memory allocations.
    /// * `t`            - Spectrum used to scale the transmitted colour.
    /// * `distribution` - Microfacet distribution.
    /// * `eta_a`        - Index of refraction above the surface (same side as
    ///                    surface normal).
    /// * `eta_b`        - Index of refraction below the surface (opposite side
    ///                    as surface normal).
    /// * `mode`         - Indicates whether incident ray started from a light
    ///                    source or from camera.
    #[allow(clippy::mut_from_ref)]
    pub fn alloc(
        arena: &'arena Bump,
        t: Spectrum,
        distribution: &'arena mut MicrofacetDistribution<'arena>,
        eta_a: Float,
        eta_b: Float,
        mode: TransportMode,
    ) -> &'arena mut BxDF<'arena> {
        let fresnel = FresnelDielectric::alloc(arena, eta_a, eta_b);
        let model = arena.alloc(Self {
            bxdf_type: BxDFType::BSDF_TRANSMISSION | BxDFType::BSDF_GLOSSY,
            distribution,
            fresnel,
            t,
            eta_a,
            eta_b,
            mode,
        });
        arena.alloc(BxDF::MicrofacetTransmission(model))
    }

    /// Clone into a newly allocated a new instance of `MicrofacetTransmission`.
    ///
    /// * `arena` - The arena for memory allocations.
    #[allow(clippy::mut_from_ref)]
    pub fn clone_alloc(&self, arena: &'arena Bump) -> &'arena mut BxDF<'arena> {
        let distribution = self.distribution.clone_alloc(arena);
        let fresnel = self.fresnel.clone_alloc(arena);
        let model = arena.alloc(Self {
            bxdf_type: self.bxdf_type,
            distribution,
            fresnel,
            t: self.t,
            eta_a: self.eta_a,
            eta_b: self.eta_b,
            mode: self.mode,
        });
        arena.alloc(BxDF::MicrofacetTransmission(model))
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
        if same_hemisphere(wo, wi) {
            Spectrum::ZERO // transmission only
        } else {
            let cos_theta_o = cos_theta(wo);
            let cos_theta_i = cos_theta(wi);
            if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
                Spectrum::ZERO
            } else {
                // Compute `wh` from `wo` and `wi` for microfacet transmission.
                let eta = if cos_theta(wo) > 0.0 {
                    self.eta_b / self.eta_a
                } else {
                    self.eta_a / self.eta_b
                };
                let mut wh = (wo + wi * eta).normalize();
                if wh.z < 0.0 {
                    wh = -wh;
                }

                // Same side?
                if wo.dot(&wh) * wi.dot(&wh) > 0.0 {
                    Spectrum::ZERO
                } else {
                    let f = self.fresnel.evaluate(wo.dot(&wh));

                    let sqrt_denom = wo.dot(&wh) + eta * wi.dot(&wh);
                    let factor = if self.mode == TransportMode::Radiance {
                        1.0 / eta
                    } else {
                        1.0
                    };

                    (Spectrum::ONE - f)
                        * self.t
                        * abs(self.distribution.d(&wh)
                            * self.distribution.g(wo, wi)
                            * eta
                            * eta
                            * wi.abs_dot(&wh)
                            * wo.abs_dot(&wh)
                            * factor
                            * factor
                            / (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom))
                }
            }
        }
    }

    /// Returns the value of the BxDF given the outgpoing direction.
    /// directions.
    ///
    /// * `wo` - Outgoing direction.
    /// * `u`  - The 2D uniform random values.
    pub fn sample_f(&self, wo: &Vector3f, u: &Point2f) -> BxDFSample {
        if wo.z == 0.0 {
            BxDFSample::from(self.bxdf_type)
        } else {
            let wh = self.distribution.sample_wh(wo, u);
            if wo.dot(&wh) < 0.0 {
                // Should be rare.
                BxDFSample::from(self.bxdf_type)
            } else {
                let eta = if cos_theta(wo) > 0.0 {
                    self.eta_a / self.eta_b
                } else {
                    self.eta_b / self.eta_a
                };
                if let Some(wi) = refract(wo, &wh.into(), eta) {
                    let pdf = self.pdf(wo, &wi);
                    BxDFSample::new(self.f(wo, &wi), pdf, wi, self.bxdf_type)
                } else {
                    BxDFSample::from(self.bxdf_type)
                }
            }
        }
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if same_hemisphere(wo, wi) {
            0.0
        } else {
            // Compute `wh` from `wo` and `wi` for microfacet transmission.
            let eta = if cos_theta(wo) > 0.0 {
                self.eta_b / self.eta_a
            } else {
                self.eta_a / self.eta_b
            };
            let wh = (wo + wi * eta).normalize();

            if wo.dot(&wh) * wi.dot(&wh) > 0.0 {
                0.0
            } else {
                // Compute change of variables dwh\dwi for microfacet transmission.
                let sqrt_denom = wo.dot(&wh) + eta * wi.dot(&wh);
                let dwh_dwi = abs((eta * eta * wi.dot(&wh)) / (sqrt_denom * sqrt_denom));
                self.distribution.pdf(wo, &wh) * dwh_dwi
            }
        }
    }
}

impl<'arena> fmt::Display for MicrofacetTransmission<'arena> {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "MicrofacetTransmission {{ bxdf_type: {}, fresnel: {}, t: {}, eta_a: {}, eta_b: {}, mode: {}, distribution: {} }}",
            self.bxdf_type, self.fresnel, self.t, self.eta_a, self.eta_b, self.mode, self.distribution
        )
    }
}
