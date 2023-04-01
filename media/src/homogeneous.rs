//! Homogeneous Medium

use core::geometry::*;
use core::interaction::MediumInteraction;
use core::medium::HenyeyGreenstein;
use core::medium::Medium;
use core::pbrt::*;
use core::sampler::*;
use core::spectrum::*;
use std::sync::Arc;

/// Implements a homogeneous medium representing a region of space with constant
/// σa and σs values throughout its extent.
pub struct HomogeneousMedium {
    /// Absorption cross section `σa` is the probability density that light is
    /// absorbed per unit distance traveled in the medium
    _sigma_a: Spectrum,

    /// Scattering coefficient `σs` is the probability of an out-scattering
    /// event occurring per unit distance
    sigma_s: Spectrum,

    /// Total reduction in radiance due to absorption and out-scattering
    /// `σt = σs + σa`. This combined effect of absorption and out-scattering is
    /// called attenuation or extinction.
    sigma_t: Spectrum,

    /// The asymmetry parameter for Henyey-Greenstein phase function.
    g: Float,
}

impl HomogeneousMedium {
    /// Create a new `HomogeneousMedium `.
    ///
    /// * `sigma_a` - Absorption cross section `σa`.
    /// * `sigma_s` - Scattering coefficient `σs`.
    /// * `g`       - The asymmetry parameter for Henyey-Greenstein phase
    ///               function.
    pub fn new(sigma_a: Spectrum, sigma_s: Spectrum, g: Float) -> Self {
        Self {
            _sigma_a: sigma_a,
            sigma_s,
            sigma_t: sigma_s + sigma_a,
            g,
        }
    }
}

impl Medium for HomogeneousMedium {
    /// Returns the beam transmittance along a given ray.
    ///
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn tr(&self, ray: &Ray, _sampler: &mut ArcSampler) -> Spectrum {
        (-self.sigma_t * min(ray.t_max * ray.d.length(), Float::MAX)).exp()
    }

    /// Samples a medium scattering interaction along a world-space ray.
    ///
    /// The ray will generally have been intersected against the scene geometry;
    /// thus, implementations of this method shouldn’t ever sample a medium
    /// interaction at a point on the ray beyond its `t_max` value.
    ///
    /// NOTE: Calling code will need to assign this medium as we cannot pass
    /// back and `ArcMedium` out of here for `Self`.
    ///
    /// * `ray`     - The ray.
    /// * `sampler` - The sampler.
    fn sample(&self, ray: &Ray, sampler: &mut ArcSampler) -> (Spectrum, Option<MediumInteraction>) {
        // Sample a channel and distance along the ray.
        let sampler = Arc::get_mut(sampler).unwrap();
        let channel = min(
            (sampler.get_1d() * SPECTRUM_SAMPLES as Float) as usize,
            SPECTRUM_SAMPLES - 1,
        );
        let dist = -(1.0 - sampler.get_1d()).ln() / self.sigma_t[channel];
        let t = min(dist / ray.d.length(), ray.t_max);
        let sampled_medium = t < ray.t_max;

        let mi = if sampled_medium {
            let phase_function = HenyeyGreenstein::new(self.g);
            Some(MediumInteraction::new(
                ray.at(t),
                -ray.d,
                ray.time,
                None,
                phase_function,
            ))
        } else {
            None
        };

        // Compute the transmittance and sampling density
        let tr = (-self.sigma_t * min(t, Float::MAX) * ray.d.length()).exp();

        // Return weighting factor for scattering from homogeneous medium
        let density = if sampled_medium { self.sigma_t * tr } else { tr };

        let mut pdf = 0.0;
        for i in 0..SPECTRUM_SAMPLES {
            pdf += density[i];
        }

        pdf *= 1.0 / SPECTRUM_SAMPLES as Float;
        if pdf == 0.0 {
            assert!(tr.is_black());
            pdf = 1.0;
        }

        let s = if sampled_medium {
            tr * self.sigma_s / pdf
        } else {
            tr / pdf
        };

        (s, mi)
    }
}
