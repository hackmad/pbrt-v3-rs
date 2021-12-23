//! BSDF

#![allow(dead_code)]
use super::*;
use crate::interaction::*;
use crate::rng::*;
use bitflags::bitflags;
use bumpalo::Bump;
use std::fmt;

bitflags! {
    /// Stores combinations of reflection models.
    pub struct BxDFType: u8 {
        const BSDF_NONE = 0;
        const BSDF_REFLECTION = 1 << 0;
        const BSDF_TRANSMISSION = 1 << 1;
        const BSDF_DIFFUSE = 1 << 2;
        const BSDF_GLOSSY = 1 << 3;
        const BSDF_SPECULAR = 1 << 4;
    }
}

impl Default for BxDFType {
    /// Returns the "default value" for a `BxDFType`.
    fn default() -> Self {
        Self::BSDF_NONE
    }
}

impl fmt::Display for BxDFType {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        if self.bits & Self::BSDF_REFLECTION.bits == Self::BSDF_REFLECTION.bits {
            s += "BSDF_REFLECTION";
        }
        if self.bits & Self::BSDF_TRANSMISSION.bits == Self::BSDF_TRANSMISSION.bits {
            if s.len() > 0 {
                s += " | ";
            }
            s += "BSDF_TRANSMISSION";
        }
        if self.bits & Self::BSDF_DIFFUSE.bits == Self::BSDF_DIFFUSE.bits {
            if s.len() > 0 {
                s += " | ";
            }
            s += "BSDF_DIFFUSE";
        }
        if self.bits & Self::BSDF_GLOSSY.bits == Self::BSDF_GLOSSY.bits {
            if s.len() > 0 {
                s += " | ";
            }
            s += "BSDF_GLOSSY";
        }
        if self.bits & Self::BSDF_SPECULAR.bits == Self::BSDF_SPECULAR.bits {
            if s.len() > 0 {
                s += " | ";
            }
            s += "BSDF_SPECULAR";
        }
        write!(f, "{}", s)
    }
}

/// Maximum number of BxDFs that can be stored in `BSDF`.
pub const MAX_BXDFS: usize = 8;

/// BSDF interface represents a collection of BRDFs and BTDFs.
#[derive(Clone)]
pub struct BSDF {
    /// The shading normal given by per-vertex normals and/or bump mapping.
    /// It is the first axis in the orthonormal coordinate system and also
    /// used to define hemispheres for integrating incident illumincation for
    /// surface reflection.
    pub ns: Normal3f,

    /// The geometric normal defined by surface geometry.
    pub ng: Normal3f,

    /// Second axis for the orthonormal coordinate system.
    pub ss: Vector3f,

    /// Third axis for the orthonormal coordinate system.
    pub ts: Vector3f,

    /// The `BxDFs`.
    pub bxdfs: Vec<BxDF>,

    /// Relative index of refraction over the surfaceboundary.
    pub eta: Float,
}

impl BSDF {
    /// Creates a new `BSDF`.
    ///
    /// * `si`  - The differential geometry at the point on a surface.
    /// * `eta` - Optional relative index of refraction over the surface
    ///           boundary. If not provided, defaults to 1.0; used for
    ///           opaque surfaces.
    pub fn new(si: &SurfaceInteraction, eta: Option<Float>) -> Self {
        let eta = eta.map_or_else(|| 1.0, |e| e);
        let ns = si.shading.n;
        let ss = si.shading.dpdu.normalize();

        Self {
            eta,
            ns,
            ng: si.hit.n,
            ss,
            ts: Vector3::from(ns).cross(&ss),
            bxdfs: Vec::with_capacity(MAX_BXDFS),
        }
    }

    /// Allocates a new `BSDF`.
    ///
    /// * `si`  - The differential geometry at the point on a surface.
    /// * `eta` - Optional relative index of refraction over the surface
    ///           boundary. If not provided, defaults to 1.0; used for
    ///           opaque surfaces.
    pub fn alloc(allocator: &Bump, si: &SurfaceInteraction, eta: Option<Float>) -> Self {
        allocator.alloc(Self::new(si, eta)).to_owned()
    }

    /// Add a `BxDF`.
    ///
    /// * `bxdf` - The `BxDF`.
    pub fn add(&mut self, bxdf: BxDF) {
        assert!(
            self.bxdfs.len() < MAX_BXDFS,
            "Cannot add BxDFs. BSDF maximum limit {} reached.",
            MAX_BXDFS
        );
        self.bxdfs.push(bxdf);
    }

    /// Returns the number of `BxDF`s that match the given type.
    ///
    /// * `bxdf_type` - The `BxdFType` to match (default to `BxDFType::all()`).
    pub fn num_components(&self, bxdf_type: BxDFType) -> usize {
        let mut num = 0;
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                num += 1;
            }
        }
        num
    }

    /// Transforms a vector from world space to local space.
    ///
    /// * `v` - The vector to transform.
    pub fn world_to_local(&self, v: &Vector3f) -> Vector3f {
        Vector3f::new(v.dot(&self.ss), v.dot(&self.ts), v.dot(&self.ns))
    }

    /// Transforms a vector from local space to world space.
    ///
    /// * `v` - The vector to transform.
    pub fn local_to_world(&self, v: &Vector3f) -> Vector3f {
        Vector3f::new(
            self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z,
            self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z,
            self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z,
        )
    }

    /// Returns the BSDF evaluated for a pair of directions.
    ///
    /// * `wo_w`      - Outgoing direction in world-space.
    /// * `wi_w`      - Incident direction in world-space.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn f(&self, wo_w: &Vector3f, wi_w: &Vector3f, bxdf_type: BxDFType) -> Spectrum {
        let wi = self.world_to_local(wi_w);
        let wo = self.world_to_local(wo_w);

        if wo.z == 0.0 {
            Spectrum::new(0.0)
        } else {
            let reflect = wi_w.dot(&self.ng) * wo_w.dot(&self.ng) > 0.0;
            let mut f = Spectrum::new(0.0);
            for bxdf in self.bxdfs.iter() {
                let curr_type = bxdf.get_type();
                if bxdf.matches_flags(bxdf_type)
                    && ((reflect && curr_type.contains(BxDFType::BSDF_REFLECTION))
                        || (!reflect && curr_type.contains(BxDFType::BSDF_TRANSMISSION)))
                {
                    f += bxdf.f(&wo, &wi);
                }
            }
            f
        }
    }

    /// Returns the value of the BSDF given the outgpoing direction.
    /// direction.
    ///
    /// * `wo_world`  - Outgoing direction in world-space.
    /// * `u`         - The 2D uniform random values.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn sample_f(&self, wo_world: &Vector3f, u: &Point2f, bxdf_type: BxDFType) -> BxDFSample {
        // Choose which `BxDF` to sample.
        let matching_comps = self.num_components(bxdf_type);
        if matching_comps == 0 {
            debug!("For wo_world = {}, matching_comps = 0", wo_world);
            return BxDFSample::default();
        }
        let comp = min(
            (u[0] * matching_comps as Float).floor() as usize,
            matching_comps - 1,
        );

        // Get BxDF for chosen component.
        let mut count = comp;
        let mut bxdf_idx: Option<usize> = None;
        for (i, bxdf) in self.bxdfs.iter().enumerate() {
            if bxdf.matches_flags(bxdf_type) {
                if count == 0 {
                    bxdf_idx = Some(i);
                    break;
                }
                count -= 1;
            }
        }
        let bxdf_idx = bxdf_idx.expect("bsdf::sample_f() did not find matching bxdf");
        let matched_bxdf = &self.bxdfs[bxdf_idx];
        debug!(
            "BSDF::Sample_f chose comp = {} / matching = {}, bxdf: {}",
            comp,
            matching_comps,
            matched_bxdf.get_type()
        );

        // Remap BxDF sample `u` to `[0,1)^2`.
        let u_remapped = Point2f::new(
            min(
                u[0] * matching_comps as Float - comp as Float,
                ONE_MINUS_EPSILON,
            ),
            u[1],
        );

        // Sample chosen `BxDF`.
        let wo = self.world_to_local(wo_world);
        if wo.z == 0.0 {
            debug!("For wo_world = {}, wo = {}, wo.z = 0", wo_world, wo);
            return BxDFSample::default();
        }

        let mut sample = matched_bxdf.sample_f(&wo, &u_remapped);
        debug!(
            "For wo_world = {}, wo = {}, sampled f = {}, pdf = {}, ratio = {}, wi = {}",
            wo_world,
            wo,
            sample.f,
            sample.pdf,
            if sample.pdf > 0.0 {
                sample.f / sample.pdf // Potential for NaN if sample.pdf == 0
            } else {
                Spectrum::new(0.0)
            },
            sample.wi
        );

        if sample.pdf == 0.0 {
            debug!("sample.pdf = 0");
            return BxDFSample::default();
        }
        sample.bxdf_type = matched_bxdf.get_type();

        let wi_world = self.local_to_world(&sample.wi);

        // Compute overall PDF with all matching BxDFs.
        if sample.bxdf_type & BxDFType::BSDF_SPECULAR == BxDFType::BSDF_NONE && matching_comps > 1 {
            for (i, bxdf) in self.bxdfs.iter().enumerate() {
                if bxdf_idx != i && bxdf.matches_flags(bxdf_type) {
                    sample.pdf += bxdf.pdf(&wo, &sample.wi);
                }
            }
        }
        if matching_comps > 1 {
            sample.pdf /= matching_comps as Float;
        }

        // Compute value of BSDF for sampled direction.
        if sample.bxdf_type & BxDFType::BSDF_SPECULAR == BxDFType::BSDF_NONE {
            let reflect = wi_world.dot(&self.ng) * wo_world.dot(&self.ng) > 0.0;
            sample.f = Spectrum::new(0.0);
            for bxdf in self.bxdfs.iter() {
                if bxdf.matches_flags(bxdf_type)
                    && ((reflect
                        && (bxdf.get_type() & BxDFType::BSDF_REFLECTION != BxDFType::BSDF_NONE))
                        || (!reflect
                            && (bxdf.get_type() & BxDFType::BSDF_TRANSMISSION
                                != BxDFType::BSDF_NONE)))
                {
                    sample.f += bxdf.f(&wo, &sample.wi);
                }
            }
        }

        debug!(
            "Overall f = {}, pdf = {}, ratio = {}",
            sample.f,
            sample.pdf,
            if sample.pdf > 0.0 {
                sample.f / sample.pdf
            } else {
                Spectrum::new(0.0)
            }
        );

        // We need to return wi in world-space.
        sample.wi = wi_world;
        sample
    }

    /// Computes the hemispherical-directional reflectance function ρ.
    ///
    /// * `wo_world`  - Outgoing direction in world-space.
    /// * `u`         - Samples used by Monte Carlo algorithm.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn rho_hd(&self, wo_world: &Vector3f, u: &[Point2f], bxdf_type: BxDFType) -> Spectrum {
        let wo = self.world_to_local(wo_world);
        let mut l = Spectrum::new(0.0);
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                l += bxdf.rho_hd(&wo, u);
            }
        }
        l
    }

    /// Computes the hemispherical-hemispherical-directional reflectance function ρ.
    ///
    /// * `u1`        - Samples used b Monte Carlo algorithm.
    /// * `u2`        - Samples used b Monte Carlo algorithm.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn rho_hh(&self, u1: &[Point2f], u2: &[Point2f], bxdf_type: BxDFType) -> Spectrum {
        let mut l = Spectrum::new(0.0);
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                l += bxdf.rho_hh(u1, u2);
            }
        }
        l
    }

    /// Evaluates the PDF for the sampling method. Default is based on the
    /// cosine-weighted sampling in `BxDF::sample_f()` default implementation.
    ///
    /// * `wo_world`  - Outgoing direction in world-space.
    /// * `wi_world`  - Incident direction in world-space.
    /// * `bxdf_type` - The `BxdFType` to evaluate.
    pub fn pdf(&self, wo_world: &Vector3f, wi_world: &Vector3f, bxdf_type: BxDFType) -> Float {
        if self.bxdfs.len() == 0 {
            return 0.0;
        }

        let wo = self.world_to_local(wo_world);
        let wi = self.world_to_local(wi_world);

        if wo.z == 0.0 {
            return 0.0;
        }

        let mut matching_comps = 0;
        let mut pdf = 0.0;
        for bxdf in self.bxdfs.iter() {
            if bxdf.matches_flags(bxdf_type) {
                matching_comps += 1;
                pdf += bxdf.pdf(&wo, &wi);
            }
        }
        if matching_comps > 0 {
            pdf / matching_comps as Float
        } else {
            0.0
        }
    }
}
