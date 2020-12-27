//! BxDF Type

/// Types of BSDF models.
pub const BSDF_NONE: u8 = 0b00000000;
pub const BSDF_REFLECTION: u8 = 0b00000001;
pub const BSDF_TRANSMISSION: u8 = 0b00000010;
pub const BSDF_DIFFUSE: u8 = 0b00000100;
pub const BSDF_GLOSSY: u8 = 0b00001000;
pub const BSDF_SPECULAR: u8 = 0b00010000;
pub const BSDF_ALL: u8 = 0b00011111;

/// Stores combinations of reflection models.
#[repr(C)]
#[derive(Copy, Clone, Default)]
pub struct BxDFType {
    t: u8,
}

impl BxDFType {
    /// Tests a single type flag and returns whether it is set or not.
    ///
    /// * `flag` - BxDFType flag.
    pub fn matches(&self, flag: u8) -> bool {
        self.t & flag > 0
    }
}

impl PartialEq for BxDFType {
    /// Returns true if the reflection models match.
    ///
    /// * `other` - The reflection model to compare.
    fn eq(&self, other: &Self) -> bool {
        self.t & other.t == self.t
    }
}

impl From<u8> for BxDFType {
    /// Convert a `u8` value to `BxDFType`.
    ///
    /// * `t` - A `u8` value containing combination of `BXDF_*` flags combined
    ///         bitwise OR operator.
    fn from(t: u8) -> Self {
        assert!(t <= BSDF_ALL, "Invalid BxDF flags {}=({:#08b})", t, t);
        Self { t }
    }
}
