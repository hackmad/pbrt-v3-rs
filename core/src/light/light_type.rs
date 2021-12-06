//! Light Types

use bitflags::bitflags;

bitflags! {
    /// Stores combination of flags for the light types.
    pub struct LightType: u8 {
        const DELTA_POSITION_LIGHT = 1;
        const DELTA_DIRECTION_LIGHT = 2;
        const AREA_LIGHT = 4;
        const INFINITE_LIGHT = 8;
    }
}

impl LightType {
    /// Tests a single light type flag and returns whether it is set or not.
    ///
    /// * `other` - Light type flag to match.
    pub fn matches(&self, other: Self) -> bool {
        self.bits & other.bits > 0
    }

    /// Returns true if the light flags has the DELTA_POSITION_LIGHT or
    /// DELTA_DIRECTION_LIGHT flag set.
    ///
    /// * `flags` - Light flags to check.
    pub fn is_delta_light(&self) -> bool {
        self.bits & Self::DELTA_POSITION_LIGHT.bits > 0
            || self.bits & Self::DELTA_DIRECTION_LIGHT.bits > 0
    }
}

impl From<u8> for LightType {
    /// Convert a `u8` value to `LightType`.
    ///
    /// * `t` - A `u8` value containing combination of `*_LIGHT` flags combined
    ///         bitwise OR operator.
    fn from(t: u8) -> Self {
        assert!(
            t <= LightType::INFINITE_LIGHT.bits,
            "Invalid light flags {}=({:#08b})",
            t,
            t
        );
        t.into()
    }
}
