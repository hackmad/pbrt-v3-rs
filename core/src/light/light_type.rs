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
        self.contains(Self::DELTA_POSITION_LIGHT) || self.contains(Self::DELTA_DIRECTION_LIGHT)
    }
}
