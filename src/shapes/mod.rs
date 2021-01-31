//! Geometry

mod cone;
mod curve;
mod cylinder;
mod disk;
mod hyperboloid;
mod loopsubdiv;
mod paraboloid;
mod sphere;
mod triangle;

// Re-export
pub use cone::*;
pub use curve::*;
pub use cylinder::*;
pub use disk::*;
pub use hyperboloid::*;
pub use loopsubdiv::*;
pub use paraboloid::*;
pub use sphere::*;
pub use triangle::*;

use crate::core::geometry::ArcTransform;
use crate::core::paramset::ParamSet;

/// Stores properties for shape creation.
#[derive(Clone)]
pub struct ShapeProps {
    /// Parameter set.
    pub params: ParamSet,

    /// Transformation from object space to world space.
    pub o2w: ArcTransform,

    /// Transformation from world space to object space.
    pub w2o: ArcTransform,

    /// Indicates whether their surface normal directions should be reversed
    /// from the default.
    pub reverse_orientation: bool,
}
