//! Reflection and surface scattering models.

use std::sync::Arc;

/// BSDF model.
pub struct BSDF {}

/// Atomic reference counted `BSDF`.
pub type ArcBSDF = Arc<BSDF>;
