//! Bidirectional scattering surface reflectance distribution function.

use std::sync::Arc;

/// BSSRDF trait provides common behavior.
pub trait BSSRDF {}

/// Atomic reference counted `BSSRDF`.
pub type ArcBSSRDF = Arc<dyn BSSRDF + Send + Sync>;
