//! File Utility Functions

use std::path::PathBuf;
use std::result::Result;

/// Returns the absolute path after resolving the given path.
///
/// * `path` - The path.
pub fn absolute_path(path: &str) -> Result<String, String> {
    match PathBuf::from(path)
        .canonicalize()
        .map(PathBuf::into_os_string)
        .map(|s| s.into_string().ok())
    {
        Ok(Some(abs_path)) => Ok(abs_path),
        Ok(None) => Err(format!("invalid path {}", path)),
        Err(err) => Err(format!("invalid path {}. {}.", path, err)),
    }
}
