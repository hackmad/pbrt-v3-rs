//! File Utility Functions

#![allow(dead_code)]

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

/// Returns the path to the parent folder; or `None` if path is root.
///
/// * `path` - The path.
pub fn parent_path(path: &str) -> Option<String> {
    PathBuf::from(path)
        .parent()
        .map_or(None, |p| p.to_str().map(|s| String::from(s)))
}

/// Returns `true` if given path is a relative path; otherwise returns `false`.
///
/// * `path` - The path.
pub fn is_relative_path(path: &str) -> bool {
    PathBuf::from(path).is_relative()
}

/// Returns `true` if given path is a absolute path; otherwise returns `false`.
///
/// * `path` - The path.
pub fn is_absolute_path(path: &str) -> bool {
    PathBuf::from(path).is_absolute()
}
