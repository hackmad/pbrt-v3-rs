//! Parameter Set Items

#![allow(dead_code)]
use std::fmt;

/// Stores a parameter set item consisting of a list of values of type `T`.
#[derive(Clone, Default)]
pub struct ParamSetItem<T: fmt::Display> {
    /// The values.
    pub values: Vec<T>,
}

impl<T: fmt::Display> ParamSetItem<T> {
    /// Create new `ParamSet<T>`.
    /// * `values`    - The values.
    pub fn new(values: Vec<T>) -> Self {
        Self { values }
    }
}

impl<T: fmt::Display> fmt::Display for ParamSetItem<T> {
    /// Formats the value using the given formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for (i, v) in self.values.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:}", v)?;
        }
        write!(f, "]")
    }
}
