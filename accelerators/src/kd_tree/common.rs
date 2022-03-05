//! KD Tree Common

use core::pbrt::{Axis, Float};
use std::cmp::{Ordering, PartialOrd};

/// Stores information about nodes.
#[derive(Clone)]
pub enum KdAccelNode {
    Interior {
        /// Position of split along a chosen split axis.
        split_pos: Float,

        /// Split axis.
        split_axis: Axis,

        /// The position in the nodes list of the child representing the space
        /// above the splitting plane.
        above_child: usize,
    },

    Leaf {
        /// Number of primitives overlapping the leaf node.
        n_prims: usize,

        /// Stores primitive id when single primitive stored in leaf.
        one_primitive: usize,

        /// Stores offset to the first index for the leaf node when more than one
        /// primitive is stored in leaf.
        primitive_indices_offset: usize,
    },
}

impl KdAccelNode {
    /// Initialize a leaf node.
    ///
    /// * `prim_nums`         - Indices into `primitive_indices`.
    /// * `primitive_indices` - Primitive indices (will be extended with prim_nums
    ///                         if prim_nums.len() > 1).
    pub fn new_leaf(prim_nums: &[usize], primitive_indices: &mut Vec<usize>) -> Self {
        // Store primitive ids for leaf node.
        let np = prim_nums.len();
        if np == 0 {
            Self::Leaf {
                n_prims: np,
                one_primitive: 0,
                primitive_indices_offset: 0,
            }
        } else if np == 1 {
            Self::Leaf {
                n_prims: np,
                one_primitive: prim_nums[0],
                primitive_indices_offset: 0,
            }
        } else {
            let primitive_indices_offset = primitive_indices.len();
            primitive_indices.extend_from_slice(prim_nums);

            Self::Leaf {
                n_prims: np,
                one_primitive: 0,
                primitive_indices_offset,
            }
        }
    }

    /// Initialize an interior node.
    ///
    /// * `axis`  - Axis of split.
    /// * `ac`    - Index of child representing the space above the splitting plane.
    /// * `split` - Split position along `axis`.
    pub fn new_interior(axis: Axis, ac: usize, split: Float) -> Self {
        Self::Interior {
            above_child: ac,
            split_axis: axis,
            split_pos: split,
        }
    }
}

impl Default for KdAccelNode {
    /// Return a default value for `KdAccelNode`.
    fn default() -> Self {
        Self::Leaf {
            n_prims: 0,
            one_primitive: 0,
            primitive_indices_offset: 0,
        }
    }
}

/// EdgeType enumeration.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub(crate) enum EdgeType {
    Start = 0,
    End = 1,
}

/// Represents projections of bounding boxes along an axis.
#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct BoundEdge {
    /// Position of edge along an axis.
    pub(crate) t: Float,

    /// Index of primitive.
    pub(crate) prim_num: usize,

    /// Edge type.
    pub(crate) edge_type: EdgeType,
}

impl BoundEdge {
    /// Create a new `BoundingEdge`.
    ///
    /// * `t`         - Position of edge along an axis.
    /// * `prim_num`  - Index of primitive.
    /// * `edge_type` - Indicates whether edge is starting or ending.
    pub fn new(t: Float, prim_num: usize, edge_type: EdgeType) -> Self {
        Self {
            t,
            prim_num,
            edge_type,
        }
    }
}

impl Default for BoundEdge {
    /// Return a default value for `BoundEdge`.
    fn default() -> Self {
        Self {
            t: 0.0,
            prim_num: 0,
            edge_type: EdgeType::Start,
        }
    }
}

impl PartialOrd for BoundEdge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.t == other.t {
            self.edge_type.partial_cmp(&other.edge_type)
        } else {
            self.t.partial_cmp(&other.t)
        }
    }
}

/// Used to record the nodes yet to be processed for the ray; it is ordered so
/// that the last active entry in the array is the next node that should be
/// considered.
#[derive(Copy, Clone, Default)]
pub(crate) struct KdToDo<'a> {
    /// The current node.
    pub(crate) node: Option<&'a KdAccelNode>,

    /// Index of node.
    pub(crate) node_idx: usize,

    /// Start of parametric range for the ray’s overlap with the current node.
    pub(crate) t_min: Float,

    /// End of parametric range for the ray’s overlap with the current node.
    pub(crate) t_max: Float,
}

/// The maximum number of entries needed in `KdToDo` array is the maximum depth
/// of the kd-tree; the array size used in the following should be more than
/// enough in practice.
pub const MAX_TO_DO: usize = 64;
