//! KD Tree Common

#![allow(dead_code)]
use crate::core::pbrt::Float;

/// Information about leaf and interior nodes.
#[derive(Copy, Clone)]
pub struct KdAccelNode {
    /// Split data.
    split: SplitData,

    /// Node data.
    node: NodeData,
}

impl KdAccelNode {
    /// Initialize a leaf node.
    ///
    /// * `prim_nums`         - Indices into `primitive_indices`.
    /// * `np`                - Number of primitives overlapping the leaf node.
    /// * `primitive_indices` - Primitive indices.
    pub fn init_leaf(&mut self, prim_nums: &[u32], np: u32, primitive_indices: &mut Vec<u32>) {
        // First set lower 2 bits to `3` for flag to indicate leaf node.
        self.node.flags = 3;

        // Set the higher 30 bits for n_prims by first shifting `np` by 2 bits.
        unsafe {
            self.node.n_prims |= np << 2;
        }

        // Store primitive ids for leaf node.
        if np == 0 {
            self.split.one_primitive = 0;
        } else if np == 1 {
            self.split.one_primitive = prim_nums[0];
        } else {
            self.split.primitive_indices_offset = primitive_indices.len() as u32;
            for i in 0..(np as usize) {
                primitive_indices.push(prim_nums[i]);
            }
        }
    }

    /// Initialize an interior node.
    ///
    /// * `axis`  - Axis of split (0 = x, 1 = y, 2 = z).
    /// * `ac`    - Index of child representing the space above the splitting plane.
    /// * `split` - Split position along `axis`.
    pub fn init_interior(&mut self, axis: u32, ac: u32, split: Float) {
        // Set split position.
        self.split.pos = split;

        // First set lower 2 bits to `axis` for flag to indicate split axis.
        self.node.flags = axis;

        // Set the higher 30 bits for `above_child` by first shifting `ac` by 2 bits.
        unsafe {
            self.node.above_child |= ac << 2;
        }
    }

    /// Returns the position along the split axis.
    pub fn split_pos(&self) -> Float {
        let pos: Float;
        unsafe {
            pos = self.split.pos;
        }
        pos
    }

    /// Returns the split axis.
    pub fn split_axis(&self) -> u32 {
        let axis: u32;
        unsafe {
            axis = self.node.flags & 3;
        }
        axis
    }

    /// Returns the number of primitives overlapping the leaf node.
    pub fn n_primitives(&self) -> u32 {
        let n_prims: u32;
        unsafe {
            n_prims = self.node.n_prims;
        }
        n_prims >> 2
    }

    /// Returns whether the node is a leaf node.
    pub fn is_leaf(&self) -> bool {
        let flags: u32;
        unsafe {
            flags = self.node.flags;
        }
        flags & 3 == 3
    }

    /// For an interior node, returns the position, in the nodes list, of the
    /// child representing the space above the splitting plane.
    pub fn above_child(&self) -> u32 {
        let ac: u32;
        unsafe {
            ac = self.node.above_child;
        }
        ac >> 2
    }

    /// For leaf nodes, returns primitive ids or 0 if there are none.
    pub fn one_primitive(&self) -> u32 {
        let op: u32;
        unsafe { op = self.split.one_primitive }
        op
    }

    /// For leaf nodes, returns the offset to the first index for the leaf node.
    pub fn primitive_indices_offset(&self) -> u32 {
        let offset: u32;
        unsafe { offset = self.split.primitive_indices_offset }
        offset
    }
}

impl Default for KdAccelNode {
    /// Return a default value for `KdAccelNode`.
    fn default() -> Self {
        Self {
            split: SplitData {
                one_primitive: 0_u32,
            },
            node: NodeData { flags: 0_u32 },
        }
    }
}

/// Efficiently stores information about splits in 32 bits.
#[derive(Copy, Clone)]
#[repr(C)]
union SplitData {
    /// For interior nodes, determines the position of split along a chosen axis
    /// (`see NodeData::flags`).
    pub pos: Float,

    /// For leaf nodes, stores primitive ids or 0 if there are none.
    pub one_primitive: u32,

    /// For leaf nodes, offset to the first index for the leaf node.
    pub primitive_indices_offset: u32,
}

/// Efficiently stores information about nodes in 32 bits.
#[derive(Copy, Clone)]
#[repr(C)]
union NodeData {
    /// For interior nodes, flags is used to differentiate the x, y andz split
    /// axis with values 0, 1, and 2 respectively.  For leaf nodes, the value is
    /// set to 3. The lower 2 bits are used.
    pub flags: u32,

    /// Number of primitives overlapping the leaf node. The upper 30 bits are
    /// used.
    pub n_prims: u32,

    /// For interior nodes, this stores the position, in the nodes list, of the
    /// child representing the space above the splitting plane. The upper 30
    /// bits are used.
    pub above_child: u32,
}

/// EdgeType enumeration.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub enum EdgeType {
    Start,
    End,
}

/// Represents projections of bounding boxes along an axis.
#[derive(Copy, Clone)]
pub struct BoundEdge {
    pub t: Float,
    pub prim_num: u32,
    pub edge_type: EdgeType,
}

impl BoundEdge {
    /// Create a new `BoundingEdge`.
    ///
    /// * `t`           -
    /// * `prim_num`    - * `is_starting` -
    pub fn new(t: Float, prim_num: u32, is_starting: bool) -> Self {
        Self {
            t,
            prim_num,
            edge_type: if is_starting {
                EdgeType::Start
            } else {
                EdgeType::End
            },
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

/// KdToDo is used to record the nodes yet to be processed for the ray; it is
/// ordered so that the last active entry in the array is the next node that
/// should be considered.
#[derive(Copy, Clone, Default)]
pub struct KdToDo<'a> {
    pub node: Option<&'a KdAccelNode>,
    pub node_idx: usize,
    pub t_min: Float,
    pub t_max: Float,
}

/// The maximum number of entries needed in KdToDo array is the maximum depth of
/// the kd-tree; the array size used in the following should be more than enough
/// in practice.
pub const MAX_TO_DO: usize = 64;
