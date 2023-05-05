//! Common

use crate::{stat_memory_counter, stat_register_fns, stats::*};

stat_memory_counter!("Memory/Primitives", PRIMITIVE_MEMORY, primitive_stats_memory);

stat_register_fns!(primitive_stats_memory);
