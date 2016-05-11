use libc::{size_t, c_int};

// Bounded package merge algorithm, based on the paper
// "A Fast and Space-Economical Algorithm for Length-Limited Coding
// Jyrki Katajainen, Alistair Moffat, Andrew Turpin".

/// Nodes forming chains. Also used to represent leaves.
#[repr(C)]
pub struct Node {
  weight: size_t,     // Total weight (symbol count) of this chain.
  tail: *const Node,  // Previous node(s) of this chain, or 0 if none.
  count: c_int,       // Leaf symbol index, or number of leaves before this chain.
}

/// Memory pool for nodes.
#[repr(C)]
pub struct NodePool {
  next: *const Node,   // Pointer to a possibly free node in the pool.
}


/// Initializes a chain node with the given values and marks it as in use.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn InitNode(weight: size_t, count: c_int, tail: *const Node, node_ptr: *mut Node) {
    let node = unsafe {
        assert!(!node_ptr.is_null());
        &mut *node_ptr
    };

    node.weight = weight;
    node.count = count;
    node.tail = tail;
}
