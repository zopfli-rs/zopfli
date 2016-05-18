use libc::{size_t, c_int, c_uint};

// Bounded package merge algorithm, based on the paper
// "A Fast and Space-Economical Algorithm for Length-Limited Coding
// Jyrki Katajainen, Alistair Moffat, Andrew Turpin".

/// Nodes forming chains.
#[repr(C)]
pub struct Node {
  weight: size_t,     // Total weight (symbol count) of this chain.
  tail: *const Node,  // Previous node(s) of this chain, or 0 if none.
  count: c_int,       // Leaf symbol index, or number of leaves before this chain.
}

#[repr(C)]
pub struct Leaf {
  weight: size_t,     // Total weight (symbol count) of this chain.
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

/// Converts result of boundary package-merge to the bitlengths. The result in the
/// last chain of the last list contains the amount of active leaves in each list.
/// chain: Chain to extract the bit length from (last chain from last list).
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ExtractBitLengths(chain: *const Node, leaves: *const Leaf, bitlengths: *mut c_uint) {
    let mut counts = [0; 16];
    let mut end = 16;
    let mut ptr = 15;
    let mut value = 1;

    let mut node_ptr = chain;
    while !node_ptr.is_null() {
        let node = unsafe {
           &*node_ptr
        };

        end -= 1;
        counts[end] = node.count;

        node_ptr = node.tail;
    }

    let mut val = counts[15];
    while ptr >= end {
        while val > counts[ptr - 1] {
            unsafe {
                let leaf = &*leaves.offset((val - 1) as isize);
                *bitlengths.offset(leaf.count as isize) = value;
            }
            val -= 1;
        }
        ptr -= 1;
        value += 1;
    }
}
