use libc::{c_int, c_ushort, c_uchar, size_t};

use symbols::{ZOPFLI_WINDOW_MASK, ZOPFLI_MIN_MATCH};

const HASH_SHIFT: c_int = 5;
const HASH_MASK: c_int = 32767;

#[repr(C)]
pub struct ZopfliHash {
    head: Vec<c_int>,  /* Hash value to index of its most recent occurrence. */
    prev: Vec<c_ushort>,  /* Index to index of prev. occurrence of same hash. */
    hashval: Vec<c_int>,  /* Index to hash value at this index. */
    val: c_int,  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    head2: Vec<c_int>,  /* Hash value to index of its most recent occurrence. */
    prev2: Vec<c_ushort>,  /* Index to index of prev. occurrence of same hash. */
    hashval2: Vec<c_int>,  /* Index to hash value at this index. */
    val2: c_int,  /* Current hash value. */

    same: Vec<c_ushort>,  /* Amount of repetitions of same byte after this .*/
}

impl ZopfliHash {
    pub fn new(window_size: size_t) -> ZopfliHash {
        ZopfliHash {
           head: vec![-1; 65536],
           prev: (0..window_size as c_ushort).collect::<Vec<_>>(),
           hashval: vec![-1; window_size],
           val: 0,

           /* Fields with similar purpose as the above hash, but for the second hash with
           a value that is calculated differently.  */
           head2: vec![-1; 65536],
           prev2: (0..window_size as c_ushort).collect::<Vec<_>>(),
           hashval2: vec![-1; window_size],
           val2: 0,

           same: vec![0; window_size],
       }
    }
}

/// Update the sliding hash value with the given byte. All calls to this function
/// must be made on consecutive input characters. Since the hash value exists out
/// of multiple input bytes, a few warmups with this function are needed initially.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn UpdateHashValue(h_ptr: *mut ZopfliHash, c: c_uchar) {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.val = ((h.val << HASH_SHIFT) ^ c as c_int) & HASH_MASK;
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliWarmupHash(array: *const c_uchar, pos: size_t, end: size_t, h_ptr: *mut ZopfliHash) {
    unsafe {
        UpdateHashValue(h_ptr, *array.offset((pos + 0) as isize));
    }
    if pos + 1 < end {
        unsafe {
            UpdateHashValue(h_ptr, *array.offset((pos + 1) as isize));
        }
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliUpdateHash(array: *const c_uchar, pos: size_t, end: size_t, h_ptr: *mut ZopfliHash) {
    let hpos = (pos & ZOPFLI_WINDOW_MASK) as usize;
    let mut amount: c_int = 0;

    let hash_value = if pos + ZOPFLI_MIN_MATCH as size_t <= end {
        unsafe { *array.offset((pos + ZOPFLI_MIN_MATCH as size_t - 1) as isize) }
    } else {
        0
    };
    UpdateHashValue(h_ptr, hash_value);

    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.hashval[hpos] = h.val;

    if h.head[h.val as usize] != -1 && h.hashval[h.head[h.val as usize] as usize] == h.val {
        h.prev[hpos] = h.head[h.val as usize] as c_ushort;
    } else {
        h.prev[hpos] = hpos as c_ushort;
    }

    h.head[h.val as usize] = hpos as c_int;

    /* Update "same". */
    if h.same[((pos - 1) & ZOPFLI_WINDOW_MASK) as usize] > 1 {
        amount = h.same[((pos - 1) & ZOPFLI_WINDOW_MASK) as usize] as c_int - 1;
    }

    unsafe {
        while pos + amount as size_t + 1 < end && *array.offset(pos as isize) == *array.offset((pos + amount as size_t + 1) as isize) && amount < -1 {
            amount += 1;
        }
    }
    h.same[hpos] = amount as c_ushort;

    h.val2 = (((h.same[hpos] - ZOPFLI_MIN_MATCH) & 255) ^ h.val as c_ushort) as c_int;
    h.hashval2[hpos] = h.val2;
    if h.head2[h.val2 as usize] != -1 as i32 && h.hashval2[h.head2[h.val2 as usize] as usize] == h.val2 {
        h.prev2[hpos] = h.head2[h.val2 as usize] as c_ushort;
    } else {
        h.prev2[hpos] = hpos as c_ushort;
    }
    h.head2[h.val2 as usize] = hpos as c_int;
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHead(h_ptr: *mut ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.head.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashPrev(h_ptr: *mut ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.prev.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHashval(h_ptr: *mut ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.hashval.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashVal(h_ptr: *mut ZopfliHash) -> c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.val
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHead2(h_ptr: *mut ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.head2.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashPrev2(h_ptr: *mut ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.prev2.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHashval2(h_ptr: *mut ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.hashval2.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashVal2(h_ptr: *mut ZopfliHash) -> c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.val2
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashSame(h_ptr: *mut ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.same.as_mut_ptr()
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCleanHash(ptr: *mut ZopfliHash) {
    if ptr.is_null() { return }
    unsafe { Box::from_raw(ptr); }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliInitHash(window_size: size_t) -> *mut ZopfliHash {
    Box::into_raw(Box::new(ZopfliHash::new(window_size)))
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliResetHash(window_size: size_t, h_ptr: *mut ZopfliHash) {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };

    h.val = 0;
    h.head = vec![-1; 65536];
    h.prev = (0..window_size as c_ushort).collect::<Vec<_>>();
    h.hashval = vec![-1; window_size];

    h.same = vec![0; window_size];

    h.val2 = 0;
    h.head2 = vec![-1; 65536];
    h.prev2 = (0..window_size as c_ushort).collect::<Vec<_>>();
    h.hashval2 = vec![-1; window_size];
}
