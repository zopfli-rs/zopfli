use libc::{c_int, c_ushort, c_uchar, size_t};

use symbols::{ZOPFLI_WINDOW_MASK, ZOPFLI_MIN_MATCH};

const HASH_SHIFT: c_int = 5;
const HASH_MASK: c_int = 32767;

#[repr(C)]
pub struct ZopfliHash {
    head: *mut c_int,  /* Hash value to index of its most recent occurrence. */
    prev: *mut c_ushort,  /* Index to index of prev. occurrence of same hash. */
    hashval: *mut c_int,  /* Index to hash value at this index. */
    val: c_int,  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    head2: *mut c_int,  /* Hash value to index of its most recent occurrence. */
    prev2: *mut c_ushort,  /* Index to index of prev. occurrence of same hash. */
    hashval2: *mut c_int,  /* Index to hash value at this index. */
    val2: c_int,  /* Current hash value. */

    same: *mut c_ushort,  /* Amount of repetitions of same byte after this .*/
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
    let hpos = (pos & ZOPFLI_WINDOW_MASK) as isize;
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
    unsafe {
        *h.hashval.offset(hpos) = h.val;

        if *h.head.offset(h.val as isize) != -1 && *h.hashval.offset(*h.head.offset(h.val as isize) as isize) == h.val {
            *h.prev.offset(hpos) = *h.head.offset(h.val as isize) as c_ushort;
        } else {
            *h.prev.offset(hpos) = hpos as c_ushort;
        }
        *h.head.offset(h.val as isize) = hpos as c_int;

        /* Update "same". */
        if *h.same.offset(((pos - 1) & ZOPFLI_WINDOW_MASK) as isize) > 1 {
            amount = *h.same.offset(((pos - 1) & ZOPFLI_WINDOW_MASK) as isize) as c_int - 1;
        }
        while pos + amount as size_t + 1 < end && *array.offset(pos as isize) == *array.offset((pos + amount as size_t + 1) as isize) && amount < -1 {
            amount += 1;
        }
        *h.same.offset(hpos) = amount as c_ushort;

        h.val2 = (((*h.same.offset(hpos) - ZOPFLI_MIN_MATCH) & 255) ^ h.val as c_ushort) as c_int;
        *h.hashval2.offset(hpos) = h.val2;
        if *h.head2.offset(h.val2 as isize) != -1 as i32 && *h.hashval2.offset(*h.head2.offset(h.val2 as isize) as isize) == h.val2 {
            *h.prev2.offset(hpos) = *h.head2.offset(h.val2 as isize) as c_ushort;
        } else {
            *h.prev2.offset(hpos) = hpos as c_ushort;
        }
        *h.head2.offset(h.val2 as isize) = hpos as c_int;
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHead(h_ptr: *const ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.head
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashPrev(h_ptr: *const ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.prev
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHashval(h_ptr: *const ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.hashval
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashVal(h_ptr: *const ZopfliHash) -> c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.val
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHead2(h_ptr: *const ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.head2
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashPrev2(h_ptr: *const ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.prev2
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashHashval2(h_ptr: *const ZopfliHash) -> *mut c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.hashval2
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashVal2(h_ptr: *const ZopfliHash) -> c_int {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.val2
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashSame(h_ptr: *const ZopfliHash) -> *mut c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &*h_ptr
    };
    h.same
}
