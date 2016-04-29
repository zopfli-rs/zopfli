use libc::{c_int, c_ushort, c_uchar};

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
    val2: *mut c_int,  /* Current hash value. */

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
