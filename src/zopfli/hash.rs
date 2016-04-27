use libc::{c_int, c_ushort};

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
