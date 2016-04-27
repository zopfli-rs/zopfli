use libc::{c_ushort, c_uchar};

// Cache used by ZopfliFindLongestMatch to remember previously found length/dist
// values.
// This is needed because the squeeze runs will ask these values multiple times for
// the same position.
// Uses large amounts of memory, since it has to remember the distance belonging
// to every possible shorter-than-the-best length (the so called "sublen" array).
#[repr(C)]
pub struct ZopfliLongestMatchCache {
    length: *mut c_ushort,
    dist: *mut c_ushort,
    sublen: *mut c_uchar,
}
