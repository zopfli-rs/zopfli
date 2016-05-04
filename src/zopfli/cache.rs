use libc::{c_ushort, c_uchar, size_t, c_uint};

use symbols::{ZOPFLI_CACHE_LENGTH};

// Cache used by ZopfliFindLongestMatch to remember previously found length/dist
// values.
// This is needed because the squeeze runs will ask these values multiple times for
// the same position.
// Uses large amounts of memory, since it has to remember the distance belonging
// to every possible shorter-than-the-best length (the so called "sublen" array).
pub struct ZopfliLongestMatchCache {
    length: Vec<c_ushort>,
    dist: Vec<c_ushort>,
    sublen: Vec<c_uchar>,
}

impl ZopfliLongestMatchCache {
    pub fn new(blocksize: size_t) -> ZopfliLongestMatchCache {
        ZopfliLongestMatchCache {
            /* length > 0 and dist 0 is invalid combination, which indicates on purpose
            that this cache value is not filled in yet. */
            length: vec![1; blocksize],
            dist: vec![0; blocksize],
            /* Rather large amount of memory. */
            sublen: vec![0; ZOPFLI_CACHE_LENGTH * blocksize * 3],
        }
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliInitCache(blocksize: size_t) -> *mut ZopfliLongestMatchCache {
    Box::into_raw(Box::new(ZopfliLongestMatchCache::new(blocksize)))
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCleanCache(ptr: *mut ZopfliLongestMatchCache) {
    if ptr.is_null() { return }
    unsafe { Box::from_raw(ptr); }
}

/// Returns the length up to which could be stored in the cache.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliMaxCachedSublen(lmc_ptr: *mut ZopfliLongestMatchCache, pos: size_t, _length: size_t) -> c_uint {

    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };

    let start = ZOPFLI_CACHE_LENGTH * pos * 3;
    if lmc.sublen[start + 1] == 0 && lmc.sublen[start + 2] == 0 {
        return 0;  // No sublen cached.
    }
    lmc.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] as c_uint + 3
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCacheToSublen(lmc_ptr: *mut ZopfliLongestMatchCache, pos: size_t, length: size_t, sublen: *mut c_ushort) {
    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };

    let maxlength = ZopfliMaxCachedSublen(lmc_ptr, pos, length);
    let mut prevlength = 0;

    if length < 3 {
        return;
    }

    let start = ZOPFLI_CACHE_LENGTH * pos * 3;

    for j in 0..ZOPFLI_CACHE_LENGTH {
        let length = lmc.sublen[start + (j * 3)] as c_uint + 3;
        let dist = lmc.sublen[start + (j * 3 + 1)] as c_ushort + 256 * lmc.sublen[start + (j * 3 + 2)] as c_ushort;

        let mut i = prevlength;
        while i <= length {
            unsafe {
                *sublen.offset(i as isize) = dist;
            }
            i += 1;
        }
        if length == maxlength {
            break;
        }
        prevlength = length + 1;
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliSublenToCache(sublen: *mut c_ushort, pos: size_t, length: size_t, lmc_ptr: *mut ZopfliLongestMatchCache) {
    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };

    let mut j: usize = 0;
    let mut bestlength: c_uint = 0;

    if length < 3 {
        return;
    }

    let start = ZOPFLI_CACHE_LENGTH * pos * 3;

    let mut i: isize = 3;
    while i <= length as isize {
        unsafe {
            if i == length as isize || *sublen.offset(i) != *sublen.offset(i + 1) {
                lmc.sublen[start + (j * 3)] = (i - 3) as c_uchar;
                lmc.sublen[start + (j * 3 + 1)] = (*sublen.offset(i)).wrapping_rem(256) as c_uchar;
                lmc.sublen[start + (j * 3 + 2)] = ((*sublen.offset(i) >> 8)).wrapping_rem(256) as c_uchar;
                bestlength = i as c_uint;
                j += 1;
                if j >= ZOPFLI_CACHE_LENGTH {
                    break;
                }
            }
        }
        i += 1;
    }

    if j < ZOPFLI_CACHE_LENGTH {
        assert!(bestlength == length as c_uint);
        lmc.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] = (bestlength - 3) as c_uchar;
    } else {
        assert!(bestlength <= length as c_uint);
    }
    assert!(bestlength == ZopfliMaxCachedSublen(lmc, pos, length));
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCacheLengthAt(lmc_ptr: *mut ZopfliLongestMatchCache, pos: size_t) -> c_ushort {
    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };
    lmc.length[pos]
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCacheDistAt(lmc_ptr: *mut ZopfliLongestMatchCache, pos: size_t) -> c_ushort {
    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };
    lmc.dist[pos]
}
