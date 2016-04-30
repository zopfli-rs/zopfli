use libc::{c_ushort, c_uchar, size_t, c_uint};

use symbols::{ZOPFLI_CACHE_LENGTH};

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

/// Returns the length up to which could be stored in the cache.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliMaxCachedSublen(lmc_ptr: *mut ZopfliLongestMatchCache, pos: size_t, _length: size_t) -> c_uint {

    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };

    unsafe {
        let cache = lmc.sublen.offset((ZOPFLI_CACHE_LENGTH * pos * 3) as isize);
        if *cache.offset(1) == 0 && *cache.offset(2) == 0 {
            return 0;  // No sublen cached.
        }
        *cache.offset(((ZOPFLI_CACHE_LENGTH - 1) * 3) as isize) as c_uint + 3
    }
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

    unsafe {
        let cache = lmc.sublen.offset((ZOPFLI_CACHE_LENGTH * pos * 3) as isize);

        for j in 0..ZOPFLI_CACHE_LENGTH {
            let length = *cache.offset((j * 3) as isize) as c_uint + 3;
            let dist = *cache.offset((j * 3 + 1) as isize) as c_ushort + 256 * *cache.offset((j * 3 + 2) as isize) as c_ushort;
            let mut i = prevlength;
            while i <= length {
                *sublen.offset(i as isize) = dist;
                i += 1;
            }
            if length == maxlength {
                break;
            }
            prevlength = length + 1;
        }
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliSublenToCache(sublen: *mut c_ushort, pos: size_t, length: size_t, lmc_ptr: *mut ZopfliLongestMatchCache) {
    let lmc = unsafe {
        assert!(!lmc_ptr.is_null());
        &mut *lmc_ptr
    };

    let mut j: isize = 0;
    let mut bestlength: c_uint = 0;

    if length < 3 {
        return;
    }

    unsafe {

        let cache = lmc.sublen.offset((ZOPFLI_CACHE_LENGTH * pos * 3) as isize);

        let mut i: isize = 3;
        while i <= length as isize {
            if i == length as isize || *sublen.offset(i) != *sublen.offset(i + 1) {
                *cache.offset(j * 3) = (i - 3) as c_uchar;
                *cache.offset(j * 3 + 1) = (*sublen.offset(i)).wrapping_rem(256) as c_uchar;
                *cache.offset(j * 3 + 2) = ((*sublen.offset(i) >> 8)).wrapping_rem(256) as c_uchar;
                bestlength = i as c_uint;
                j += 1;
                if j >= ZOPFLI_CACHE_LENGTH as isize {
                    break;
                }
            }
            i += 1;
        }

        if j < ZOPFLI_CACHE_LENGTH as isize {
            assert!(bestlength == length as c_uint);
            *cache.offset(((ZOPFLI_CACHE_LENGTH - 1) * 3) as isize) = (bestlength - 3) as c_uchar;
        } else {
            assert!(bestlength <= length as c_uint);
        }
        assert!(bestlength == ZopfliMaxCachedSublen(lmc, pos, length));
    }
}
