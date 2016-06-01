use libc::{c_ushort, c_uchar, size_t, c_uint};

use util::{ZOPFLI_CACHE_LENGTH};

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

    pub fn length_at(&self, pos: size_t) -> c_ushort {
        self.length[pos]
    }

    pub fn dist_at(&self, pos: size_t) -> c_ushort {
        self.dist[pos]
    }

    /// Returns the length up to which could be stored in the cache.
    pub fn max_sublen(&self, pos: size_t) -> c_uint {
        let start = ZOPFLI_CACHE_LENGTH * pos * 3;
        if self.sublen[start + 1] == 0 && self.sublen[start + 2] == 0 {
            return 0;  // No sublen cached.
        }
        self.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] as c_uint + 3
    }

    /// Stores sublen array in the cache.
    pub fn store_sublen(&mut self, sublen: *mut c_ushort, pos: size_t, length: size_t) {
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
                    self.sublen[start + (j * 3)] = (i - 3) as c_uchar;
                    self.sublen[start + (j * 3 + 1)] = (*sublen.offset(i)).wrapping_rem(256) as c_uchar;
                    self.sublen[start + (j * 3 + 2)] = ((*sublen.offset(i) >> 8)).wrapping_rem(256) as c_uchar;
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
            self.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] = (bestlength - 3) as c_uchar;
        } else {
            assert!(bestlength <= length as c_uint);
        }
        assert!(bestlength == self.max_sublen(pos));
    }

    /// Extracts sublen array from the cache.
    pub fn fetch_sublen(&self, pos: size_t, length: size_t, sublen: *mut c_ushort) {
        let maxlength = self.max_sublen(pos);
        let mut prevlength = 0;

        if length < 3 {
            return;
        }

        let start = ZOPFLI_CACHE_LENGTH * pos * 3;

        for j in 0..ZOPFLI_CACHE_LENGTH {
            let length = self.sublen[start + (j * 3)] as c_uint + 3;
            let dist = self.sublen[start + (j * 3 + 1)] as c_ushort + 256 * self.sublen[start + (j * 3 + 2)] as c_ushort;

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
}
