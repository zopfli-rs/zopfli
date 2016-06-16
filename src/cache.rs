use util::{ZOPFLI_CACHE_LENGTH};

// Cache used by ZopfliFindLongestMatch to remember previously found length/dist
// values.
// This is needed because the squeeze runs will ask these values multiple times for
// the same position.
// Uses large amounts of memory, since it has to remember the distance belonging
// to every possible shorter-than-the-best length (the so called "sublen" array).
pub struct ZopfliLongestMatchCache {
    length: Vec<u16>,
    dist: Vec<u16>,
    sublen: Vec<u8>,
}

impl ZopfliLongestMatchCache {
    pub fn new(blocksize: usize) -> ZopfliLongestMatchCache {
        ZopfliLongestMatchCache {
            /* length > 0 and dist 0 is invalid combination, which indicates on purpose
            that this cache value is not filled in yet. */
            length: vec![1; blocksize],
            dist: vec![0; blocksize],
            /* Rather large amount of memory. */
            sublen: vec![0; ZOPFLI_CACHE_LENGTH * blocksize * 3],
        }
    }

    pub fn length_at(&self, pos: usize) -> u16 {
        self.length[pos]
    }

    pub fn dist_at(&self, pos: usize) -> u16 {
        self.dist[pos]
    }

    pub fn store_length_at(&mut self, pos: usize, val: u16) {
        self.length[pos] = val;
    }

    pub fn store_dist_at(&mut self, pos: usize, val: u16) {
        self.dist[pos] = val;
    }

    /// Returns the length up to which could be stored in the cache.
    pub fn max_sublen(&self, pos: usize) -> u32 {
        let start = ZOPFLI_CACHE_LENGTH * pos * 3;
        if self.sublen[start + 1] == 0 && self.sublen[start + 2] == 0 {
            return 0;  // No sublen cached.
        }
        self.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] as u32 + 3
    }

    /// Stores sublen array in the cache.
    pub fn store_sublen(&mut self, sublen: &[u16], pos: usize, length: usize) {
        if length < 3 {
            return;
        }

        let start = ZOPFLI_CACHE_LENGTH * pos * 3;
        let mut i = 3;
        let mut j = 0;
        let mut bestlength = 0;
        while i <= length {
            if i == length || sublen[i] != sublen[i + 1] {
                self.sublen[start + (j * 3)] = (i - 3) as u8;
                self.sublen[start + (j * 3 + 1)] = sublen[i].wrapping_rem(256) as u8;
                self.sublen[start + (j * 3 + 2)] = ((sublen[i] >> 8)).wrapping_rem(256) as u8;
                bestlength = i as u32;
                j += 1;
                if j >= ZOPFLI_CACHE_LENGTH {
                    break;
                }
            }
            i += 1;
        }

        if j < ZOPFLI_CACHE_LENGTH {
            assert!(bestlength == length as u32);
            self.sublen[start + ((ZOPFLI_CACHE_LENGTH - 1) * 3)] = (bestlength - 3) as u8;
        } else {
            assert!(bestlength <= length as u32);
        }
        assert!(bestlength == self.max_sublen(pos));
    }

    /// Extracts sublen array from the cache.
    pub fn fetch_sublen(&self, pos: usize, length: usize) -> Option<Vec<u16>> {
        if length < 3 {
            return None;
        }

        let start = ZOPFLI_CACHE_LENGTH * pos * 3;
        let maxlength = self.max_sublen(pos) as usize;
        let mut sublen = vec![0; maxlength + 1];
        let mut prevlength = 0;

        for j in 0..ZOPFLI_CACHE_LENGTH {
            let length = self.sublen[start + (j * 3)] as usize + 3;
            let dist = self.sublen[start + (j * 3 + 1)] as u16 + 256 * self.sublen[start + (j * 3 + 2)] as u16;

            let mut i = prevlength;
            while i <= length {
                sublen[i] = dist;
                i += 1;
            }
            if length == maxlength {
                break;
            }
            prevlength = length + 1;
        }
        Some(sublen)
    }
}
