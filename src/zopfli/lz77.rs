use std::{slice, ptr};

use libc::{size_t, c_ushort, c_uchar, c_int, c_uint};

use cache::{ZopfliLongestMatchCache};
use symbols::{ZopfliGetLengthSymbol, ZopfliGetDistSymbol, ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_MAX_MATCH, ZOPFLI_MIN_MATCH};
use zopfli::ZopfliOptions;

// Comment from C:
// Stores lit/length and dist pairs for LZ77.
// Parameter litlens: Contains the literal symbols or length values.
// Parameter dists: Contains the distances. A value is 0 to indicate that there is
// no dist and the corresponding litlens value is a literal instead of a length.
// Parameter size: The size of both the litlens and dists arrays.
// The memory can best be managed by using ZopfliInitLZ77Store to initialize it,
// ZopfliCleanLZ77Store to destroy it, and ZopfliStoreLitLenDist to append values.

#[repr(C)]
pub struct ZopfliLZ77Store {
  litlens: *mut c_ushort,  /* Lit or len. */
  dists: *mut c_ushort,  /* If 0: indicates literal in corresponding litlens,
      if > 0: length in corresponding litlens, this is the distance. */
  size: size_t,

  data: *mut c_uchar,  /* original data */
  pos: *mut size_t,  /* position in data where this LZ77 command begins */

  ll_symbol: *mut c_ushort,
  d_symbol: *mut c_ushort,

  /* Cumulative histograms wrapping around per chunk. Each chunk has the amount
  of distinct symbols as length, so using 1 value per LZ77 symbol, we have a
  precise histogram at every N symbols, and the rest can be calculated by
  looping through the actual symbols of this chunk. */
  ll_counts: *mut size_t,
  d_counts: *mut size_t,
}

#[derive(Debug)]
#[repr(C)]
pub struct Lz77Store {
   pub litlens: Vec<c_ushort>,
   pub dists: Vec<c_ushort>,

   pos: Vec<size_t>,

   ll_symbol: Vec<c_ushort>,
   d_symbol: Vec<c_ushort>,

   ll_counts: Vec<size_t>,
   d_counts: Vec<size_t>,
}

impl Lz77Store {
    pub fn new() -> Lz77Store {
        Lz77Store {
          litlens: vec![],
          dists: vec![],

          pos: vec![],

          ll_symbol: vec![],
          d_symbol: vec![],

          ll_counts: vec![],
          d_counts: vec![],
       }
    }

    pub fn lit_len_dist(&mut self, length: c_ushort, dist: c_ushort, pos: size_t) {
        let origsize = self.litlens.len();
        let llstart = ZOPFLI_NUM_LL * (origsize / ZOPFLI_NUM_LL);
        let dstart = ZOPFLI_NUM_D * (origsize / ZOPFLI_NUM_D);

        if origsize % ZOPFLI_NUM_LL == 0 {
            for i in 0..ZOPFLI_NUM_LL {
                if origsize == 0 {
                    self.ll_counts.push(0);
                } else {
                    let last_histogram_value = self.ll_counts[origsize - ZOPFLI_NUM_LL + i];
                    self.ll_counts.push(last_histogram_value);
                }
            }
        }

        if origsize % ZOPFLI_NUM_D == 0 {
            for i in 0..ZOPFLI_NUM_D {
                if origsize == 0 {
                    self.d_counts.push(0);
                } else {
                    let last_histogram_value = self.d_counts[origsize - ZOPFLI_NUM_D + i];
                    self.d_counts.push(last_histogram_value);
                }
            }
        }

        self.litlens.push(length);
        self.dists.push(dist);
        self.pos.push(pos);

        // Why isn't this at the beginning of this function?
        // assert(length < 259);

        if dist == 0 {
            self.ll_symbol.push(length);
            self.d_symbol.push(0);
            self.ll_counts[llstart + length as usize] += 1;
        } else {
            self.ll_symbol.push(ZopfliGetLengthSymbol(length as c_int) as c_ushort);
            self.d_symbol.push(ZopfliGetDistSymbol(dist as c_int) as c_ushort);
            self.ll_counts[llstart + ZopfliGetLengthSymbol(length as c_int) as usize] += 1;
            self.d_counts[dstart + ZopfliGetDistSymbol(dist as c_int) as usize] += 1;
        }
    }
}

#[no_mangle]
pub extern fn lz77_store_new() -> *mut Lz77Store {
    Box::into_raw(Box::new(Lz77Store::new()))
}

#[no_mangle]
pub extern fn lz77_store_free(ptr: *mut Lz77Store) {
    if ptr.is_null() { return }
    unsafe { Box::from_raw(ptr); }
}

#[no_mangle]
pub extern fn lz77_store_lit_len_dist(ptr: *mut Lz77Store, length: c_ushort, dist: c_ushort, pos: size_t) {
    let store = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    store.lit_len_dist(length, dist, pos);
}

#[no_mangle]
pub extern fn lz77_store_from_c(store: *const ZopfliLZ77Store) -> *mut Lz77Store {
    Box::into_raw(Box::new(store.into()))
}

impl From<*const ZopfliLZ77Store> for Lz77Store {
    fn from(ptr: *const ZopfliLZ77Store) -> Lz77Store {
        let store = unsafe {
            assert!(!ptr.is_null());
            &*ptr
        };

        let len = store.size as usize;
        let ll_len = (ZOPFLI_NUM_LL * (store.size / ZOPFLI_NUM_LL) + ZOPFLI_NUM_LL) as usize;
        let d_len = (ZOPFLI_NUM_D * (store.size / ZOPFLI_NUM_D) + ZOPFLI_NUM_D) as usize;

        Lz77Store {
            litlens: ptr_to_vec(store.litlens, len),
            dists: ptr_to_vec(store.dists, len),

            pos: ptr_to_vec(store.pos, len),

            ll_symbol: ptr_to_vec(store.ll_symbol, len),
            d_symbol: ptr_to_vec(store.d_symbol, len),

            ll_counts: ptr_to_vec(store.ll_counts, ll_len),
            d_counts: ptr_to_vec(store.d_counts, d_len),
        }
    }
}

fn ptr_to_vec<T: Clone>(ptr: *mut T, length: usize) -> Vec<T> {
    if ptr.is_null() {
        vec![]
    } else {
        unsafe { slice::from_raw_parts(ptr, length).to_vec() }
    }
}

#[no_mangle]
pub extern fn lz77_store_result(ptr: *mut Lz77Store, store: &mut ZopfliLZ77Store) {
    let lz77 = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let len = lz77.litlens.len();

    store.litlens = lz77.litlens.as_mut_ptr();
    store.dists = lz77.dists.as_mut_ptr();
    store.size = len;
    store.pos = lz77.pos.as_mut_ptr();
    store.ll_symbol = lz77.ll_symbol.as_mut_ptr();
    store.d_symbol = lz77.d_symbol.as_mut_ptr();
    store.ll_counts = lz77.ll_counts.as_mut_ptr();
    store.d_counts = lz77.d_counts.as_mut_ptr();
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliInitLZ77Store(data: *mut c_uchar, store_ptr: *mut ZopfliLZ77Store) {
    let store = unsafe {
        assert!(!store_ptr.is_null());
        &mut *store_ptr
    };

    store.size = 0;
    store.litlens = ptr::null_mut();
    store.dists = ptr::null_mut();
    store.pos = ptr::null_mut();
    store.data = data;
    store.ll_symbol = ptr::null_mut();
    store.d_symbol = ptr::null_mut();
    store.ll_counts = ptr::null_mut();
    store.d_counts = ptr::null_mut();
}

/// Some state information for compressing a block.
/// This is currently a bit under-used (with mainly only the longest match cache),
/// but is kept for easy future expansion.
#[repr(C)]
pub struct ZopfliBlockState {
    options: *const ZopfliOptions,
    /* Cache for length/distance pairs found so far. */
    lmc: *mut ZopfliLongestMatchCache,
    /* The start (inclusive) and end (not inclusive) of the current block. */
    blockstart: size_t,
    blockend: size_t,
}

#[repr(C)]
pub struct LongestMatch {
    distance: c_ushort,
    length: c_ushort,
    from_cache: c_int,
    limit: size_t,
}

/// Gets distance, length and sublen values from the cache if possible.
/// Returns 1 if it got the values from the cache, 0 if not.
/// Updates the limit value to a smaller one if possible with more limited
/// information from the cache.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn TryGetFromLongestMatchCache(s_ptr: *mut ZopfliBlockState, pos: size_t, mut limit: size_t, sublen: *mut c_ushort) -> LongestMatch {
    let s = unsafe {
        assert!(!s_ptr.is_null());
        &mut *s_ptr
    };

    let mut longest_match = LongestMatch {
        distance: 0,
        length: 0,
        from_cache: 0,
        limit: limit,
    };

    if s.lmc.is_null() {
        return longest_match;
    }

    /* The LMC cache starts at the beginning of the block rather than the
     beginning of the whole array. */
    let lmcpos = pos - s.blockstart;

    /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
      that this cache value is not filled in yet. */
    let length_lmcpos = unsafe { (*s.lmc).length_at(lmcpos) };
    let dist_lmcpos = unsafe { (*s.lmc).dist_at(lmcpos) };
    let cache_available = length_lmcpos == 0 || dist_lmcpos != 0;
    let max_sublen = unsafe { (*s.lmc).max_sublen(lmcpos) };
    let limit_ok_for_cache = cache_available &&
       (limit == ZOPFLI_MAX_MATCH || length_lmcpos <= limit as c_ushort ||
       (!sublen.is_null() && max_sublen >= limit as c_uint));

    if limit_ok_for_cache && cache_available {
        if sublen.is_null() || length_lmcpos as c_uint <= max_sublen {
            let mut length = length_lmcpos;
            if length > limit as c_ushort {
                length = limit as c_ushort;
            }
            let distance;
            if !sublen.is_null() {
                unsafe {
                    (*s.lmc).fetch_sublen(lmcpos, length as usize, sublen);
                    distance = *sublen.offset(length as isize);
                }

                if limit == ZOPFLI_MAX_MATCH && length >= ZOPFLI_MIN_MATCH {
                    unsafe {
                        assert!(*sublen.offset(length as isize) == dist_lmcpos);
                    }
                }
            } else {
                distance = dist_lmcpos;
            }
            longest_match.distance = distance;
            longest_match.length = length;
            longest_match.from_cache = 1;
            return longest_match;
        }
        /* Can't use much of the cache, since the "sublens" need to be calculated,
        but at least we already know when to stop. */
        limit = length_lmcpos as size_t;
        longest_match.limit = limit;
    }

    longest_match
}

/// Stores the found sublen, distance and length in the longest match cache, if
/// possible.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn StoreInLongestMatchCache(s_ptr: *mut ZopfliBlockState, pos: size_t, limit: size_t, sublen: *mut c_ushort, distance: c_ushort, length: c_ushort) {
    let s = unsafe {
        assert!(!s_ptr.is_null());
        &mut *s_ptr
    };

    /* The LMC cache starts at the beginning of the block rather than the
    beginning of the whole array. */
    let lmcpos = pos - s.blockstart;

    if s.lmc.is_null() {
        return;
    }

    /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
    that this cache value is not filled in yet. */
    let mut length_lmcpos = unsafe { (*s.lmc).length_at(lmcpos) };
    let mut dist_lmcpos = unsafe { (*s.lmc).dist_at(lmcpos) };

    let cache_available = length_lmcpos == 0 || dist_lmcpos != 0;

    if limit == ZOPFLI_MAX_MATCH && !sublen.is_null() && !cache_available {
        assert!(length_lmcpos == 1 && dist_lmcpos == 0);
        if length < ZOPFLI_MIN_MATCH {
            dist_lmcpos = 0;
            length_lmcpos = 0;
        } else {
            dist_lmcpos = distance;
            length_lmcpos = length;
        }
        assert!(!(length_lmcpos == 1 && dist_lmcpos == 0));
        unsafe {
            (*s.lmc).store_sublen(sublen, lmcpos, length as size_t);
        }
    }
}

/// Finds how long the match of scan and match is. Can be used to find how many
/// bytes starting from scan, and from match, are equal. Returns the last byte
/// after scan, which is still equal to the corresponding byte after match.
/// scan is the position to compare match is the earlier position to compare.
/// end is the last possible byte, beyond which to stop looking.
/// safe_end is a few (8) bytes before end, for comparing multiple bytes at once.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetMatch(mut scan: *mut c_uchar, mut match_char: *mut c_uchar, end: *mut c_uchar, safe_end: *mut c_uchar) -> *mut c_uchar {
    unsafe {
        /* 8 checks at once per array bounds check (size_t is 64-bit). */
        // C code has other options if size_t is not 64-bit, but this is all I'm supporting
        while scan < safe_end && *(scan as *const u64) == *(match_char as *const u64) {
            scan = scan.offset(8);
            match_char = match_char.offset(8);
        }

        /* The remaining few bytes. */
        while scan != end && *scan == *match_char {
            scan = scan.offset(1);
            match_char = match_char.offset(1);
        }

        scan
    }
}
