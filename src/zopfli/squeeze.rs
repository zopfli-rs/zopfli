use std::{mem, slice, ptr, cmp};

use libc::{c_void, c_uint, c_double, c_int, size_t, c_uchar, c_ushort, malloc, c_float};

use hash::ZopfliHash;
use lz77::{ZopfliLZ77Store, Lz77Store, ZopfliBlockState, find_longest_match, lz77_store_from_c,     lz77_store_result, verify_len_dist};
use symbols::{ZopfliGetDistExtraBits, ZopfliGetLengthExtraBits, ZopfliGetLengthSymbol, ZopfliGetDistSymbol, ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_LARGE_FLOAT, ZOPFLI_WINDOW_SIZE, ZOPFLI_WINDOW_MASK, ZOPFLI_MAX_MATCH, ZOPFLI_MIN_MATCH};

const K_INV_LOG2: c_double = 1.4426950408889;  // 1.0 / log(2.0)

/// Cost model which should exactly match fixed tree.
/// type: CostModelFun
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetCostFixed(litlen: c_uint, dist: c_uint, _unused: c_void) -> c_double {
    let result = if dist == 0 {
        if litlen <= 143 {
            8
        } else {
            9
        }
    } else {
        let dbits = ZopfliGetDistExtraBits(dist as c_int);
        let lbits = ZopfliGetLengthExtraBits(litlen as c_int);
        let lsym = ZopfliGetLengthSymbol(litlen as c_int);
        let mut cost = 0;
        if lsym <= 279 {
            cost += 7;
        } else {
            cost += 8;
        }
        cost += 5;  // Every dist symbol has length 5.
        cost + dbits + lbits
    };
    result as c_double
}

/// Cost model based on symbol statistics.
/// type: CostModelFun
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetCostStat(litlen: c_uint, dist: c_uint, context: *const c_void) -> c_double {
    let stats = unsafe {
        assert!(!context.is_null());
        &*(context as *const SymbolStats)
    };
    if dist == 0 {
        stats.ll_symbols[litlen as usize]
    } else {
        let lsym = ZopfliGetLengthSymbol(litlen as c_int) as usize;
        let lbits = ZopfliGetLengthExtraBits(litlen as c_int) as c_double;
        let dsym = ZopfliGetDistSymbol(dist as c_int) as usize;
        let dbits = ZopfliGetDistExtraBits(dist as c_int) as c_double;
        stats.ll_symbols[lsym] + lbits + stats.d_symbols[dsym] + dbits
    }
}

pub struct RanState {
    m_w: u32,
    m_z: u32,
}

impl RanState {
    pub fn new() -> RanState {
        RanState {
            m_w: 1,
            m_z: 2,
        }
    }

    /// Get random number: "Multiply-With-Carry" generator of G. Marsaglia
    pub fn random_marsaglia(&mut self) -> u32 {
        self.m_z = 36969 * (self.m_z & 65535) + (self.m_z >> 16);
        self.m_w = 18000 * (self.m_w & 65535) + (self.m_w >> 16);
        (self.m_z << 16).wrapping_add(self.m_w) // 32-bit result.
    }
}

#[no_mangle]
pub extern fn ran_state_new() -> *mut RanState {
    Box::into_raw(Box::new(RanState::new()))
}

#[derive(Copy)]
pub struct SymbolStats {
  /* The literal and length symbols. */
  litlens: [size_t; ZOPFLI_NUM_LL],
  /* The 32 unique dist symbols, not the 32768 possible dists. */
  dists: [size_t; ZOPFLI_NUM_D],

  /* Length of each lit/len symbol in bits. */
  ll_symbols: [c_double; ZOPFLI_NUM_LL],
  /* Length of each dist symbol in bits. */
  d_symbols: [c_double; ZOPFLI_NUM_D],
}

impl Clone for SymbolStats {
    fn clone(&self) -> Self {
        *self
    }
}

impl SymbolStats {
    pub fn new() -> SymbolStats {
        SymbolStats {
            litlens: [0; ZOPFLI_NUM_LL],
            dists: [0; ZOPFLI_NUM_D],
            ll_symbols: [0.0; ZOPFLI_NUM_LL],
            d_symbols: [0.0; ZOPFLI_NUM_D],
        }
    }

    pub fn randomize_stat_freqs(&mut self, state: &mut RanState) {
        fn randomize_freqs(freqs: &mut [size_t], state: &mut RanState) {
            let n = freqs.len();
            let mut i: usize = 0;
            let end = n as usize;

            while i < end {
                if (state.random_marsaglia() >> 4) % 3 == 0 {
                    let index = state.random_marsaglia() % n as c_uint;
                    freqs[i] = freqs[index as usize];
                }
                i += 1;
            }
        }
        randomize_freqs(&mut self.litlens, state);
        randomize_freqs(&mut self.dists, state);
        self.litlens[256] = 1; // End symbol.
    }

    /// Calculates the entropy of each symbol, based on the counts of each symbol. The
    /// result is similar to the result of ZopfliCalculateBitLengths, but with the
    /// actual theoritical bit lengths according to the entropy. Since the resulting
    /// values are fractional, they cannot be used to encode the tree specified by
    /// DEFLATE.
    pub fn calculate_entropy(&mut self) {
        fn calculate_and_store_entropy(count: &[size_t], bitlengths: &mut [c_double]) {
            let n = count.len();

            let mut sum = 0;
            for i in 0..n {
                sum += count[i];
            }

            let log2sum = (if sum == 0 { n } else { sum } as c_double).ln() * K_INV_LOG2;

            for i in 0..n {
                // When the count of the symbol is 0, but its cost is requested anyway, it
                // means the symbol will appear at least once anyway, so give it the cost as if
                // its count is 1.
                if count[i] == 0 {
                    bitlengths[i] = log2sum;
                } else {
                    bitlengths[i] = log2sum - (count[i] as c_double).ln() * K_INV_LOG2;
                }

                // Depending on compiler and architecture, the above subtraction of two
                // floating point numbers may give a negative result very close to zero
                // instead of zero (e.g. -5.973954e-17 with gcc 4.1.2 on Ubuntu 11.4). Clamp
                // it to zero. These floating point imprecisions do not affect the cost model
                // significantly so this is ok.
                if bitlengths[i] < 0.0 && bitlengths[i] > -1E-5 {
                    bitlengths[i] = 0.0;
                }
                assert!(bitlengths[i] >= 0.0);
            }
        }

        calculate_and_store_entropy(&self.litlens, &mut self.ll_symbols);
        calculate_and_store_entropy(&self.dists, &mut self.d_symbols);
    }

    /// Appends the symbol statistics from the store.
    pub fn get_statistics(&mut self, store: Lz77Store) {
        for i in 0..store.dists.len() {
            if store.dists[i] == 0 {
                self.litlens[store.litlens[i] as usize] += 1;
            } else {
                self.litlens[ZopfliGetLengthSymbol(store.litlens[i] as c_int) as usize] +=1 ;
                self.dists[ZopfliGetDistSymbol(store.dists[i] as c_int) as usize] += 1;
            }
        }
        self.litlens[256] = 1;  /* End symbol. */

        self.calculate_entropy();
    }
}

#[no_mangle]
pub extern fn symbol_stats_new() -> *mut SymbolStats {
    Box::into_raw(Box::new(SymbolStats::new()))
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ClearStatFreqs(stats_ptr: *mut SymbolStats) {
    let stats = unsafe {
        assert!(!stats_ptr.is_null());
        &mut *stats_ptr
    };
    stats.litlens = [0; ZOPFLI_NUM_LL];
    stats.dists = [0; ZOPFLI_NUM_D];
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CopyStats(source_ptr: *mut SymbolStats, dest_ptr: *mut SymbolStats) {
    let source = unsafe {
        assert!(!source_ptr.is_null());
        &mut *source_ptr
    };
    let dest = unsafe {
        assert!(!dest_ptr.is_null());
        &mut *dest_ptr
    };
    *dest = *source;
}

/// Adds the bit lengths.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn AddWeighedStatFreqs(stats1_ptr: *mut SymbolStats, w1: c_double, stats2_ptr: *mut SymbolStats, w2: c_double, result_ptr: *mut SymbolStats) {
    let stats1 = unsafe {
        assert!(!stats1_ptr.is_null());
        &mut *stats1_ptr
    };
    let stats2 = unsafe {
        assert!(!stats2_ptr.is_null());
        &mut *stats2_ptr
    };
    let result = unsafe {
        assert!(!result_ptr.is_null());
        &mut *result_ptr
    };

    for i in 0..ZOPFLI_NUM_LL {
        result.litlens[i] = (stats1.litlens[i] as c_double * w1 + stats2.litlens[i] as c_double * w2) as size_t;
    }
    for i in 0..ZOPFLI_NUM_D {
        result.dists[i] = (stats1.dists[i] as c_double * w1 + stats2.dists[i] as c_double * w2) as size_t;
    }
    result.litlens[256] = 1; // End symbol.
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn RandomizeStatFreqs(state_ptr: *mut RanState, stats_ptr: *mut SymbolStats) {
    let stats = unsafe {
        assert!(!stats_ptr.is_null());
        &mut *stats_ptr
    };
    let state = unsafe {
        assert!(!state_ptr.is_null());
        &mut *state_ptr
    };

    stats.randomize_stat_freqs(state);
}

/// Calculates the entropy of the statistics
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CalculateStatistics(stats_ptr: *mut SymbolStats) {
    let stats = unsafe {
        assert!(!stats_ptr.is_null());
        &mut *stats_ptr
    };

    stats.calculate_entropy();
}

/// Appends the symbol statistics from the store.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetStatistics(store_ptr: *const ZopfliLZ77Store, stats_ptr: *mut SymbolStats) {
    let store: Lz77Store = store_ptr.into();
    let stats = unsafe {
        assert!(!stats_ptr.is_null());
        &mut *stats_ptr
    };
    stats.get_statistics(store);
}

/// Finds the minimum possible cost this cost model can return for valid length and
/// distance symbols.
pub fn get_cost_model_min_cost(costmodel: fn(c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void) -> c_double {
    let mut bestlength: c_int = 0; // length that has lowest cost in the cost model
    let mut bestdist: c_int = 0; // distance that has lowest cost in the cost model

    // Table of distances that have a different distance symbol in the deflate
    // specification. Each value is the first distance that has a new symbol. Only
    // different symbols affect the cost model so only these need to be checked.
    // See RFC 1951 section 3.2.5. Compressed blocks (length and distance codes).

    let dsymbols: [c_int; 30] = [
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513,
        769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
    ];

    let mut mincost: c_double = ZOPFLI_LARGE_FLOAT;
    for i in 3..259 {
        let c = costmodel(i, 1, costcontext);
        if c < mincost {
            bestlength = i as c_int;
            mincost = c;
        }
    }

    mincost = ZOPFLI_LARGE_FLOAT;
    for i in 0..30 {
        let c = costmodel(3, dsymbols[i] as c_uint, costcontext);
        if c < mincost {
            bestdist = dsymbols[i];
            mincost = c;
        }
    }
    costmodel(bestlength as c_uint, bestdist as c_uint, costcontext)
}

/// Performs the forward pass for "squeeze". Gets the most optimal length to reach
/// every byte from a previous byte, using cost calculations.
/// s: the ZopfliBlockState
/// in: the input data array
/// instart: where to start
/// inend: where to stop (not inclusive)
/// costmodel: function to calculate the cost of some lit/len/dist pair.
/// costcontext: abstract context for the costmodel function
/// length_array: output array of size (inend - instart) which will receive the best
///     length to reach this byte from a previous byte.
/// returns the cost that was, according to the costmodel, needed to get to the end.
// TODO: upstream is now reusing an already allocated hash; we're ignoring it
pub fn get_best_lengths(s_ptr: *mut ZopfliBlockState, in_data: *mut c_uchar, instart: size_t, inend: size_t, costmodel: fn (c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void, length_array: *mut c_ushort, h_ptr: *mut ZopfliHash, costs: *mut c_float) -> c_double {
    let s = unsafe {
        assert!(!s_ptr.is_null());
        &mut *s_ptr
    };

    // Best cost to get here so far.
    let blocksize = inend - instart;
    let mut leng;
    let mut longest_match;
    let sublen = unsafe { malloc(mem::size_of::<c_ushort>() as size_t * 259) as *mut c_ushort };
    let windowstart = if instart > ZOPFLI_WINDOW_SIZE {
        instart - ZOPFLI_WINDOW_SIZE
    } else {
        0
    };

    let mincost = get_cost_model_min_cost(costmodel, costcontext);

    if instart == inend {
        return 0.0;
    }

    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.reset(ZOPFLI_WINDOW_SIZE);
    let arr = unsafe { slice::from_raw_parts(in_data, inend) };
    h.warmup(arr, windowstart, inend);
    for i in windowstart..instart {
        h.update(arr, i);
    }

    unsafe {
        for i in 1..(blocksize + 1) {
            *costs.offset(i as isize) = ZOPFLI_LARGE_FLOAT as c_float;
        }
        *costs.offset(0) = 0.0; /* Because it's the start. */
        *length_array.offset(0) = 0;
    }

    let mut i = instart;
    while i < inend {
        let mut j = i - instart;  // Index in the costs array and length_array.
        h.update(arr, i);

        // If we're in a long repetition of the same character and have more than
        // ZOPFLI_MAX_MATCH characters before and after our position.
        if h.same[i & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH as c_ushort * 2
            && i > instart + ZOPFLI_MAX_MATCH + 1
            && i + ZOPFLI_MAX_MATCH * 2 + 1 < inend
            && h.same[(i - ZOPFLI_MAX_MATCH) & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH as c_ushort {

            let symbolcost = costmodel(ZOPFLI_MAX_MATCH as c_uint, 1, costcontext);
            // Set the length to reach each one to ZOPFLI_MAX_MATCH, and the cost to
            // the cost corresponding to that length. Doing this, we skip
            // ZOPFLI_MAX_MATCH values to avoid calling ZopfliFindLongestMatch.

            for _ in 0..ZOPFLI_MAX_MATCH {
                unsafe {
                    *costs.offset((j + ZOPFLI_MAX_MATCH) as isize) = *costs.offset(j as isize) + symbolcost as c_float;
                    *length_array.offset((j + ZOPFLI_MAX_MATCH) as isize) = ZOPFLI_MAX_MATCH as c_ushort;
                }
                i += 1;
                j += 1;
                h.update(arr, i);
            }
        }

        longest_match = find_longest_match(s, h, arr, i, inend, ZOPFLI_MAX_MATCH, sublen);
        leng = longest_match.length;

        // Literal.
        if i + 1 <= inend {
            let new_cost = costmodel(arr[i] as c_uint, 0, costcontext) + unsafe { *costs.offset(j as isize) } as c_double;
            assert!(new_cost >= 0.0);
            if new_cost < unsafe { *costs.offset(j as isize + 1) } as c_double {
                unsafe {
                    *costs.offset(j as isize + 1) = new_cost as c_float;
                    *length_array.offset((j + 1) as isize) = 1;
                }
            }
        }
        // Lengths.
        let kend = cmp::min(leng as size_t, inend - i) as usize;
        let mincostaddcostj = mincost + unsafe { *costs.offset(j as isize) } as c_double;

        for k in 3..(kend + 1) {
            // Calling the cost model is expensive, avoid this if we are already at
            // the minimum possible cost that it can return.
            if unsafe { *costs.offset((j + k) as isize) } as c_double <= mincostaddcostj {
                continue;
            }

            let new_cost = costmodel(k as c_uint, unsafe { *sublen.offset(k as isize) } as c_uint, costcontext) + unsafe { *costs.offset(j as isize) } as c_double;
            assert!(new_cost >= 0.0);
            if new_cost < unsafe { *costs.offset((j + k) as isize) } as c_double {
                assert!(k as usize <= ZOPFLI_MAX_MATCH);
                unsafe {
                    *costs.offset((j + k) as isize) = new_cost as c_float;
                    *length_array.offset((j + k) as isize) = k as c_ushort;
                }
            }
        }
        i += 1;
    }

    assert!(unsafe { *costs.offset(blocksize as isize) } >= 0.0);
    unsafe { *costs.offset(blocksize as isize) as c_double }
}

// TODO: upstream is now reusing an already allocated hash; we're ignoring it
pub fn follow_path(s_ptr: *mut ZopfliBlockState, in_data: *const c_uchar, instart: size_t, inend: size_t, path: Vec<c_ushort>, store_ptr: *mut ZopfliLZ77Store, _h_ptr: *mut ZopfliHash) {
    let s = unsafe {
        assert!(!s_ptr.is_null());
        &mut *s_ptr
    };
    let store = unsafe {
        assert!(!store_ptr.is_null());
        &mut *store_ptr
    };
    let rust_store = lz77_store_from_c(store_ptr);

    let windowstart = if instart > ZOPFLI_WINDOW_SIZE {
        instart - ZOPFLI_WINDOW_SIZE
    } else {
        0
    };

    if instart == inend {
        return;
    }

    let mut h = ZopfliHash::new(ZOPFLI_WINDOW_SIZE);

    let arr = unsafe { slice::from_raw_parts(in_data, inend) };
    h.warmup(arr, windowstart, inend);

    for i in windowstart..instart {
        h.update(arr, i);
    }

    let mut pos = instart;
    let pathsize = path.len();
    for i in 0..pathsize {
        let mut length = path[i];
        assert!(pos < inend);

        h.update(arr, pos);

        // Add to output.
        if length >= ZOPFLI_MIN_MATCH as c_ushort {
            // Get the distance by recalculating longest match. The found length
            // should match the length from the path.
            let longest_match = find_longest_match(s, &mut h, arr, pos, inend, length as size_t, ptr::null_mut());
            let dist = longest_match.distance;
            let dummy_length = longest_match.length;
            assert!(!(dummy_length != length && length > 2 && dummy_length > 2));
            verify_len_dist(arr, pos, dist, length);
            unsafe {
                (&mut *rust_store).lit_len_dist(length, dist, pos);
            }
        } else {
            length = 1;
            unsafe {
                (&mut *rust_store).lit_len_dist(arr[pos] as c_ushort, 0, pos);
            }
        }

        assert!(pos + (length as size_t) <= inend);
        for j in 1..(length as size_t) {
            h.update(arr, pos + j);
        }

        pos += length as size_t;
    }
    lz77_store_result(rust_store, store);
}

/// Calculates the optimal path of lz77 lengths to use, from the calculated
/// length_array. The length_array must contain the optimal length to reach that
/// byte. The path will be filled with the lengths to use, so its data size will be
/// the amount of lz77 symbols.
pub fn trace_backwards(size: size_t, length_array: *const c_ushort) -> Vec<c_ushort> {
    let mut index = size;
    if size == 0 {
        return vec![];
    }
    let mut path = vec![]; // TODO: with capacity
    loop {
        let lai = unsafe { *length_array.offset(index as isize) };
        let laiu = lai as usize;
        path.push(lai);
        assert!(laiu <= index);
        assert!(laiu <= ZOPFLI_MAX_MATCH);
        assert!(lai != 0);
        index -= laiu;
        if index == 0 {
            break;
        }
    }

    /* Mirror result. */
    let pathsize = path.len();
    for index in 0..(pathsize / 2) {
        let temp = path[index];
        path[index] = path[pathsize - index - 1];
        path[pathsize - index - 1] = temp;
    }
    path
}

/// Does a single run for ZopfliLZ77Optimal. For good compression, repeated runs
/// with updated statistics should be performed.
/// s: the block state
/// in: the input data array
/// instart: where to start
/// inend: where to stop (not inclusive)
/// length_array: array of size (inend - instart) used to store lengths
/// costmodel: function to use as the cost model for this squeeze run
/// costcontext: abstract context for the costmodel function
/// store: place to output the LZ77 data
/// returns the cost that was, according to the costmodel, needed to get to the end.
///     This is not the actual cost.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn LZ77OptimalRun(s_ptr: *mut ZopfliBlockState, in_data: *mut c_uchar, instart: size_t, inend: size_t, length_array: *mut c_ushort, costmodel: fn (c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void, store_ptr: *mut ZopfliLZ77Store, h_ptr: *mut ZopfliHash, costs: *mut c_float) {

    let cost = get_best_lengths(s_ptr, in_data, instart, inend, costmodel, costcontext, length_array, h_ptr, costs);
    let path = trace_backwards(inend - instart, length_array);
    follow_path(s_ptr, in_data, instart, inend, path, store_ptr, h_ptr);
    assert!(cost < ZOPFLI_LARGE_FLOAT);
}
