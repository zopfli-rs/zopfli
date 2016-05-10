use libc::{c_void, c_uint, c_double, c_int, size_t};

use lz77::{ZopfliLZ77Store, Lz77Store};
use symbols::{ZopfliGetDistExtraBits, ZopfliGetLengthExtraBits, ZopfliGetLengthSymbol, ZopfliGetDistSymbol, ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_LARGE_FLOAT};

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
    for i in 0..store.dists.len() {
        if store.dists[i] == 0 {
            stats.litlens[store.litlens[i] as usize] += 1;
        } else {
            stats.litlens[ZopfliGetLengthSymbol(store.litlens[i] as c_int) as usize] +=1 ;
            stats.dists[ZopfliGetDistSymbol(store.dists[i] as c_int) as usize] += 1;
        }
    }
    stats.litlens[256] = 1;  /* End symbol. */

    CalculateStatistics(stats);
}

/// Finds the minimum possible cost this cost model can return for valid length and
/// distance symbols.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetCostModelMinCost(costmodel: fn(c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void) -> c_double {
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
