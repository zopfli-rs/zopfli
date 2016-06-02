//! The squeeze functions do enhanced LZ77 compression by optimal parsing with a
//! cost model, rather than greedily choosing the longest length or using a single
//! step of lazy matching like regular implementations.
//!
//! Since the cost model is based on the Huffman tree that can only be calculated
//! after the LZ77 data is generated, there is a chicken and egg problem, and
//! multiple runs are done with updated cost models to converge to a better
//! solution.

use std::{mem, ptr, cmp};

use libc::{c_void, c_uint, c_double, c_int, size_t, c_ushort, malloc, c_float};

use deflate::calculate_block_size;
use hash::ZopfliHash;
use lz77::{Lz77Store, ZopfliBlockState, find_longest_match};
use symbols::{get_dist_extra_bits, get_dist_symbol, get_length_extra_bits, get_length_symbol};
use util::{ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_LARGE_FLOAT, ZOPFLI_WINDOW_SIZE, ZOPFLI_WINDOW_MASK, ZOPFLI_MAX_MATCH};

const K_INV_LOG2: c_double = 1.4426950408889;  // 1.0 / log(2.0)

/// Cost model which should exactly match fixed tree.
/// type: CostModelFun
#[allow(non_snake_case)]
pub fn GetCostFixed(litlen: c_uint, dist: c_uint, _unused: *const c_void) -> c_double {
    let result = if dist == 0 {
        if litlen <= 143 {
            8
        } else {
            9
        }
    } else {
        let dbits = get_dist_extra_bits(dist as c_int);
        let lbits = get_length_extra_bits(litlen as c_int);
        let lsym = get_length_symbol(litlen as c_int);
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
#[allow(non_snake_case)]
pub fn GetCostStat(litlen: c_uint, dist: c_uint, context: *const c_void) -> c_double {
    let stats = unsafe {
        assert!(!context.is_null());
        &*(context as *const SymbolStats)
    };
    if dist == 0 {
        stats.ll_symbols[litlen as usize]
    } else {
        let lsym = get_length_symbol(litlen as c_int) as usize;
        let lbits = get_length_extra_bits(litlen as c_int) as c_double;
        let dsym = get_dist_symbol(dist as c_int) as usize;
        let dbits = get_dist_extra_bits(dist as c_int) as c_double;
        lbits + dbits + stats.ll_symbols[lsym] + stats.d_symbols[dsym]
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
            let end = n;

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
    /// result is similar to the result of length_limited_code_lengths, but with the
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
    pub fn get_statistics(&mut self, store: &Lz77Store) {
        for i in 0..store.dists.len() {
            if store.dists[i] == 0 {
                self.litlens[store.litlens[i] as usize] += 1;
            } else {
                self.litlens[get_length_symbol(store.litlens[i] as c_int) as usize] +=1 ;
                self.dists[get_dist_symbol(store.dists[i] as c_int) as usize] += 1;
            }
        }
        self.litlens[256] = 1;  /* End symbol. */

        self.calculate_entropy();
    }

    pub fn clear_freqs(&mut self) {
        self.litlens = [0; ZOPFLI_NUM_LL];
        self.dists = [0; ZOPFLI_NUM_D];
    }
}

pub fn add_weighed_stat_freqs(stats1: &SymbolStats, w1: c_double, stats2: &SymbolStats, w2: c_double) -> SymbolStats {
    let mut result = SymbolStats::new();

    for i in 0..ZOPFLI_NUM_LL {
        result.litlens[i] = (stats1.litlens[i] as c_double * w1 + stats2.litlens[i] as c_double * w2) as size_t;
    }
    for i in 0..ZOPFLI_NUM_D {
        result.dists[i] = (stats1.dists[i] as c_double * w1 + stats2.dists[i] as c_double * w2) as size_t;
    }
    result.litlens[256] = 1; // End symbol.
    result
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
pub fn get_best_lengths(s: &mut ZopfliBlockState, in_data: &[u8], instart: size_t, inend: size_t, costmodel: fn (c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void, h: &mut ZopfliHash, costs: &mut Vec<c_float>) -> (c_double, Vec<c_ushort>) {
    // Best cost to get here so far.
    let blocksize = inend - instart;
    let mut length_array = vec![0; blocksize + 1];
    let mut leng;
    let mut longest_match;
    let sublen = unsafe { malloc(mem::size_of::<c_ushort>() * 259) as *mut c_ushort };
    let windowstart = if instart > ZOPFLI_WINDOW_SIZE {
        instart - ZOPFLI_WINDOW_SIZE
    } else {
        0
    };

    let mincost = get_cost_model_min_cost(costmodel, costcontext);

    if instart == inend {
        return (0.0, length_array);
    }

    h.reset(ZOPFLI_WINDOW_SIZE);
    let arr = &in_data[..inend];
    h.warmup(arr, windowstart, inend);
    for i in windowstart..instart {
        h.update(arr, i);
    }

    costs.resize(blocksize + 1, 0.0);

    for i in 1..(blocksize + 1) {
        costs[i] = ZOPFLI_LARGE_FLOAT as c_float;
    }
    costs[0] = 0.0; /* Because it's the start. */
    length_array[0] = 0;

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
                costs[j + ZOPFLI_MAX_MATCH] = costs[j] + symbolcost as c_float;
                length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH as c_ushort;
                i += 1;
                j += 1;
                h.update(arr, i);
            }
        }

        longest_match = find_longest_match(s, h, arr, i, inend, ZOPFLI_MAX_MATCH, sublen);
        leng = longest_match.length;

        // Literal.
        if i + 1 <= inend {
            let new_cost = costmodel(arr[i] as c_uint, 0, costcontext) + costs[j] as c_double;
            assert!(new_cost >= 0.0);
            if new_cost < costs[j + 1] as c_double {
                costs[j + 1] = new_cost as c_float;
                length_array[j + 1] = 1;
            }
        }
        // Lengths.
        let kend = cmp::min(leng as size_t, inend - i);
        let mincostaddcostj = mincost + costs[j] as c_double;

        for k in 3..(kend + 1) {
            // Calling the cost model is expensive, avoid this if we are already at
            // the minimum possible cost that it can return.
            if costs[j + k] as c_double <= mincostaddcostj {
                continue;
            }

            let new_cost = costmodel(k as c_uint, unsafe { *sublen.offset(k as isize) } as c_uint, costcontext) + costs[j] as c_double;
            assert!(new_cost >= 0.0);
            if new_cost < costs[j + k] as c_double {
                assert!(k <= ZOPFLI_MAX_MATCH);
                costs[j + k] = new_cost as c_float;
                length_array[j + k] = k as c_ushort;
            }
        }
        i += 1;
    }

    assert!(costs[blocksize] >= 0.0);
    (costs[blocksize] as c_double, length_array)
}


/// Calculates the optimal path of lz77 lengths to use, from the calculated
/// length_array. The length_array must contain the optimal length to reach that
/// byte. The path will be filled with the lengths to use, so its data size will be
/// the amount of lz77 symbols.
pub fn trace_backwards(size: size_t, length_array: Vec<c_ushort>) -> Vec<c_ushort> {
    let mut index = size;
    if size == 0 {
        return vec![];
    }
    let mut path = vec![]; // TODO: with capacity
    loop {
        let lai = length_array[index];
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

/// Does a single run for lz77_optimal. For good compression, repeated runs
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
pub fn lz77_optimal_run(s: &mut ZopfliBlockState, in_data: &[u8], instart: size_t, inend: size_t, costmodel: fn (c_uint, c_uint, *const c_void) -> c_double, costcontext: *const c_void, store: &mut Lz77Store, h: &mut ZopfliHash, costs: &mut Vec<c_float>) {
    let (cost, length_array) = get_best_lengths(s, in_data, instart, inend, costmodel, costcontext, h, costs);
    let path = trace_backwards(inend - instart, length_array);
    store.follow_path(in_data, instart, inend, path, s);
    assert!(cost < ZOPFLI_LARGE_FLOAT);
}


/// Does the same as lz77_optimal, but optimized for the fixed tree of the
/// deflate standard.
/// The fixed tree never gives the best compression. But this gives the best
/// possible LZ77 encoding possible with the fixed tree.
/// This does not create or output any fixed tree, only LZ77 data optimized for
/// using with a fixed tree.
/// If instart is larger than 0, it uses values before instart as starting
/// dictionary.
pub fn lz77_optimal_fixed(s: &mut ZopfliBlockState, in_data: &[u8], instart: size_t, inend: size_t, store: &mut Lz77Store) {
    s.blockstart = instart;
    s.blockend = inend;
    let mut h = ZopfliHash::new(ZOPFLI_WINDOW_SIZE);
    let mut costs = Vec::with_capacity(inend - instart - 1);
    lz77_optimal_run(s, in_data, instart, inend, GetCostFixed, ptr::null(), store, &mut h, &mut costs);
}

/// Calculates lit/len and dist pairs for given data.
/// If instart is larger than 0, it uses values before instart as starting
/// dictionary.
pub fn lz77_optimal(s: &mut ZopfliBlockState, in_data: &[u8], instart: size_t, inend: size_t, numiterations: c_int) -> Lz77Store {

    let mut h = ZopfliHash::new(ZOPFLI_WINDOW_SIZE);
    let mut costs = Vec::with_capacity(inend - instart + 1);

    /* Dist to get to here with smallest cost. */
    let mut currentstore = Lz77Store::new();
    let mut outputstore = currentstore.clone();

    let mut stats = SymbolStats::new();
    let mut beststats = SymbolStats::new();

    let mut bestcost = ZOPFLI_LARGE_FLOAT;
    let mut lastcost = 0.0;
    /* Try randomizing the costs a bit once the size stabilizes. */
    let mut ran_state = RanState::new();
    let mut lastrandomstep = -1;

    /* Do regular deflate, then loop multiple shortest path runs, each time using
    the statistics of the previous run. */

    /* Initial run. */
    currentstore.greedy(s, in_data, instart, inend);
    stats.get_statistics(&currentstore);

    /* Repeat statistics with each time the cost model from the previous stat
    run. */
    for i in 0..numiterations {
        currentstore.reset();
        let stats_ptr: *const SymbolStats = &stats;
        lz77_optimal_run(s, in_data, instart, inend, GetCostStat, stats_ptr as *const c_void, &mut currentstore, &mut h, &mut costs);
        let cost = calculate_block_size(&currentstore, 0, currentstore.size(), 2);

        if s.options.verbose_more != 0 || (s.options.verbose != 0 && cost < bestcost) {
              println!("Iteration {}: {} bit", i, cost);
        }
        if cost < bestcost {
            /* Copy to the output store. */
            outputstore = currentstore.clone();
            beststats = stats;
            bestcost = cost;
        }
        let laststats = stats;
        stats.clear_freqs();
        stats.get_statistics(&currentstore);
        if lastrandomstep != -1 {
            /* This makes it converge slower but better. Do it only once the
            randomness kicks in so that if the user does few iterations, it gives a
            better result sooner. */
            stats = add_weighed_stat_freqs(&stats, 1.0, &laststats, 0.5);
            stats.calculate_entropy();
        }
        if i > 5 && cost == lastcost {
            stats = beststats;
            stats.randomize_stat_freqs(&mut ran_state);
            stats.calculate_entropy();
            lastrandomstep = i;
        }
        lastcost = cost;
    }
    outputstore
}
