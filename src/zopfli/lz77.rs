use std::{ptr, cmp};

use libc::{size_t, c_ushort, c_uchar, c_int, c_uint};

use cache::{ZopfliLongestMatchCache};
use hash::ZopfliHash;
use symbols::{get_dist_symbol, get_length_symbol};
use util::{ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_MAX_MATCH, ZOPFLI_MIN_MATCH, ZOPFLI_WINDOW_MASK, ZOPFLI_MAX_CHAIN_HITS, ZOPFLI_WINDOW_SIZE};
use zopfli::ZopfliOptions;

/// Stores lit/length and dist pairs for LZ77.
/// Parameter litlens: Contains the literal symbols or length values.
/// Parameter dists: Contains the distances. A value is 0 to indicate that there is
/// no dist and the corresponding litlens value is a literal instead of a length.
/// Parameter size: The size of both the litlens and dists arrays.
#[derive(Debug, Clone)]
pub struct Lz77Store {
   pub litlens: Vec<c_ushort>,
   pub dists: Vec<c_ushort>,

   pub pos: Vec<size_t>,

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

    pub fn reset(&mut self) {
        self.litlens.clear();
        self.dists.clear();
        self.pos.clear();
        self.ll_symbol.clear();
        self.d_symbol.clear();
        self.ll_counts.clear();
        self.d_counts.clear();
    }

    pub fn size(&self) -> usize {
        self.litlens.len()
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
            self.ll_symbol.push(get_length_symbol(length as c_int) as c_ushort);
            self.d_symbol.push(get_dist_symbol(dist as c_int) as c_ushort);
            self.ll_counts[llstart + get_length_symbol(length as c_int) as usize] += 1;
            self.d_counts[dstart + get_dist_symbol(dist as c_int) as usize] += 1;
        }
    }

    /// Does LZ77 using an algorithm similar to gzip, with lazy matching, rather than
    /// with the slow but better "squeeze" implementation.
    /// The result is placed in the Lz77Store.
    /// If instart is larger than 0, it uses values before instart as starting
    /// dictionary.
    pub fn greedy(&mut self, s: &mut ZopfliBlockState, in_data: &[u8], instart: size_t, inend: size_t) {
        let mut leng: c_ushort;
        let mut dist: c_ushort;
        let mut lengthscore: c_int;
        let windowstart = if instart > ZOPFLI_WINDOW_SIZE {
            instart - ZOPFLI_WINDOW_SIZE
        } else {
            0
        };

        let mut longest_match;

        /* Lazy matching. */
        let mut prev_length: c_uint = 0;
        let mut prev_match: c_uint = 0;
        let mut prevlengthscore: c_int;
        let mut match_available = false;

        if instart == inend {
            return;
        }

        let mut h = ZopfliHash::new(ZOPFLI_WINDOW_SIZE);

        let arr = &in_data[..inend];
        h.warmup(arr, windowstart, inend);

        for i in windowstart..instart {
            h.update(arr, i);
        }

        let mut i = instart;
        while i < inend {
            h.update(arr, i);

            longest_match = find_longest_match(s, &mut h, arr, i, inend, ZOPFLI_MAX_MATCH, ptr::null_mut());
            dist = longest_match.distance;
            leng = longest_match.length;
            lengthscore = get_length_score(leng as c_int, dist as c_int);

            /* Lazy matching. */
            prevlengthscore = get_length_score(prev_length as c_int, prev_match as c_int);
            if match_available {
                match_available = false;
                if lengthscore > prevlengthscore + 1 {
                    self.lit_len_dist(arr[i - 1] as c_ushort, 0, i - 1);
                    if (lengthscore as size_t) >= ZOPFLI_MIN_MATCH && (leng as size_t) < ZOPFLI_MAX_MATCH {
                        match_available = true;
                        prev_length = leng as c_uint;
                        prev_match = dist as c_uint;
                        i += 1;
                        continue;
                    }
                } else {
                    /* Add previous to output. */
                    leng = prev_length as c_ushort;
                    dist = prev_match as c_ushort;
                    /* Add to output. */
                    verify_len_dist(arr, i - 1, dist, leng);
                    self.lit_len_dist(leng, dist, i - 1);
                    for _ in 2..leng {
                        assert!(i < inend);
                        i += 1;
                        h.update(arr, i);
                     }
                     i += 1;
                     continue;
                }
            } else if (lengthscore as size_t) >= ZOPFLI_MIN_MATCH && (leng as size_t) < ZOPFLI_MAX_MATCH {
                match_available = true;
                prev_length = leng as c_uint;
                prev_match = dist as c_uint;
                i += 1;
                continue;
            }
            /* End of lazy matching. */

            /* Add to output. */
            if (lengthscore as size_t) >= ZOPFLI_MIN_MATCH {
                verify_len_dist(arr, i, dist, leng);
                self.lit_len_dist(leng, dist, i);
            } else {
                leng = 1;
                self.lit_len_dist(arr[i] as c_ushort, 0, i);
            }
            for _ in 1..leng {
                assert!(i < inend);
                i += 1;
                h.update(arr, i);
            }
            i += 1;
        }
    }

    pub fn follow_path(&mut self, in_data: &[u8], instart: size_t, inend: size_t, path: Vec<c_ushort>, s: &mut ZopfliBlockState) {
        let windowstart = if instart > ZOPFLI_WINDOW_SIZE {
            instart - ZOPFLI_WINDOW_SIZE
        } else {
            0
        };

        if instart == inend {
            return;
        }

        let mut h = ZopfliHash::new(ZOPFLI_WINDOW_SIZE);

        let arr = &in_data[..inend];
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
                self.lit_len_dist(length, dist, pos);
            } else {
                length = 1;
                self.lit_len_dist(arr[pos] as c_ushort, 0, pos);
            }

            assert!(pos + (length as size_t) <= inend);
            for j in 1..(length as size_t) {
                h.update(arr, pos + j);
            }

            pos += length as size_t;
        }
    }
}

/// Some state information for compressing a block.
/// This is currently a bit under-used (with mainly only the longest match cache),
/// but is kept for easy future expansion.
pub struct ZopfliBlockState<'a> {
    pub options: &'a ZopfliOptions,
    /* Cache for length/distance pairs found so far. */
    lmc: Option<ZopfliLongestMatchCache>,
    /* The start (inclusive) and end (not inclusive) of the current block. */
    pub blockstart: size_t,
    pub blockend: size_t,
}

impl<'a> ZopfliBlockState<'a> {
    pub fn new(options: &ZopfliOptions, blockstart: size_t, blockend: size_t, add_lmc: c_int) -> ZopfliBlockState {
        ZopfliBlockState {
            options: options,
            blockstart: blockstart,
            blockend: blockend,
            lmc: if add_lmc > 0 {
                Some(ZopfliLongestMatchCache::new(blockend - blockstart))
            } else {
                None
            },
        }
    }

    /// Gets distance, length and sublen values from the cache if possible.
    /// Returns 1 if it got the values from the cache, 0 if not.
    /// Updates the limit value to a smaller one if possible with more limited
    /// information from the cache.
    pub fn try_get_from_longest_match_cache(&self, pos: size_t, mut limit: size_t, sublen: *mut c_ushort) -> LongestMatch {
        let mut longest_match = LongestMatch {
            distance: 0,
            length: 0,
            from_cache: 0,
            limit: limit,
        };

        match self.lmc {
            None => return longest_match,
            Some(ref lmc) => {
                /* The LMC cache starts at the beginning of the block rather than the
                 beginning of the whole array. */
                let lmcpos = pos - self.blockstart;

                /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
                  that this cache value is not filled in yet. */
                let length_lmcpos = lmc.length_at(lmcpos);
                let dist_lmcpos = lmc.dist_at(lmcpos);
                let cache_available = length_lmcpos == 0 || dist_lmcpos != 0;
                let max_sublen = lmc.max_sublen(lmcpos);
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
                            lmc.fetch_sublen(lmcpos, length as usize, sublen);
                            unsafe {
                                distance = *sublen.offset(length as isize);
                            }

                            if limit == ZOPFLI_MAX_MATCH && length >= ZOPFLI_MIN_MATCH as c_ushort {
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
        }
    }

    /// Stores the found sublen, distance and length in the longest match cache, if
    /// possible.
    pub fn store_in_longest_match_cache(&mut self, pos: size_t, limit: size_t, sublen: *mut c_ushort, distance: c_ushort, length: c_ushort) {
        /* The LMC cache starts at the beginning of the block rather than the
        beginning of the whole array. */
        let lmcpos = pos - self.blockstart;

        match self.lmc {
            None => return,
            Some(ref mut lmc) => {
                /* Length > 0 and dist 0 is invalid combination, which indicates on purpose
                that this cache value is not filled in yet. */
                let mut length_lmcpos = lmc.length_at(lmcpos);
                let mut dist_lmcpos = lmc.dist_at(lmcpos);

                let cache_available = length_lmcpos == 0 || dist_lmcpos != 0;

                if limit == ZOPFLI_MAX_MATCH && !sublen.is_null() && !cache_available {
                    assert!(length_lmcpos == 1 && dist_lmcpos == 0);
                    if length < ZOPFLI_MIN_MATCH as c_ushort {
                        dist_lmcpos = 0;
                        length_lmcpos = 0;
                    } else {
                        dist_lmcpos = distance;
                        length_lmcpos = length;
                    }
                    assert!(!(length_lmcpos == 1 && dist_lmcpos == 0));
                    lmc.store_sublen(sublen, lmcpos, length as size_t);
                }
            }
        }
    }
}

pub struct LongestMatch {
    pub distance: c_ushort,
    pub length: c_ushort,
    from_cache: c_int,
    limit: size_t,
}

/// Finds how long the match of scan and match is. Can be used to find how many
/// bytes starting from scan, and from match, are equal. Returns the last byte
/// after scan, which is still equal to the corresponding byte after match.
/// scan is the position to compare match is the earlier position to compare.
/// end is the last possible byte, beyond which to stop looking.
/// safe_end is a few (8) bytes before end, for comparing multiple bytes at once.
pub fn get_match(array: &[c_uchar], scan_offset: usize, match_offset: usize, end: usize, _safe_end: usize) -> usize {
    let mut scan_offset = scan_offset;
    let mut match_offset = match_offset;
    // /* 8 checks at once per array bounds check (size_t is 64-bit). */
    // // C code has other options if size_t is not 64-bit, but this is all I'm supporting
    // while scan_offset < safe_end && array[scan_offset] as *const u64 == array[match_offset] as *const u64 {
    //     scan_offset += 8;
    //     match_offset += 8;
    // }

    /* The remaining few bytes. */
    while scan_offset != end && array[scan_offset] == array[match_offset] {
        scan_offset += 1;
        match_offset += 1;
    }

    scan_offset
}

pub fn find_longest_match(s: &mut ZopfliBlockState, h: &mut ZopfliHash, array: &[c_uchar], pos: size_t, size: size_t, limit: size_t, sublen: *mut c_ushort) -> LongestMatch {
    let mut limit = limit;

    let hpos = pos & ZOPFLI_WINDOW_MASK;
    let mut bestdist = 0;
    let mut bestlength = 1;
    let mut which_hash = 1;
    let mut chain_counter = ZOPFLI_MAX_CHAIN_HITS;  /* For quitting early. */
    let mut dist;  /* Not unsigned short on purpose. */

    let mut longest_match = s.try_get_from_longest_match_cache(pos, limit, sublen);

    if longest_match.from_cache == 1 {
        assert!(pos + (longest_match.length as size_t) <= size);
        return longest_match;
    }

    limit = longest_match.limit;

    assert!(limit <= ZOPFLI_MAX_MATCH);
    assert!(limit >= ZOPFLI_MIN_MATCH);
    assert!(pos < size);

    if size - pos < ZOPFLI_MIN_MATCH {
        /* The rest of the code assumes there are at least ZOPFLI_MIN_MATCH bytes to
        try. */
        longest_match.distance = 0;
        longest_match.length = 0;
        longest_match.from_cache = 0;
        longest_match.limit = 0;
        return longest_match;
    }

    if pos + limit > size {
        limit = size - pos;
    }
    let arrayend = pos + limit;
    let arrayend_safe = arrayend - 8;

    assert!(h.val(which_hash) < 65536);

    let mut pp = h.head_at(h.val(which_hash) as usize, which_hash);  /* During the whole loop, p == hprev[pp]. */
    let mut p = h.prev_at(pp as usize, which_hash);

    assert!(pp as size_t == hpos);

    dist = if (p as c_int) < pp {
        pp - (p as c_int)
    } else {
        (ZOPFLI_WINDOW_SIZE - (p as size_t)) as c_int + pp
    };

    let mut scan_offset;
    let mut match_offset;

    /* Go through all distances. */
    while (dist as size_t) < ZOPFLI_WINDOW_SIZE {
        let mut currentlength = 0;

        assert!((p as size_t) < ZOPFLI_WINDOW_SIZE);
        assert!(p == h.prev_at(pp as usize, which_hash));
        assert!(h.hash_val_at(p as usize, which_hash) == h.val(which_hash));

        if dist > 0 {
            assert!(pos < size);
            assert!((dist as size_t) <= pos);
            scan_offset = pos;
            match_offset = pos - (dist as size_t);

            /* Testing the byte at position bestlength first, goes slightly faster. */
            if pos + bestlength >= size || array[scan_offset + bestlength] == array[match_offset + bestlength] {

                let same0 = h.same[pos & ZOPFLI_WINDOW_MASK];
                if same0 > 2 && array[scan_offset] == array[match_offset] {
                    let same1 = h.same[(pos - (dist as size_t)) & ZOPFLI_WINDOW_MASK];
                    let mut same = if same0 < same1 {
                        same0
                    } else {
                        same1
                    };
                    if (same as size_t) > limit {
                        same = limit as c_ushort;
                    }
                    scan_offset += same as size_t;
                    match_offset += same as size_t;
                }
                scan_offset = get_match(array, scan_offset, match_offset, arrayend, arrayend_safe);
                currentlength = scan_offset - pos;  /* The found length. */
            }

            if currentlength > bestlength {
                if !sublen.is_null() {
                    for j in (bestlength + 1)..(currentlength + 1) {
                        unsafe {
                            *sublen.offset(j as isize) = dist as c_ushort;
                        }
                    }
                }
                bestdist = dist;
                bestlength = currentlength;
                if currentlength >= limit {
                    break;
                }
            }
        }

        /* Switch to the other hash once this will be more efficient. */
        if which_hash == 1 && bestlength >= h.same[hpos] as size_t && h.val(2) == h.hash_val_at(p as usize, 2) {
            /* Now use the hash that encodes the length and first byte. */
            which_hash = 2;
        }

        pp = p as c_int;
        p = h.prev_at(p as usize, which_hash);
        if (p as c_int) == pp {
            break;  /* Uninited prev value. */
        }

        dist += if (p as c_int) < pp {
            pp - (p as c_int)
        } else {
            (ZOPFLI_WINDOW_SIZE - (p as usize)) as c_int + pp
        };

        chain_counter -= 1;
        if chain_counter <= 0 {
            break;
        }
    }

    s.store_in_longest_match_cache(pos, limit, sublen, bestdist as c_ushort, bestlength as c_ushort);

    assert!(bestlength <= limit);

    assert!(pos + bestlength <= size);
    longest_match.distance = bestdist as c_ushort;
    longest_match.length = bestlength as c_ushort;
    longest_match.from_cache = 0;
    longest_match.limit = limit;
    longest_match
}

/// Gets a score of the length given the distance. Typically, the score of the
/// length is the length itself, but if the distance is very long, decrease the
/// score of the length a bit to make up for the fact that long distances use large
/// amounts of extra bits.
///
/// This is not an accurate score, it is a heuristic only for the greedy LZ77
/// implementation. More accurate cost models are employed later. Making this
/// heuristic more accurate may hurt rather than improve compression.
///
/// The two direct uses of this heuristic are:
/// -avoid using a length of 3 in combination with a long distance. This only has
///  an effect if length == 3.
/// -make a slightly better choice between the two options of the lazy matching.
///
/// Indirectly, this affects:
/// -the block split points if the default of block splitting first is used, in a
///  rather unpredictable way
/// -the first zopfli run, so it affects the chance of the first run being closer
///  to the optimal output
pub fn get_length_score(length: c_int, distance: c_int) -> c_int {
    // At 1024, the distance uses 9+ extra bits and this seems to be the sweet spot
    // on tested files.
    if distance > 1024 {
        length - 1
    } else {
        length
    }
}

pub fn verify_len_dist(data: &[c_uchar], pos: size_t, dist: c_ushort, length: c_ushort) {
    for i in 0..length {
        let d1 = data[pos - (dist as usize) + (i as usize)];
        let d2 = data[pos + (i as usize)];
        if d1 != d2 {
            assert!(d1 == d2);
            break;
        }
    }
}

pub fn get_byte_range(lz77: &Lz77Store, lstart: size_t, lend: size_t) -> size_t {
    let l = lend - 1;
    if lstart == lend {
        return 0;
    }

    lz77.pos[l] + (if lz77.dists[l] == 0 { 1 } else { lz77.litlens[l] }) as size_t - lz77.pos[lstart]
}

pub fn get_histogram_at(lz77: &Lz77Store, lpos: size_t) -> (Vec<size_t>, Vec<size_t>) {
    let mut ll = vec![0; ZOPFLI_NUM_LL];
    let mut d = vec![0; ZOPFLI_NUM_D];

    /* The real histogram is created by using the histogram for this chunk, but
    all superfluous values of this chunk subtracted. */
    let llpos = ZOPFLI_NUM_LL * (lpos / ZOPFLI_NUM_LL);
    let dpos = ZOPFLI_NUM_D * (lpos / ZOPFLI_NUM_D);

    for i in 0..ZOPFLI_NUM_LL {
        ll[i] = lz77.ll_counts[llpos + i];
    }
    let end = cmp::min(llpos + ZOPFLI_NUM_LL, lz77.size());
    for i in (lpos + 1)..end {
        ll[lz77.ll_symbol[i] as usize] -= 1;
    }
    for i in 0..ZOPFLI_NUM_D {
        d[i] = lz77.d_counts[dpos + i];
    }
    let end = cmp::min(dpos + ZOPFLI_NUM_D, lz77.size());
    for i in (lpos + 1)..end {
        if lz77.dists[i] != 0 {
            d[lz77.d_symbol[i] as usize] -= 1;
        }
    }

    (ll, d)
}

/// Gets the histogram of lit/len and dist symbols in the given range, using the
/// cumulative histograms, so faster than adding one by one for large range. Does
/// not add the one end symbol of value 256.
pub fn get_histogram(lz77: &Lz77Store, lstart: size_t, lend: size_t) -> (Vec<size_t>, Vec<size_t>) {
    if lstart + ZOPFLI_NUM_LL * 3 > lend {
        let mut ll_counts = vec![0; ZOPFLI_NUM_LL];
        let mut d_counts = vec![0; ZOPFLI_NUM_D];
        for i in lstart..lend  {
            ll_counts[lz77.ll_symbol[i] as usize] += 1;
            if lz77.dists[i] != 0 {
                d_counts[lz77.d_symbol[i] as usize] += 1;
            }
        }
        (ll_counts, d_counts)
    } else {
        /* Subtract the cumulative histograms at the end and the start to get the
        histogram for this range. */
        let (ll, d) = get_histogram_at(lz77, lend - 1);

        if lstart > 0 {
            let (ll2, d2) = get_histogram_at(lz77, lstart - 1);

            (
                ll.iter().zip(ll2.iter()).map(|(&ll_item1, &ll_item2)|
                    ll_item1 - ll_item2
                ).collect(),
                d.iter().zip(d2.iter()).map(|(&d_item1, &d_item2)|
                    d_item1 - d_item2
                ).collect(),
            )
        } else {
            (ll, d)
        }
    }
}
