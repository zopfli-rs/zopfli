use std::cmp;

use cache::{ZopfliLongestMatchCache, Cache, NoCache};
use hash::{ZopfliHash, Which};
use symbols::{get_dist_symbol, get_length_symbol};
use util::{ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_MAX_MATCH, ZOPFLI_MIN_MATCH, ZOPFLI_WINDOW_MASK, ZOPFLI_MAX_CHAIN_HITS, ZOPFLI_WINDOW_SIZE};
use Options;

#[derive(Clone, Debug, Copy)]
pub enum LitLen {
    Literal(u16),
    LengthDist(u16, u16),
}

impl LitLen {
    pub fn size(&self) -> usize {
        match *self {
            LitLen::Literal(_) => 1,
            LitLen::LengthDist(len, _) => len as usize,
        }
    }
}

/// Stores lit/length and dist pairs for LZ77.
/// Parameter litlens: Contains the literal symbols or length values.
/// Parameter dists: Contains the distances. A value is 0 to indicate that there is
/// no dist and the corresponding litlens value is a literal instead of a length.
/// Parameter size: The size of both the litlens and dists arrays.
#[derive(Debug, Clone, Default)]
pub struct Lz77Store {
   pub litlens: Vec<LitLen>,

   pub pos: Vec<usize>,

   ll_symbol: Vec<u16>,
   d_symbol: Vec<u16>,

   ll_counts: Vec<usize>,
   d_counts: Vec<usize>,
}

impl Lz77Store {
    pub fn new() -> Lz77Store {
        Lz77Store {
          litlens: vec![],

          pos: vec![],

          ll_symbol: vec![],
          d_symbol: vec![],

          ll_counts: vec![],
          d_counts: vec![],
       }
    }

    pub fn reset(&mut self) {
        self.litlens.clear();
        self.pos.clear();
        self.ll_symbol.clear();
        self.d_symbol.clear();
        self.ll_counts.clear();
        self.d_counts.clear();
    }

    pub fn size(&self) -> usize {
        self.litlens.len()
    }

    pub fn append_store_item(&mut self, litlen: LitLen, pos: usize) {
        let origsize = self.litlens.len();
        let llstart = ZOPFLI_NUM_LL * (origsize / ZOPFLI_NUM_LL);
        let dstart = ZOPFLI_NUM_D * (origsize / ZOPFLI_NUM_D);

        if origsize % ZOPFLI_NUM_LL == 0 {
            if origsize == 0 {
                self.ll_counts.resize(origsize + ZOPFLI_NUM_LL, 0);
            } else {
                let mut last_histogram = (&self.ll_counts[(origsize - ZOPFLI_NUM_LL)..origsize]).to_vec();
                self.ll_counts.append(&mut last_histogram);
            }
        }

        if origsize % ZOPFLI_NUM_D == 0 {
            if origsize == 0 {
                self.d_counts.resize(ZOPFLI_NUM_D, 0);
            } else {
                let mut last_histogram = (&self.d_counts[(origsize - ZOPFLI_NUM_D)..origsize]).to_vec();
                self.d_counts.append(&mut last_histogram);
            }
        }

        self.pos.push(pos);

        // Why isn't this at the beginning of this function?
        // assert(length < 259);

        self.litlens.push(litlen);
        match litlen {
            LitLen::Literal(length) => {
                self.ll_symbol.push(length);
                self.d_symbol.push(0);
                self.ll_counts[llstart + length as usize] += 1;
            },
            LitLen::LengthDist(length, dist) => {
                let len_sym = get_length_symbol(length as usize);
                self.ll_symbol.push(len_sym as u16);
                self.d_symbol.push(get_dist_symbol(dist as i32) as u16);
                self.ll_counts[llstart + len_sym as usize] += 1;
                self.d_counts[dstart + get_dist_symbol(dist as i32) as usize] += 1;
            },
        }
    }

    pub fn lit_len_dist(&mut self, length: u16, dist: u16, pos: usize) {
        let litlen = if dist == 0 {
            LitLen::Literal(length)
        } else {
            LitLen::LengthDist(length, dist)
        };

        self.append_store_item(litlen, pos);
    }

    /// Does LZ77 using an algorithm similar to gzip, with lazy matching, rather than
    /// with the slow but better "squeeze" implementation.
    /// The result is placed in the Lz77Store.
    /// If instart is larger than 0, it uses values before instart as starting
    /// dictionary.
    pub fn greedy<C>(&mut self, s: &mut ZopfliBlockState<C>, in_data: &[u8], instart: usize, inend: usize)
        where C: Cache,
    {
        if instart == inend {
            return;
        }
        let windowstart = instart.saturating_sub(ZOPFLI_WINDOW_SIZE);
        let mut h = ZopfliHash::new();

        let arr = &in_data[..inend];
        h.warmup(arr, windowstart, inend);

        for i in windowstart..instart {
            h.update(arr, i);
        }

        let mut i = instart;
        let mut leng;
        let mut dist;
        let mut lengthscore;

        /* Lazy matching. */
        let mut prev_length = 0;
        let mut prev_match = 0;
        let mut prevlengthscore;
        let mut match_available = false;
        while i < inend {
            h.update(arr, i);

            let longest_match = find_longest_match(s, &mut h, arr, i, inend, ZOPFLI_MAX_MATCH, &mut None);
            dist = longest_match.distance;
            leng = longest_match.length;
            lengthscore = get_length_score(leng as i32, dist as i32);

            /* Lazy matching. */
            prevlengthscore = get_length_score(prev_length as i32, prev_match as i32);
            if match_available {
                match_available = false;
                if lengthscore > prevlengthscore + 1 {
                    self.lit_len_dist(arr[i - 1] as u16, 0, i - 1);
                    if (lengthscore as usize) >= ZOPFLI_MIN_MATCH && (leng as usize) < ZOPFLI_MAX_MATCH {
                        match_available = true;
                        prev_length = leng as u32;
                        prev_match = dist as u32;
                        i += 1;
                        continue;
                    }
                } else {
                    /* Add previous to output. */
                    leng = prev_length as u16;
                    dist = prev_match as u16;
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
            } else if (lengthscore as usize) >= ZOPFLI_MIN_MATCH && (leng as usize) < ZOPFLI_MAX_MATCH {
                match_available = true;
                prev_length = leng as u32;
                prev_match = dist as u32;
                i += 1;
                continue;
            }
            /* End of lazy matching. */

            /* Add to output. */
            if (lengthscore as usize) >= ZOPFLI_MIN_MATCH {
                verify_len_dist(arr, i, dist, leng);
                self.lit_len_dist(leng, dist, i);
            } else {
                leng = 1;
                self.lit_len_dist(arr[i] as u16, 0, i);
            }
            for _ in 1..leng {
                assert!(i < inend);
                i += 1;
                h.update(arr, i);
            }
            i += 1;
        }
    }

    pub fn follow_path<C>(&mut self, in_data: &[u8], instart: usize, inend: usize, path: Vec<u16>, s: &mut ZopfliBlockState<C>)
        where C: Cache,
    {
        let windowstart = instart.saturating_sub(ZOPFLI_WINDOW_SIZE);

        if instart == inend {
            return;
        }

        let mut h = ZopfliHash::new();

        let arr = &in_data[..inend];
        h.warmup(arr, windowstart, inend);

        for i in windowstart..instart {
            h.update(arr, i);
        }

        let mut pos = instart;
        for &item in &path {
            let mut length = item;
            assert!(pos < inend);

            h.update(arr, pos);

            // Add to output.
            if length >= ZOPFLI_MIN_MATCH as u16 {
                // Get the distance by recalculating longest match. The found length
                // should match the length from the path.
                let longest_match = find_longest_match(s, &mut h, arr, pos, inend, length as usize, &mut None);
                let dist = longest_match.distance;
                let dummy_length = longest_match.length;
                assert!(!(dummy_length != length && length > 2 && dummy_length > 2));
                verify_len_dist(arr, pos, dist, length);
                self.lit_len_dist(length, dist, pos);
            } else {
                length = 1;
                self.lit_len_dist(arr[pos] as u16, 0, pos);
            }

            assert!(pos + (length as usize) <= inend);
            for j in 1..(length as usize) {
                h.update(arr, pos + j);
            }

            pos += length as usize;
        }
    }

    fn get_histogram_at(&self, lpos: usize) -> (Vec<usize>, Vec<usize>) {
        let mut ll = vec![0; ZOPFLI_NUM_LL];
        let mut d = vec![0; ZOPFLI_NUM_D];

        /* The real histogram is created by using the histogram for this chunk, but
        all superfluous values of this chunk subtracted. */
        let llpos = ZOPFLI_NUM_LL * (lpos / ZOPFLI_NUM_LL);
        let dpos = ZOPFLI_NUM_D * (lpos / ZOPFLI_NUM_D);

        for (i, item) in ll.iter_mut().enumerate() {
            *item = self.ll_counts[llpos + i];
        }
        let end = cmp::min(llpos + ZOPFLI_NUM_LL, self.size());
        for i in (lpos + 1)..end {
            ll[self.ll_symbol[i] as usize] -= 1;
        }

        for (i, item) in d.iter_mut().enumerate() {
            *item = self.d_counts[dpos + i];
        }
        let end = cmp::min(dpos + ZOPFLI_NUM_D, self.size());
        for i in (lpos + 1)..end {
            match self.litlens[i] {
                LitLen::LengthDist(_, _) => d[self.d_symbol[i] as usize] -= 1,
                _ => {},
            }
        }

        (ll, d)
    }

    /// Gets the histogram of lit/len and dist symbols in the given range, using the
    /// cumulative histograms, so faster than adding one by one for large range. Does
    /// not add the one end symbol of value 256.
    pub fn get_histogram(&self, lstart: usize, lend: usize) -> (Vec<usize>, Vec<usize>) {
        if lstart + ZOPFLI_NUM_LL * 3 > lend {
            let mut ll_counts = vec![0; ZOPFLI_NUM_LL];
            let mut d_counts = vec![0; ZOPFLI_NUM_D];
            for i in lstart..lend  {
                ll_counts[self.ll_symbol[i] as usize] += 1;
                match self.litlens[i] {
                    LitLen::LengthDist(_, _) => d_counts[self.d_symbol[i] as usize] += 1,
                    _ => {},
                }
            }
            (ll_counts, d_counts)
        } else {
            /* Subtract the cumulative histograms at the end and the start to get the
            histogram for this range. */
            let (ll, d) = self.get_histogram_at(lend - 1);

            if lstart > 0 {
                let (ll2, d2) = self.get_histogram_at(lstart - 1);

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

    pub fn get_byte_range(&self, lstart: usize, lend: usize) -> usize {
        if lstart == lend {
            return 0;
        }

        let l = lend - 1;
        self.pos[l] + self.litlens[l].size() - self.pos[lstart]
    }
}

/// Some state information for compressing a block.
/// This is currently a bit under-used (with mainly only the longest match cache),
/// but is kept for easy future expansion.
pub struct ZopfliBlockState<'a, C> {
    pub options: &'a Options,
    /* Cache for length/distance pairs found so far. */
    lmc: C,
    /* The start (inclusive) and end (not inclusive) of the current block. */
    pub blockstart: usize,
    pub blockend: usize,
}

impl<'a> ZopfliBlockState<'a, ZopfliLongestMatchCache> {
    pub fn new(options: &'a Options, blockstart: usize, blockend: usize) -> Self {
        ZopfliBlockState {
            options: options,
            blockstart: blockstart,
            blockend: blockend,
            lmc: ZopfliLongestMatchCache::new(blockend - blockstart),
        }
    }
}

impl<'a> ZopfliBlockState<'a, NoCache> {
    pub fn new_without_cache(options: &'a Options, blockstart: usize, blockend: usize) -> Self {
        ZopfliBlockState {
            options: options,
            blockstart: blockstart,
            blockend: blockend,
            lmc: NoCache,
        }
    }
}

impl<'a, C> ZopfliBlockState<'a, C>
    where C: Cache,
{
    /// Gets distance, length and sublen values from the cache if possible.
    /// Returns 1 if it got the values from the cache, 0 if not.
    /// Updates the limit value to a smaller one if possible with more limited
    /// information from the cache.
    fn try_get_from_longest_match_cache(&self, pos: usize, limit: usize, sublen: &mut Option<&mut [u16]>) -> LongestMatch {
        self.lmc.try_get(pos, limit, sublen, self.blockstart)
    }

    /// Stores the found sublen, distance and length in the longest match cache, if
    /// possible.
    fn store_in_longest_match_cache(&mut self, pos: usize, limit: usize, sublen: &mut Option<&mut [u16]>, distance: u16, length: u16) {
        self.lmc.store(pos, limit, sublen, distance, length, self.blockstart)
    }
}

pub struct LongestMatch {
    pub distance: u16,
    pub length: u16,
    pub from_cache: bool,
    pub limit: usize,
}

impl LongestMatch {
    pub fn new(limit: usize) -> Self {
        LongestMatch {
            distance: 0,
            length: 0,
            from_cache: false,
            limit: limit,
        }
    }
}

/// Finds how long the match of `scan` and `match` is. Can be used to find how many
/// bytes starting from `scan`, and from `match`, are equal. Returns the last byte
/// after `scan`, which is still equal to the corresponding byte after `match`.
/// `scan` is the position to compare; `match` is the earlier position to compare.
/// `end` is the last possible byte, beyond which to stop looking.
/// `safe_end` is a few (8) bytes before end, for comparing multiple bytes at once.
fn get_match(array: &[u8], scan_offset: usize, match_offset: usize, end: usize) -> usize {
    let mut scan_offset = scan_offset;
    let mut match_offset = match_offset;
    // /* 8 checks at once per array bounds check (usize is 64-bit). */
    // // C code has other options if usize is not 64-bit, but this is all I'm supporting
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

pub fn find_longest_match<C>(s: &mut ZopfliBlockState<C>, h: &mut ZopfliHash, array: &[u8], pos: usize, size: usize, limit: usize, mut sublen: &mut Option<&mut [u16]>) -> LongestMatch
    where C: Cache,
{
    let mut longest_match = s.try_get_from_longest_match_cache(pos, limit, sublen);

    if longest_match.from_cache {
        assert!(pos + (longest_match.length as usize) <= size);
        return longest_match;
    }

    let mut limit = longest_match.limit;

    assert!(limit <= ZOPFLI_MAX_MATCH);
    assert!(limit >= ZOPFLI_MIN_MATCH);
    assert!(pos < size);

    if size - pos < ZOPFLI_MIN_MATCH {
        /* The rest of the code assumes there are at least ZOPFLI_MIN_MATCH bytes to
        try. */
        longest_match.distance = 0;
        longest_match.length = 0;
        longest_match.from_cache = false;
        longest_match.limit = 0;
        return longest_match;
    }

    if pos + limit > size {
        limit = size - pos;
    }

    let (bestdist, bestlength) = find_longest_match_loop(h, array, pos, size, limit, sublen);

    s.store_in_longest_match_cache(pos, limit, sublen, bestdist as u16, bestlength as u16);

    assert!(bestlength <= limit);

    assert!(pos + bestlength <= size);
    longest_match.distance = bestdist as u16;
    longest_match.length = bestlength as u16;
    longest_match.from_cache = false;
    longest_match.limit = limit;
    longest_match
}

fn find_longest_match_loop(h: &mut ZopfliHash, array: &[u8], pos: usize, size: usize, limit: usize, sublen: &mut Option<&mut [u16]>) -> (i32, usize) {
    let mut which_hash = Which::Hash1;
    assert!(h.val(which_hash) < 65536);
    let mut pp = h.head_at(h.val(which_hash) as usize, which_hash);  /* During the whole loop, p == hprev[pp]. */
    let mut p = h.prev_at(pp as usize, which_hash);

    let hpos = pos & ZOPFLI_WINDOW_MASK;
    assert!(pp as usize == hpos);

    let mut dist = if (p as i32) < pp {
        pp - (p as i32)
    } else {
        (ZOPFLI_WINDOW_SIZE - (p as usize)) as i32 + pp
    };

    let mut bestlength = 1;
    let mut bestdist = 0;
    let mut chain_counter = ZOPFLI_MAX_CHAIN_HITS;  /* For quitting early. */
    let arrayend = pos + limit;
    let mut scan_offset;
    let mut match_offset;

    /* Go through all distances. */
    while (dist as usize) < ZOPFLI_WINDOW_SIZE {
        let mut currentlength = 0;

        assert!((p as usize) < ZOPFLI_WINDOW_SIZE);
        assert!(p == h.prev_at(pp as usize, which_hash));
        assert!(h.hash_val_at(p as usize, which_hash) == h.val(which_hash));

        if dist > 0 {
            assert!(pos < size);
            assert!((dist as usize) <= pos);
            scan_offset = pos;
            match_offset = pos - (dist as usize);

            /* Testing the byte at position bestlength first, goes slightly faster. */
            if pos + bestlength >= size || array[scan_offset + bestlength] == array[match_offset + bestlength] {

                let same0 = h.same[pos & ZOPFLI_WINDOW_MASK];
                if same0 > 2 && array[scan_offset] == array[match_offset] {
                    let same1 = h.same[(pos - (dist as usize)) & ZOPFLI_WINDOW_MASK];
                    let same = [same0, same1, limit as u16].iter().min().map(|x| *x as usize).unwrap();
                    scan_offset += same;
                    match_offset += same;
                }
                scan_offset = get_match(array, scan_offset, match_offset, arrayend);
                currentlength = scan_offset - pos;  /* The found length. */
            }

            if currentlength > bestlength {
                if let &mut Some(ref mut subl) = sublen {
                    for j in (bestlength + 1)..(currentlength + 1) {
                        subl[j] = dist as u16;
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
        if which_hash == Which::Hash1 && bestlength >= h.same[hpos] as usize && h.val(Which::Hash2) == h.hash_val_at(p as usize, Which::Hash2) {
            /* Now use the hash that encodes the length and first byte. */
            which_hash = Which::Hash2;
        }

        pp = p as i32;
        p = h.prev_at(p as usize, which_hash);
        if (p as i32) == pp {
            break;  /* Uninited prev value. */
        }

        dist += if (p as i32) < pp {
            pp - (p as i32)
        } else {
            (ZOPFLI_WINDOW_SIZE - (p as usize)) as i32 + pp
        };

        chain_counter -= 1;
        if chain_counter == 0 {
            break;
        }
    }
    (bestdist, bestlength)
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
fn get_length_score(length: i32, distance: i32) -> i32 {
    // At 1024, the distance uses 9+ extra bits and this seems to be the sweet spot
    // on tested files.
    if distance > 1024 {
        length - 1
    } else {
        length
    }
}

fn verify_len_dist(data: &[u8], pos: usize, dist: u16, length: u16) {
    for i in 0..length {
        let d1 = data[pos - (dist as usize) + (i as usize)];
        let d2 = data[pos + (i as usize)];
        if d1 != d2 {
            assert!(d1 == d2);
            break;
        }
    }
}
