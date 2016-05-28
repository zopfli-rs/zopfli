use std::slice;

use libc::{c_uint, c_int, size_t, c_uchar};

use katajainen::length_limited_code_lengths;
use lz77::{ZopfliLZ77Store, lz77_store_from_c, get_histogram};
use symbols::{ZopfliGetLengthSymbol, ZopfliGetDistSymbol, ZopfliGetLengthSymbolExtraBits, ZopfliGetDistSymbolExtraBits, ZOPFLI_NUM_LL};
use tree::lengths_to_symbols;

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn GetFixedTree(ll_lengths: *mut c_uint, d_lengths: *mut c_uint) {
    unsafe {
        for i in 0..144 {
            *ll_lengths.offset(i) = 8;
        }
        for i in 144..256 {
            *ll_lengths.offset(i) = 9;
        }
        for i in 256..280 {
            *ll_lengths.offset(i) = 7;
        }
        for i in 280..288 {
            *ll_lengths.offset(i) = 8;
        }
        for i in 0..32 {
            *d_lengths.offset(i) = 5;
        }
    }
}


/// Change the population counts in a way that the consequent Huffman tree
/// compression, especially its rle-part will be more likely to compress this data
/// more efficiently. length containts the size of the histogram.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn OptimizeHuffmanForRle(length: c_int, counts: *mut size_t) {
    let mut length = length as usize;
    let c = unsafe { slice::from_raw_parts(counts, length) };

    // 1) We don't want to touch the trailing zeros. We may break the
    // rules of the format by adding more data in the distance codes.
    loop {
        if length == 0 {
            return;
        }
        if c[length - 1] != 0 {
            // Now counts[0..length - 1] does not have trailing zeros.
            break;
        }
        length -= 1;
    }

    // 2) Let's mark all population counts that already can be encoded
    // with an rle code.
    let mut good_for_rle = vec![0; length as size_t];

    // Let's not spoil any of the existing good rle codes.
    // Mark any seq of 0's that is longer than 5 as a good_for_rle.
    // Mark any seq of non-0's that is longer than 7 as a good_for_rle.
    let mut symbol = c[0];
    let mut stride = 0;
    for i in 0..(length + 1) {
        if i == length || c[i] != symbol {
            if (symbol == 0 && stride >= 5) || (symbol != 0 && stride >= 7) {
                for k in 0..stride {
                    good_for_rle[(i - k - 1) as usize] = 1;
                }
            }
            stride = 1;
            if i != length {
                symbol = c[i];
            }
        } else {
            stride += 1;
        }
    }

    // 3) Let's replace those population counts that lead to more rle codes.
    stride = 0;
    let mut limit = c[0];
    let mut sum = 0;
    for i in 0..(length + 1) {
        // Heuristic for selecting the stride ranges to collapse.
        if i == length || good_for_rle[i as usize] != 0 || (c[i] as i32 - limit as i32).abs() >= 4 {
            if stride >= 4 || (stride >= 3 && sum == 0) {
                // The stride must end, collapse what we have, if we have enough (4).
                let mut count = (sum + stride / 2) / stride;
                if count < 1 {
                    count = 1;
                }
                if sum == 0 {
                    // Don't make an all zeros stride to be upgraded to ones.
                    count = 0;
                }
                for k in 0..stride {
                    // We don't want to change value at counts[i],
                    // that is already belonging to the next stride. Thus - 1.
                    unsafe { *counts.offset((i - k - 1) as isize) = count as size_t };
                }
            }
            stride = 0;
            sum = 0;
            if length > 2 && i < length - 3 {
                // All interesting strides have a count of at least 4,
                // at least when non-zeros.
                limit = (c[i] + c[i + 1] + c[i + 2] + c[i + 3] + 2) / 4;
            } else if i < length {
                limit = c[i];
            } else {
                limit = 0;
            }
        }
        stride += 1;
        if i != length {
            sum += c[i];
        }
    }
}

// Ensures there are at least 2 distance codes to support buggy decoders.
// Zlib 1.2.1 and below have a bug where it fails if there isn't at least 1
// distance code (with length > 0), even though it's valid according to the
// deflate spec to have 0 distance codes. On top of that, some mobile phones
// require at least two distance codes. To support these decoders too (but
// potentially at the cost of a few bytes), add dummy code lengths of 1.
// References to this bug can be found in the changelog of
// Zlib 1.2.2 and here: http://www.jonof.id.au/forum/index.php?topic=515.0.
//
// d_lengths: the 32 lengths of the distance codes.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn PatchDistanceCodesForBuggyDecoders(d_lengths: *mut c_uint) {
    let mut num_dist_codes = 0; // Amount of non-zero distance codes
    // Ignore the two unused codes from the spec
    for i in 0..30 {
        if unsafe { *d_lengths.offset(i as isize) } != 0 {
            num_dist_codes += 1;
        }
        // Two or more codes is fine.
        if num_dist_codes >= 2 {
            return;
        }
    }

    if num_dist_codes == 0 {
        unsafe {
            *d_lengths.offset(0) = 1;
            *d_lengths.offset(1) = 1;
        }
    } else if num_dist_codes == 1 {
        unsafe {
            let index = if *d_lengths.offset(0) == 0 {
                0
            } else {
                1
            };
            *d_lengths.offset(index) = 1;
        }
    }
}

/// Same as CalculateBlockSymbolSize, but for block size smaller than histogram
/// size.
pub fn calculate_block_symbol_size_small(ll_lengths: *const c_uint, d_lengths: *const c_uint, lz77_ptr: *const ZopfliLZ77Store, lstart: size_t, lend: size_t) -> size_t {
    let rust_store = lz77_store_from_c(lz77_ptr);
    let rs = unsafe { &*rust_store };
    let mut result = 0;

    for i in lstart..lend {
        assert!(i < rs.size());
        let litlens_i = rs.litlens[i];
        let dists_i = rs.dists[i];
        assert!(litlens_i < 259);
        if dists_i == 0 {
            result += unsafe { *ll_lengths.offset(litlens_i as isize) };
        } else {
            let ll_symbol = ZopfliGetLengthSymbol(litlens_i as c_int);
            let d_symbol = ZopfliGetDistSymbol(dists_i as c_int);
            result += unsafe { *ll_lengths.offset(ll_symbol as isize) };
            result += unsafe { *d_lengths.offset(d_symbol as isize) };
            result += ZopfliGetLengthSymbolExtraBits(ll_symbol) as c_uint;
            result += ZopfliGetDistSymbolExtraBits(d_symbol) as c_uint;
        }
    }
    result += unsafe { *ll_lengths.offset(256) }; // end symbol
    result as size_t
}

/// Same as CalculateBlockSymbolSize, but with the histogram provided by the caller.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CalculateBlockSymbolSizeGivenCounts(ll_counts: *const size_t, d_counts: *const size_t, ll_lengths: *const c_uint, d_lengths: *const c_uint, lz77: *const ZopfliLZ77Store, lstart: size_t, lend: size_t) -> size_t {
    let mut result = 0;

    if lstart + ZOPFLI_NUM_LL * 3 > lend {
        return calculate_block_symbol_size_small(ll_lengths, d_lengths, lz77, lstart, lend);
    } else {
        for i in 0..256 {
            result += unsafe { *ll_lengths.offset(i as isize) * *ll_counts.offset(i as isize) as c_uint };
        }
        for i in 257..286 {
            result += unsafe { *ll_lengths.offset(i as isize) * *ll_counts.offset(i as isize) as c_uint };
            result += (ZopfliGetLengthSymbolExtraBits(i) * unsafe { *ll_counts.offset(i as isize) as c_int }) as c_uint;
        }
        for i in 0..30 {
            result += unsafe { *d_lengths.offset(i as isize) * *d_counts.offset(i as isize) as c_uint };
            result += (ZopfliGetDistSymbolExtraBits(i) * unsafe { *d_counts.offset(i as isize) as c_int }) as c_uint;
        }
        result += unsafe { *ll_lengths.offset(256) }; // end symbol
        result as size_t
    }
}

/// Calculates size of the part after the header and tree of an LZ77 block, in bits.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CalculateBlockSymbolSize(ll_lengths: *const c_uint, d_lengths: *const c_uint, lz77: *const ZopfliLZ77Store, lstart: size_t, lend: size_t) -> size_t {
    if lstart + ZOPFLI_NUM_LL * 3 > lend {
        calculate_block_symbol_size_small(ll_lengths, d_lengths, lz77, lstart, lend)
    } else {
        let (ll_counts, d_counts) = get_histogram(unsafe { &*lz77 }, lstart, lend);
        CalculateBlockSymbolSizeGivenCounts(ll_counts.as_ptr(), d_counts.as_ptr(), ll_lengths, d_lengths, lz77, lstart, lend)
    }
}

/// Encodes the Huffman tree and returns how many bits its encoding takes; only returns the size
/// and runs faster.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn EncodeTreeNoOutput(ll_lengths: *const c_uint, d_lengths: *const c_uint, use_16: c_int, use_17: c_int, use_18: c_int) -> size_t {
    let mut hlit = 29;  /* 286 - 257 */
    let mut hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/

    let mut clcounts = [0 as usize; 19];
    /* The order in which code length code lengths are encoded as per deflate. */
    let order = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    ];
    let mut result_size = 0;

    /* Trim zeros. */
    while hlit > 0 && unsafe { *ll_lengths.offset(257 + hlit - 1) } == 0 {
        hlit -= 1;
    }
    while hdist > 0 && unsafe { *d_lengths.offset(1 + hdist - 1) } == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;

    let lld_total = hlit2 + hdist + 1; /* Total amount of literal, length, distance codes. */

    let mut i = 0;

    while i < lld_total {
        /* This is an encoding of a huffman tree, so now the length is a symbol */
        let symbol = if i < hlit2 {
            unsafe { *ll_lengths.offset(i) }
        } else {
            unsafe { *d_lengths.offset(i - hlit2) }
        } as c_uchar;

        let mut count = 1;
        if use_16 != 0 || (symbol == 0 && (use_17 != 0 || use_18 != 0)) {
            let mut j = i + 1;
            let mut symbol_calc = if j < hlit2 {
                unsafe { *ll_lengths.offset(j) }
            } else {
                unsafe { *d_lengths.offset(j - hlit2) }
            } as c_uchar;

            while j < lld_total && symbol == symbol_calc {
                count += 1;
                j += 1;
                symbol_calc = if j < hlit2 {
                    unsafe { *ll_lengths.offset(j) }
                } else {
                    unsafe { *d_lengths.offset(j - hlit2) }
                } as c_uchar;
            }
        }

        i += count - 1;

        /* Repetitions of zeroes */
        if symbol == 0 && count >= 3 {
            if use_18 != 0 {
                while count >= 11 {
                    let count2 = if count > 138 {
                        138
                    } else {
                        count
                    };
                    clcounts[18] += 1;
                    count -= count2;
                }
            }
            if use_17 != 0 {
                while count >= 3 {
                    let count2 = if count > 10 {
                        10
                    } else {
                        count
                    };
                    clcounts[17] += 1;
                    count -= count2;
                }
            }
        }

        /* Repetitions of any symbol */
        if use_16 != 0 && count >= 4 {
            count -= 1;  /* Since the first one is hardcoded. */
            clcounts[symbol as usize] += 1;
            while count >= 3 {
                let count2 = if count > 6 {
                    6
                } else {
                    count
                };
                clcounts[16] += 1;
                count -= count2;
            }
        }

        /* No or insufficient repetition */
        clcounts[symbol as usize] += count as usize;
        while count > 0 {
            count -= 1;
        }
        i += 1;
    }

    let clcl = length_limited_code_lengths(&clcounts, 7);

    let mut hclen = 15;
    /* Trim zeros. */
    while hclen > 0 && clcounts[order[hclen + 4 - 1]] == 0 {
        hclen -= 1;
    }

    result_size += 14;  /* hlit, hdist, hclen bits */
    result_size += (hclen + 4) * 3;  /* clcl bits */
    for i in 0..19 {
        result_size += clcl[i] * clcounts[i];
    }
    /* Extra bits. */
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    result_size
}

/// Gives the exact size of the tree, in bits, as it will be encoded in DEFLATE.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CalculateTreeSize(ll_lengths: *const c_uint, d_lengths: *const c_uint) -> size_t {
    let mut result = 0;
    for i in 0..8 {
        let size = EncodeTreeNoOutput(ll_lengths, d_lengths, i & 1, i & 2, i & 4);
        if result == 0 || size < result {
            result = size;
        }
    }
    result
}

#[link(name = "zopfli")]
extern {
    fn AddBits(symbol: c_uint, length: c_uint, bp: *const c_uchar, out: *const *const c_uint, outsize: *const size_t);
    fn AddHuffmanBits(symbol: c_uint, length: c_uint, bp: *const c_uchar, out: *const *const c_uint, outsize: *const size_t);
}

/// Encodes the Huffman tree and returns how many bits its encoding takes and returns output.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn EncodeTree(ll_lengths: *const c_uint, d_lengths: *const c_uint, use_16: c_int, use_17: c_int, use_18: c_int, bp: *const c_uchar, out: *const *const c_uint, outsize: *const size_t) -> size_t {
    let mut hlit = 29;  /* 286 - 257 */
    let mut hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/

    let mut clcounts = [0 as usize; 19];
    /* The order in which code length code lengths are encoded as per deflate. */
    let order = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    ];
    let mut result_size = 0;


    let mut rle = vec![];
    let mut rle_bits = vec![];

    /* Trim zeros. */
    while hlit > 0 && unsafe { *ll_lengths.offset(257 + hlit - 1) } == 0 {
        hlit -= 1;
    }
    while hdist > 0 && unsafe { *d_lengths.offset(1 + hdist - 1) } == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;

    let lld_total = hlit2 + hdist + 1; /* Total amount of literal, length, distance codes. */

    let mut i = 0;

    while i < lld_total {
        /* This is an encoding of a huffman tree, so now the length is a symbol */
        let symbol = if i < hlit2 {
            unsafe { *ll_lengths.offset(i) }
        } else {
            unsafe { *d_lengths.offset(i - hlit2) }
        } as c_uchar;

        let mut count = 1;
        if use_16 != 0 || (symbol == 0 && (use_17 != 0 || use_18 != 0)) {
            let mut j = i + 1;
            let mut symbol_calc = if j < hlit2 {
                unsafe { *ll_lengths.offset(j) }
            } else {
                unsafe { *d_lengths.offset(j - hlit2) }
            } as c_uchar;

            while j < lld_total && symbol == symbol_calc {
                count += 1;
                j += 1;
                symbol_calc = if j < hlit2 {
                    unsafe { *ll_lengths.offset(j) }
                } else {
                    unsafe { *d_lengths.offset(j - hlit2) }
                } as c_uchar;
            }
        }

        i += count - 1;

        /* Repetitions of zeroes */
        if symbol == 0 && count >= 3 {
            if use_18 != 0 {
                while count >= 11 {
                    let count2 = if count > 138 {
                        138
                    } else {
                        count
                    };
                    rle.push(18);
                    rle_bits.push(count2 - 11);
                    clcounts[18] += 1;
                    count -= count2;
                }
            }
            if use_17 != 0 {
                while count >= 3 {
                    let count2 = if count > 10 {
                        10
                    } else {
                        count
                    };
                    rle.push(17);
                    rle_bits.push(count2 - 3);
                    clcounts[17] += 1;
                    count -= count2;
                }
            }
        }

        /* Repetitions of any symbol */
        if use_16 != 0 && count >= 4 {
            count -= 1;  /* Since the first one is hardcoded. */
            clcounts[symbol as usize] += 1;
            rle.push(symbol);
            rle_bits.push(0);

            while count >= 3 {
                let count2 = if count > 6 {
                    6
                } else {
                    count
                };
                rle.push(16);
                rle_bits.push(count2 - 3);
                clcounts[16] += 1;
                count -= count2;
            }
        }

        /* No or insufficient repetition */
        clcounts[symbol as usize] += count as usize;
        while count > 0 {
            rle.push(symbol);
            rle_bits.push(0);
            count -= 1;
        }
        i += 1;
    }

    let clcl = length_limited_code_lengths(&clcounts, 7);
    let clsymbols = lengths_to_symbols(&clcl, 7);

    let mut hclen = 15;
    /* Trim zeros. */
    while hclen > 0 && clcounts[order[hclen + 4 - 1]] == 0 {
        hclen -= 1;
    }

    unsafe {
        AddBits(hlit as c_uint, 5, bp, out, outsize);
        AddBits(hdist as c_uint, 5, bp, out, outsize);
        AddBits(hclen as c_uint, 4, bp, out, outsize);

        for i in 0..(hclen + 4) {
            AddBits(clcl[order[i]] as c_uint, 3, bp, out, outsize);
        }

        for i in 0..rle.len() {
            let rle_i = rle[i] as usize;
            let rle_bits_i = rle_bits[i] as c_uint;
            let sym = clsymbols[rle_i];
            AddHuffmanBits(sym, clcl[rle_i] as c_uint, bp, out, outsize);
            /* Extra bits. */
            if rle_i == 16 {
                AddBits(rle_bits_i, 2, bp, out, outsize);
            } else if rle_i == 17 {
                AddBits(rle_bits_i, 3, bp, out, outsize);
            } else if rle_i == 18 {
                AddBits(rle_bits_i, 7, bp, out, outsize);
            }
        }
    }

    result_size += 14;  /* hlit, hdist, hclen bits */
    result_size += (hclen + 4) * 3;  /* clcl bits */
    for i in 0..19 {
        result_size += clcl[i] * clcounts[i];
    }
    /* Extra bits. */
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    result_size
}