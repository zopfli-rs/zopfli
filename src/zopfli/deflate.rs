use libc::{c_uint, c_int, size_t, c_uchar, c_double};

use blocksplitter::{blocksplit, blocksplit_lz77};
use katajainen::length_limited_code_lengths;
use lz77::{get_histogram, get_byte_range, ZopfliBlockState, Lz77Store};
use squeeze::{lz77_optimal_fixed, lz77_optimal};
use symbols::{get_length_symbol, get_dist_symbol, get_length_symbol_extra_bits, get_dist_symbol_extra_bits, get_length_extra_bits_value, get_length_extra_bits, get_dist_extra_bits_value, get_dist_extra_bits};
use tree::{lengths_to_symbols};
use util::{ZOPFLI_NUM_LL, ZOPFLI_NUM_D, ZOPFLI_MASTER_BLOCK_SIZE};
use zopfli::ZopfliOptions;

pub fn fixed_tree() -> (Vec<c_uint>, Vec<c_uint>) {
    let mut ll = vec![8; ZOPFLI_NUM_LL];
    for i in 144..256 {
        ll[i] = 9;
    }
    for i in 256..280 {
        ll[i] = 7;
    }
    let d = vec![5; ZOPFLI_NUM_D];
    (ll, d)
}

/// Changes the population counts in a way that the consequent Huffman tree
/// compression, especially its rle-part, will be more likely to compress this data
/// more efficiently. length contains the size of the histogram.
pub fn optimize_huffman_for_rle(counts: &mut [size_t]) {
    let mut length = counts.len();
    // 1) We don't want to touch the trailing zeros. We may break the
    // rules of the format by adding more data in the distance codes.
    loop {
        if length == 0 {
            return;
        }
        if counts[length - 1] != 0 {
            // Now counts[0..length - 1] does not have trailing zeros.
            break;
        }
        length -= 1;
    }

    // 2) Let's mark all population counts that already can be encoded
    // with an rle code.
    let mut good_for_rle = vec![0; length];

    // Let's not spoil any of the existing good rle codes.
    // Mark any seq of 0's that is longer than 5 as a good_for_rle.
    // Mark any seq of non-0's that is longer than 7 as a good_for_rle.
    let mut symbol = counts[0];
    let mut stride = 0;
    for i in 0..(length + 1) {
        if i == length || counts[i] != symbol {
            if (symbol == 0 && stride >= 5) || (symbol != 0 && stride >= 7) {
                for k in 0..stride {
                    good_for_rle[i - k - 1] = 1;
                }
            }
            stride = 1;
            if i != length {
                symbol = counts[i];
            }
        } else {
            stride += 1;
        }
    }

    // 3) Let's replace those population counts that lead to more rle codes.
    stride = 0;
    let mut limit = counts[0];
    let mut sum = 0;
    for i in 0..(length + 1) {
        // Heuristic for selecting the stride ranges to collapse.
        if i == length || good_for_rle[i] != 0 || (counts[i] as i32 - limit as i32).abs() >= 4 {
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
                    counts[i - k - 1] = count;
                }
            }
            stride = 0;
            sum = 0;
            if length > 2 && i < length - 3 {
                // All interesting strides have a count of at least 4,
                // at least when non-zeros.
                limit = (counts[i] + counts[i + 1] + counts[i + 2] + counts[i + 3] + 2) / 4;
            } else if i < length {
                limit = counts[i];
            } else {
                limit = 0;
            }
        }
        stride += 1;
        if i != length {
            sum += counts[i];
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
pub fn patch_distance_codes_for_buggy_decoders(d_lengths: &mut[c_uint]) {
    let mut num_dist_codes = 0; // Amount of non-zero distance codes
    // Ignore the two unused codes from the spec
    for i in 0..30 {
        if d_lengths[i] != 0 {
            num_dist_codes += 1;
        }
        // Two or more codes is fine.
        if num_dist_codes >= 2 {
            return;
        }
    }

    if num_dist_codes == 0 {
        d_lengths[0] = 1;
        d_lengths[1] = 1;
    } else if num_dist_codes == 1 {
        let index = if d_lengths[0] == 0 {
            0
        } else {
            1
        };
        d_lengths[index] = 1;
    }
}

/// Same as calculate_block_symbol_size, but for block size smaller than histogram
/// size.
pub fn calculate_block_symbol_size_small(ll_lengths: &[c_uint], d_lengths: &[c_uint], lz77: &Lz77Store, lstart: size_t, lend: size_t) -> size_t {
    let mut result = 0;

    for i in lstart..lend {
        assert!(i < lz77.size());
        let litlens_i = lz77.litlens[i];
        let dists_i = lz77.dists[i];
        assert!(litlens_i < 259);
        if dists_i == 0 {
            result += ll_lengths[litlens_i as usize];
        } else {
            let ll_symbol = get_length_symbol(litlens_i as c_int);
            let d_symbol = get_dist_symbol(dists_i as c_int);
            result += ll_lengths[ll_symbol as usize];
            result += d_lengths[d_symbol as usize];
            result += get_length_symbol_extra_bits(ll_symbol) as c_uint;
            result += get_dist_symbol_extra_bits(d_symbol) as c_uint;
        }
    }
    result += ll_lengths[256]; // end symbol
    result as size_t
}

/// Same as calculate_block_symbol_size, but with the histogram provided by the caller.
pub fn calculate_block_symbol_size_given_counts(ll_counts: &Vec<size_t>, d_counts: &Vec<size_t>, ll_lengths: &[c_uint], d_lengths: &[c_uint], lz77: &Lz77Store, lstart: size_t, lend: size_t) -> size_t {
    let mut result = 0;

    if lstart + ZOPFLI_NUM_LL * 3 > lend {
        return calculate_block_symbol_size_small(ll_lengths, d_lengths, lz77, lstart, lend);
    } else {
        for i in 0..256 {
            result += ll_lengths[i] * ll_counts[i] as c_uint;
        }
        for i in 257..286 {
            result += ll_lengths[i] * ll_counts[i] as c_uint;
            result += (get_length_symbol_extra_bits(i as c_int) * ll_counts[i] as c_int) as c_uint;
        }
        for i in 0..30 {
            result += d_lengths[i] * d_counts[i] as c_uint;
            result += (get_dist_symbol_extra_bits(i as c_int) * d_counts[i] as c_int) as c_uint;
        }
        result += ll_lengths[256]; // end symbol
        result as size_t
    }
}

/// Calculates size of the part after the header and tree of an LZ77 block, in bits.
pub fn calculate_block_symbol_size(ll_lengths: &[c_uint], d_lengths: &[c_uint], lz77: &Lz77Store, lstart: size_t, lend: size_t) -> size_t {
    if lstart + ZOPFLI_NUM_LL * 3 > lend {
        calculate_block_symbol_size_small(ll_lengths, d_lengths, lz77, lstart, lend)
    } else {
        let (ll_counts, d_counts) = get_histogram(lz77, lstart, lend);
        calculate_block_symbol_size_given_counts(&ll_counts, &d_counts, ll_lengths, d_lengths, lz77, lstart, lend)
    }
}

/// Encodes the Huffman tree and returns how many bits its encoding takes; only returns the size
/// and runs faster.
pub fn encode_tree_no_output(ll_lengths: &[c_uint], d_lengths: &[c_uint], use_16: bool, use_17: bool, use_18: bool) -> size_t {
    let mut hlit = 29;  /* 286 - 257 */
    let mut hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/

    let mut clcounts = [0; 19];
    /* The order in which code length code lengths are encoded as per deflate. */
    let order = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    ];
    let mut result_size = 0;

    /* Trim zeros. */
    while hlit > 0 && ll_lengths[257 + hlit - 1] == 0 {
        hlit -= 1;
    }
    while hdist > 0 && d_lengths[1 + hdist - 1] == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;

    let lld_total = hlit2 + hdist + 1; /* Total amount of literal, length, distance codes. */

    let mut i = 0;

    while i < lld_total {
        /* This is an encoding of a huffman tree, so now the length is a symbol */
        let symbol = if i < hlit2 {
            ll_lengths[i]
        } else {
            d_lengths[i - hlit2]
        } as c_uchar;

        let mut count = 1;
        if use_16 || (symbol == 0 && (use_17 || use_18)) {
            let mut j = i + 1;
            let mut symbol_calc = if j < hlit2 {
                ll_lengths[j]
            } else {
                d_lengths[j - hlit2]
            } as c_uchar;

            while j < lld_total && symbol == symbol_calc {
                count += 1;
                j += 1;
                symbol_calc = if j < hlit2 {
                    ll_lengths[j]
                } else {
                    d_lengths[j - hlit2]
                } as c_uchar;
            }
        }

        i += count - 1;

        /* Repetitions of zeroes */
        if symbol == 0 && count >= 3 {
            if use_18 {
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
            if use_17 {
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
        if use_16 && count >= 4 {
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
        clcounts[symbol as usize] += count;
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
        result_size += clcl[i] as size_t * clcounts[i];
    }
    /* Extra bits. */
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    result_size
}

/// Gives the exact size of the tree, in bits, as it will be encoded in DEFLATE.
pub fn calculate_tree_size(ll_lengths: &[c_uint], d_lengths: &[c_uint]) -> size_t {
    let mut result = 0;
    for i in 0..8 {
        let size = encode_tree_no_output(ll_lengths, d_lengths, i & 1 > 0, i & 2 > 0, i & 4 > 0);
        if result == 0 || size < result {
            result = size;
        }
    }
    result
}

/// Encodes the Huffman tree and returns how many bits its encoding takes and returns output.
pub fn encode_tree(ll_lengths: &[c_uint], d_lengths: &[c_uint], use_16: bool, use_17: bool, use_18: bool, bp: *mut c_uchar, out: &mut Vec<u8>) -> size_t {
    let mut hlit = 29;  /* 286 - 257 */
    let mut hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/

    let mut clcounts = [0; 19];
    /* The order in which code length code lengths are encoded as per deflate. */
    let order = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    ];
    let mut result_size = 0;


    let mut rle = vec![];
    let mut rle_bits = vec![];

    /* Trim zeros. */
    while hlit > 0 && ll_lengths[257 + hlit - 1] == 0 {
        hlit -= 1;
    }
    while hdist > 0 && d_lengths[1 + hdist - 1] == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;

    let lld_total = hlit2 + hdist + 1; /* Total amount of literal, length, distance codes. */

    let mut i = 0;

    while i < lld_total {
        /* This is an encoding of a huffman tree, so now the length is a symbol */
        let symbol = if i < hlit2 {
            ll_lengths[i]
        } else {
            d_lengths[i - hlit2]
        } as c_uchar;

        let mut count = 1;
        if use_16 || (symbol == 0 && (use_17 || use_18 )) {
            let mut j = i + 1;
            let mut symbol_calc = if j < hlit2 {
                ll_lengths[j]
            } else {
                d_lengths[j - hlit2]
            } as c_uchar;

            while j < lld_total && symbol == symbol_calc {
                count += 1;
                j += 1;
                symbol_calc = if j < hlit2 {
                    ll_lengths[j]
                } else {
                    d_lengths[j - hlit2]
                } as c_uchar;
            }
        }

        i += count - 1;

        /* Repetitions of zeroes */
        if symbol == 0 && count >= 3 {
            if use_18 {
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
            if use_17 {
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
        if use_16 && count >= 4 {
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
        clcounts[symbol as usize] += count;
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

    add_bits(hlit as c_uint, 5, bp, out);
    add_bits(hdist as c_uint, 5, bp, out);
    add_bits(hclen as c_uint, 4, bp, out);

    for i in 0..(hclen + 4) {
        add_bits(clcl[order[i]], 3, bp, out);
    }

    for i in 0..rle.len() {
        let rle_i = rle[i] as usize;
        let rle_bits_i = rle_bits[i] as c_uint;
        let sym = clsymbols[rle_i];
        add_huffman_bits(sym, clcl[rle_i], bp, out);
        /* Extra bits. */
        if rle_i == 16 {
            add_bits(rle_bits_i, 2, bp, out);
        } else if rle_i == 17 {
            add_bits(rle_bits_i, 3, bp, out);
        } else if rle_i == 18 {
            add_bits(rle_bits_i, 7, bp, out);
        }
    }

    result_size += 14;  /* hlit, hdist, hclen bits */
    result_size += (hclen + 4) * 3;  /* clcl bits */
    for i in 0..19 {
        result_size += clcl[i] as size_t * clcounts[i];
    }
    /* Extra bits. */
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    result_size
}

pub fn add_dynamic_tree(ll_lengths: &[c_uint], d_lengths: &[c_uint], bp: *mut c_uchar, out: &mut Vec<u8>) {
    let mut best = 0;
    let mut bestsize = 0;

    for i in 0..8 {
        let size = encode_tree_no_output(ll_lengths, d_lengths, i & 1 > 0, i & 2 > 0, i & 4 > 0);
        if bestsize == 0 || size < bestsize {
            bestsize = size;
            best = i;
        }
    }

    encode_tree(ll_lengths, d_lengths, best & 1 > 0, best & 2 > 0, best & 4 > 0, bp, out);
}

/// Adds a deflate block with the given LZ77 data to the output.
/// options: global program options
/// btype: the block type, must be 1 or 2
/// final: whether to set the "final" bit on this block, must be the last block
/// litlens: literal/length array of the LZ77 data, in the same format as in
///     Lz77Store.
/// dists: distance array of the LZ77 data, in the same format as in
///     Lz77Store.
/// lstart: where to start in the LZ77 data
/// lend: where to end in the LZ77 data (not inclusive)
/// expected_data_size: the uncompressed block size, used for assert, but you can
///   set it to 0 to not do the assertion.
/// bp: output bit pointer
/// out: dynamic output array to append to
/// outsize: dynamic output array size
pub fn add_lz77_block(options: &ZopfliOptions, btype: c_int, final_block: c_int, in_data: &[u8], lz77: &Lz77Store, lstart: size_t, lend: size_t, expected_data_size: size_t, bp: *mut c_uchar, out: &mut Vec<u8>) {
    let compressed_size;
    let mut uncompressed_size: size_t = 0;

    if btype == 0 {
        let length = get_byte_range(lz77, lstart, lend);
        let pos = if lstart == lend {
            0
        } else {
            lz77.pos[lstart]
        };
        let end = pos + length;
        add_non_compressed_block(options, final_block, in_data, pos, end, bp, out);
        return;
    }

    add_bit(final_block, bp, out);
    add_bit(btype & 1, bp, out);
    add_bit((btype & 2) >> 1, bp, out);

    let (ll_lengths, d_lengths) = if btype == 1 {
        /* Fixed block. */
        fixed_tree()
    } else {
        /* Dynamic block. */
        assert!(btype == 2);
        let (_, ll_lengths, d_lengths) = get_dynamic_lengths(lz77, lstart, lend);

        let detect_tree_size = out.len();
        add_dynamic_tree(&ll_lengths, &d_lengths, bp, out);
        if options.verbose > 0 {
            println!("treesize: {}", out.len() - detect_tree_size);
        }
        (ll_lengths, d_lengths)
    };

    let ll_symbols = lengths_to_symbols(&ll_lengths, 15);
    let d_symbols = lengths_to_symbols(&d_lengths, 15);

    let detect_block_size = out.len();
    add_lz77_data(lz77, lstart, lend, expected_data_size, &ll_symbols, &ll_lengths, &d_symbols, &d_lengths, bp, out);

    /* End symbol. */
    add_huffman_bits(ll_symbols[256], ll_lengths[256], bp, out);

    for i in lstart..lend {
        uncompressed_size += if lz77.dists[i] == 0 {
            1
        } else {
            lz77.litlens[i] as size_t
        };
    }
    compressed_size = out.len() - detect_block_size;
    if options.verbose > 0 {
        println!("compressed block size: {} ({}k) (unc: {})", compressed_size, compressed_size / 1024, uncompressed_size);
    }
}

/// Calculates block size in bits.
/// litlens: lz77 lit/lengths
/// dists: ll77 distances
/// lstart: start of block
/// lend: end of block (not inclusive)
pub fn calculate_block_size(lz77: &Lz77Store, lstart: size_t, lend: size_t, btype: c_int) -> c_double {
    if btype == 0 {
        let length = get_byte_range(lz77, lstart, lend);
        let rem = length % 65535;
        let blocks = length / 65535 + (if rem > 0 { 1 } else { 0 });
        /* An uncompressed block must actually be split into multiple blocks if it's
           larger than 65535 bytes long. Eeach block header is 5 bytes: 3 bits,
           padding, LEN and NLEN (potential less padding for first one ignored). */
        (blocks * 5 * 8 + length * 8) as c_double
    } else if btype == 1 {
        let fixed_tree = fixed_tree();
        let ll_lengths = fixed_tree.0;
        let d_lengths = fixed_tree.1;

        let mut result: c_double = 3.0; /* bfinal and btype bits */
        result += calculate_block_symbol_size(&ll_lengths, &d_lengths, lz77, lstart, lend) as c_double;
        result
    } else {
        get_dynamic_lengths(lz77, lstart, lend).0 + 3.0
    }
}

/// Tries out OptimizeHuffmanForRle for this block, if the result is smaller,
/// uses it, otherwise keeps the original. Returns size of encoded tree and data in
/// bits, not including the 3-bit block header.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn try_optimize_huffman_for_rle(lz77: &Lz77Store, lstart: size_t, lend: size_t, ll_counts: Vec<size_t>, d_counts: Vec<size_t>, ll_lengths: Vec<c_uint>, d_lengths: Vec<c_uint>) -> (c_double, Vec<c_uint>, Vec<c_uint>) {
    let mut ll_counts2 = ll_counts.clone();
    let mut d_counts2 = d_counts.clone();

    let treesize = calculate_tree_size(&ll_lengths, &d_lengths);
    let datasize = calculate_block_symbol_size_given_counts(&ll_counts, &d_counts, &ll_lengths, &d_lengths, lz77, lstart, lend);

    optimize_huffman_for_rle(&mut ll_counts2);
    optimize_huffman_for_rle(&mut d_counts2);

    let ll_lengths2 = length_limited_code_lengths(&ll_counts2, 15);
    let mut d_lengths2 = length_limited_code_lengths(&d_counts2, 15);
    patch_distance_codes_for_buggy_decoders(&mut d_lengths2[..]);

    let treesize2 = calculate_tree_size(&ll_lengths2, &d_lengths2);
    let datasize2 = calculate_block_symbol_size_given_counts(&ll_counts, &d_counts, &ll_lengths2, &d_lengths2, lz77, lstart, lend);

    if treesize2 + datasize2 < treesize + datasize {
        (((treesize2 + datasize2) as c_double), ll_lengths2, d_lengths2)
    } else {
        ((treesize + datasize) as c_double, ll_lengths, d_lengths)
    }
}

/// Calculates the bit lengths for the symbols for dynamic blocks. Chooses bit
/// lengths that give the smallest size of tree encoding + encoding of all the
/// symbols to have smallest output size. This are not necessarily the ideal Huffman
/// bit lengths. Returns size of encoded tree and data in bits, not including the
/// 3-bit block header.
pub fn get_dynamic_lengths(lz77: &Lz77Store, lstart: size_t, lend: size_t) -> (c_double, Vec<c_uint>, Vec<c_uint>) {

    let (mut ll_counts, d_counts) = get_histogram(&lz77, lstart, lend);
    ll_counts[256] = 1;  /* End symbol. */

    let ll_lengths = length_limited_code_lengths(&ll_counts, 15);
    let mut d_lengths = length_limited_code_lengths(&d_counts, 15);

    patch_distance_codes_for_buggy_decoders(&mut d_lengths[..]);

    try_optimize_huffman_for_rle(lz77, lstart, lend, ll_counts, d_counts, ll_lengths, d_lengths)
}

/// Adds all lit/len and dist codes from the lists as huffman symbols. Does not add
/// end code 256. expected_data_size is the uncompressed block size, used for
/// assert, but you can set it to 0 to not do the assertion.
pub fn add_lz77_data(lz77: &Lz77Store, lstart: size_t, lend: size_t, expected_data_size: size_t , ll_symbols: &[c_uint], ll_lengths: &[c_uint], d_symbols: &[c_uint], d_lengths: &[c_uint], bp: *mut c_uchar, out: &mut Vec<u8>) {
    let mut testlength: size_t = 0;

    for i in lstart..lend {
        let dist: c_uint = lz77.dists[i] as c_uint;
        let litlen: c_uint = lz77.litlens[i] as c_uint;
        if dist == 0 {
            assert!(litlen < 256);
            assert!(ll_lengths[litlen as usize] > 0);
            add_huffman_bits(ll_symbols[litlen as usize], ll_lengths[litlen as usize], bp, out);
            testlength += 1;
        } else {
            let lls: c_uint = get_length_symbol(litlen as c_int) as c_uint;
            let ds: c_uint = get_dist_symbol(dist as c_int) as c_uint;
            assert!(litlen >= 3 && litlen <= 288);
            assert!(ll_lengths[lls as usize] > 0);
            assert!(d_lengths[ds as usize] > 0);
            add_huffman_bits(ll_symbols[lls as usize], ll_lengths[lls as usize], bp, out);
            add_bits(get_length_extra_bits_value(litlen as c_int) as c_uint, get_length_extra_bits(litlen as c_int) as c_uint, bp, out);
            add_huffman_bits(d_symbols[ds as usize], d_lengths[ds as usize], bp, out);
            add_bits(get_dist_extra_bits_value(dist as c_int) as c_uint, get_dist_extra_bits(dist as c_int) as c_uint, bp, out);
            testlength += litlen as size_t;
        }
    }
    assert!(expected_data_size == 0 || testlength == expected_data_size);
}

pub fn add_lz77_block_auto_type(options: &ZopfliOptions, final_block: c_int, in_data: &[u8], lz77: &Lz77Store, lstart: size_t, lend: size_t, expected_data_size: size_t, bp: *mut c_uchar, out: &mut Vec<u8>) {
    let uncompressedcost = calculate_block_size(lz77, lstart, lend, 0);
    let mut fixedcost = calculate_block_size(lz77, lstart, lend, 1);
    let dyncost = calculate_block_size(lz77, lstart, lend, 2);

    /* Whether to perform the expensive calculation of creating an optimal block
    with fixed huffman tree to check if smaller. Only do this for small blocks or
    blocks which already are pretty good with fixed huffman tree. */
    let expensivefixed = (lz77.size() < 1000) || fixedcost <= dyncost * 1.1;

    let mut fixedstore = Lz77Store::new();
    if lstart == lend {
        /* Smallest empty block is represented by fixed block */
        add_bits(final_block as c_uint, 1, bp, out);
        add_bits(1, 2, bp, out);  /* btype 01 */
        add_bits(0, 7, bp, out);  /* end symbol has code 0000000 */
        return;
    }
    if expensivefixed {
        /* Recalculate the LZ77 with lz77_optimal_fixed */
        let instart = lz77.pos[lstart];
        let inend = instart + get_byte_range(lz77, lstart, lend);

        let mut s = ZopfliBlockState::new(options, instart, inend, 1);
        lz77_optimal_fixed(&mut s, in_data, instart, inend, &mut fixedstore);
        fixedcost = calculate_block_size(&fixedstore, 0, fixedstore.size(), 1);
    }

    if uncompressedcost < fixedcost && uncompressedcost < dyncost {
        add_lz77_block(options, 0, final_block, in_data, lz77, lstart, lend, expected_data_size, bp, out);
    } else if fixedcost < dyncost {
        if expensivefixed {
            add_lz77_block(options, 1, final_block, in_data, &fixedstore, 0, fixedstore.size(), expected_data_size, bp, out);
        } else {
            add_lz77_block(options, 1, final_block, in_data, lz77, lstart, lend, expected_data_size, bp, out);
        }
    } else {
        add_lz77_block(options, 2, final_block, in_data, lz77, lstart, lend, expected_data_size, bp, out);
    }
}

/// Calculates block size in bits, automatically using the best btype.
pub fn calculate_block_size_auto_type(lz77: &Lz77Store, lstart: size_t, lend: size_t) -> c_double {
    let uncompressedcost = calculate_block_size(lz77, lstart, lend, 0);
    /* Don't do the expensive fixed cost calculation for larger blocks that are
     unlikely to use it. */
    let fixedcost = if lz77.size() > 1000 { uncompressedcost } else { calculate_block_size(lz77, lstart, lend, 1) };
    let dyncost = calculate_block_size(lz77, lstart, lend, 2);
    if uncompressedcost < fixedcost && uncompressedcost < dyncost {
        uncompressedcost
    } else if fixedcost < dyncost {
        fixedcost
    } else {
        dyncost
    }
}

pub fn add_all_blocks(splitpoints: &Vec<size_t>, lz77: &Lz77Store, options: &ZopfliOptions, final_block: c_int, in_data: &[u8], bp: *mut c_uchar, out: &mut Vec<u8>) {
    let npoints = splitpoints.len();
    for i in 0..(npoints + 1) {
        let start = if i == 0 { 0 } else { splitpoints[i - 1] };
        let end = if i == npoints { lz77.size() } else { splitpoints[i] };
        let final_block_i = if i == npoints && final_block > 0 { 1 } else { 0 };
        add_lz77_block_auto_type(options, final_block_i, in_data, lz77, start, end, 0, bp, out);
    }
}

pub fn blocksplit_attempt(options: &ZopfliOptions, final_block: c_int, in_data: &[u8], instart: size_t, inend: size_t, bp: *mut c_uchar, out: &mut Vec<u8>) {
    let mut totalcost = 0.0;
    let mut lz77 = Lz77Store::new();

    /* byte coordinates rather than lz77 index */
    let mut splitpoints_uncompressed = Vec::with_capacity(options.blocksplittingmax as usize);

    blocksplit(options, in_data, instart, inend, options.blocksplittingmax as usize, &mut splitpoints_uncompressed);
    let npoints = splitpoints_uncompressed.len();
    let mut splitpoints = Vec::with_capacity(npoints);

    for i in 0..(npoints + 1) {
        let start = if i == 0 { instart } else { splitpoints_uncompressed[i - 1] };
        let end = if i == npoints { inend } else { splitpoints_uncompressed[i] };
        let mut s = ZopfliBlockState::new(options, start, end, 1);

        let store = lz77_optimal(&mut s, in_data, start, end, options.numiterations);
        totalcost += calculate_block_size_auto_type(&store, 0, store.size());

        // ZopfliAppendLZ77Store(&store, &lz77);
        assert!(store.size() > 0);
        for j in 0..store.size() {
            lz77.lit_len_dist(store.litlens[j], store.dists[j], store.pos[j]);
        }

        if i < npoints {
            splitpoints.push(lz77.size());
        }
    }

    /* Second block splitting attempt */
    if npoints > 1 {
        let mut splitpoints2 = Vec::with_capacity(splitpoints_uncompressed.len());
        let mut totalcost2 = 0.0;

        blocksplit_lz77(options, &lz77, options.blocksplittingmax as usize, &mut splitpoints2);
        let npoints2 = splitpoints2.len();

        for i in 0..(npoints2 + 1) {
            let start = if i == 0 { 0 } else { splitpoints2[i - 1] };
            let end = if i == npoints2 { lz77.size() } else { splitpoints2[i] };
            totalcost2 += calculate_block_size_auto_type(&lz77, start, end);
        }

        if totalcost2 < totalcost {
            splitpoints = splitpoints2;
        }
    }

    add_all_blocks(&splitpoints, &lz77, &options, final_block, in_data, bp, out);
}


/// Deflate a part, to allow deflate() to use multiple master blocks if
/// needed.
/// It is possible to call this function multiple times in a row, shifting
/// instart and inend to next bytes of the data. If instart is larger than 0, then
/// previous bytes are used as the initial dictionary for LZ77.
/// This function will usually output multiple deflate blocks. If final is 1, then
/// the final bit will be set on the last block.
/// Like deflate, but allows to specify start and end byte with instart and
/// inend. Only that part is compressed, but earlier bytes are still used for the
/// back window.
pub fn deflate_part(options: &ZopfliOptions, btype: c_int, final_block: c_int, in_data: &[u8], instart: size_t, inend: size_t, bp: *mut c_uchar, out: &mut Vec<u8>) {
    /* If btype=2 is specified, it tries all block types. If a lesser btype is
    given, then however it forces that one. Neither of the lesser types needs
    block splitting as they have no dynamic huffman trees. */
    if btype == 0 {
        add_non_compressed_block(options, final_block, in_data, instart, inend, bp, out);
    } else if btype == 1 {
        let mut store = Lz77Store::new();
        let mut s = ZopfliBlockState::new(options, instart, inend, 1);

        lz77_optimal_fixed(&mut s, in_data, instart, inend, &mut store);
        add_lz77_block(options, btype, final_block, in_data, &store, 0, store.size(), 0, bp, out);
    } else {
        blocksplit_attempt(options, final_block, in_data, instart, inend, bp, out);
    }
}

/// Compresses according to the deflate specification and append the compressed
/// result to the output.
/// This function will usually output multiple deflate blocks. If final is 1, then
/// the final bit will be set on the last block.
///
/// options: global program options
/// btype: the deflate block type. Use 2 for best compression.
///   -0: non compressed blocks (00)
///   -1: blocks with fixed tree (01)
///   -2: blocks with dynamic tree (10)
/// final: whether this is the last section of the input, sets the final bit to the
///   last deflate block.
/// in: the input bytes
/// insize: number of input bytes
/// bp: bit pointer for the output array. This must initially be 0, and for
///   consecutive calls must be reused (it can have values from 0-7). This is
///   because deflate appends blocks as bit-based data, rather than on byte
///   boundaries.
/// out: pointer to the dynamic output array to which the result is appended. Must
///   be freed after use.
/// outsize: pointer to the dynamic output array size.
pub fn deflate(options_ptr: *const ZopfliOptions, btype: c_int, final_block: c_int, in_data: &[u8], bp: *mut c_uchar, out: &mut Vec<u8>) {
    let options = unsafe {
        assert!(!options_ptr.is_null());
        &*options_ptr
    };

    let offset = out.len();
    let mut i = 0;
    let insize = in_data.len();
    while i < insize {
        let masterfinal = i + ZOPFLI_MASTER_BLOCK_SIZE >= insize;
        let final2 = final_block != 0 && masterfinal;
        let size = if masterfinal { insize - i } else { ZOPFLI_MASTER_BLOCK_SIZE };
        let final2_as_int = if final2 { 1 } else { 0 };
        deflate_part(options, btype, final2_as_int, in_data, i, i + size, bp, out);
        i += size;
    }
    if options.verbose != 0 {
        let outsize = out.len();
        println!("Original Size: {}, Deflate: {}, Compression: {}% Removed", insize, outsize - offset, 100.0 * (insize - (outsize - offset)) as c_double / insize as c_double);
    }
}

/// bp = bitpointer, always in range [0, 7].
/// The outsize is number of necessary bytes to encode the bits.
/// Given the value of bp and the amount of bytes, the amount of bits represented
/// is not simply bytesize * 8 + bp because even representing one bit requires a
/// whole byte. It is: (bp == 0) ? (bytesize * 8) : ((bytesize - 1) * 8 + bp)
pub fn add_bit(bit: c_int, bp: *mut c_uchar, out: &mut Vec<u8>) {
    if unsafe { *bp } == 0 {
        out.push(0);
    }
    let outsize = out.len();
    out[outsize - 1] |= (bit << unsafe { *bp }) as u8;
    unsafe { *bp = (*bp + 1) & 7};
}

pub fn add_bits(symbol: c_uint, length: c_uint, bp: *mut c_uchar, out: &mut Vec<u8>) {
    /* TODO(lode): make more efficient (add more bits at once). */
    for i in 0..length {
        let bit = (symbol >> i) & 1;
        if unsafe { *bp } == 0 {
            out.push(0);
        }
        let outsize = out.len();
        out[outsize - 1] |= (bit << unsafe { *bp }) as u8;
        unsafe { *bp = (*bp + 1) & 7};
    }
}

/// Adds bits, like add_bits, but the order is inverted. The deflate specification
/// uses both orders in one standard.
pub fn add_huffman_bits(symbol: c_uint, length: c_uint, bp: *mut c_uchar, out: &mut Vec<u8>) {
    /* TODO(lode): make more efficient (add more bits at once). */
    for i in 0..length {
        let bit = (symbol >> (length - i - 1)) & 1;
        if unsafe { *bp } == 0 {
            out.push(0);
        }
        let outsize = out.len();
        out[outsize - 1] |= (bit << unsafe { *bp }) as u8;
        unsafe { *bp = (*bp + 1) & 7};
    }
}

/// Since an uncompressed block can be max 65535 in size, it actually adds
/// multible blocks if needed.
pub fn add_non_compressed_block(_options: &ZopfliOptions, final_block: c_int, in_data: &[u8], instart: size_t, inend: size_t, bp: *mut c_uchar, out: &mut Vec<u8>) {
    let mut pos = instart;
    loop {
        let mut blocksize = 65535;

        if pos + blocksize > inend {
            blocksize = inend - pos;
        }
        let currentfinal = if pos + blocksize >= inend { 1 } else { 0 };

        let nlen = !blocksize;

        add_bit(final_block & currentfinal, bp, out);
        /* BTYPE 00 */
        add_bit(0, bp, out);
        add_bit(0, bp, out);

        /* Any bits of input up to the next byte boundary are ignored. */
        unsafe { *bp = 0; }

        out.push((blocksize % 256) as u8);
        out.push(((blocksize / 256) % 256) as u8);
        out.push((nlen % 256) as u8);
        out.push(((nlen / 256) % 256) as u8);

        for i in 0..blocksize {
            out.push(in_data[pos + i]);
        }

        if currentfinal != 0 {
            break;
        }
        pos += blocksize;
    }
}
