use std::slice;

use libc::{c_uint, c_int, size_t};

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
