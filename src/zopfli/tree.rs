use std::slice;

use libc::{c_double, size_t};

const K_INV_LOG2: c_double = 1.4426950408889;  // 1.0 / log(2.0)

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCalculateEntropy(count_ptr: *const size_t, n: size_t, bitlengths_ptr: *mut c_double) {
    let count = unsafe {
        assert!(!count_ptr.is_null());
        slice::from_raw_parts(count_ptr, n as usize)
    };
    let bitlengths = unsafe {
        assert!(!bitlengths_ptr.is_null());
        slice::from_raw_parts_mut(bitlengths_ptr, n as usize)
    };

    let mut sum = 0;
    for i in 0..(n as usize) {
        sum += count[i];
    }

    let log2sum = (if sum == 0 { n } else { sum } as c_double).ln() * K_INV_LOG2;

    for i in 0..n as usize {
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
