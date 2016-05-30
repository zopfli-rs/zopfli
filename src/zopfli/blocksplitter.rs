use libc::{size_t, c_void, c_double, c_uchar, c_int};

use deflate::ZopfliCalculateBlockSizeAutoType;
use lz77::ZopfliLZ77Store;
use symbols::{ZOPFLI_LARGE_FLOAT};

/// Finds minimum of function f(i) where is is of type size_t, f(i) is of type
/// double, i is in range start-end (excluding end).
/// Outputs the minimum value in *smallest and returns the index of this value.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn FindMinimum(f: fn(i: size_t, context: *const c_void) -> c_double, context: *const c_void, start: size_t, end: size_t, smallest: *mut c_double) -> size_t {
    let mut start = start;
    let mut end = end;
    if end - start < 1024 {
        let mut best = ZOPFLI_LARGE_FLOAT;
        let mut result = start;
        for i in start..end {
            let v = f(i, context);
            if v < best {
                best = v;
                result = i;
            }
        }
        unsafe { *smallest = best };
        result
    } else {
        /* Try to find minimum faster by recursively checking multiple points. */
        let num = 9;  /* Good value: 9. ?!?!?!?! */
        let mut p = vec![0; num];
        let mut vp = vec![0.0; num];
        let mut besti;
        let mut best;
        let mut lastbest = ZOPFLI_LARGE_FLOAT;
        let mut pos = start;

        loop {
            if end - start <= num {
                break;
            }

            for i in 0..num {
                p[i] = start + (i + 1) * ((end - start) / (num + 1));
                vp[i] = f(p[i], context);
            }

            besti = 0;
            best = vp[0];

            for i in 1..num {
                if vp[i] < best {
                  best = vp[i];
                  besti = i;
                }
            }
            if best > lastbest {
                break;
            }

            start = if besti == 0 { start } else { p[besti - 1] };
            end = if besti == num - 1 { end } else { p[besti + 1] };

            pos = p[besti];
            lastbest = best;
        }
        unsafe { *smallest = lastbest };
        pos
    }
}

/// Returns estimated cost of a block in bits.  It includes the size to encode the
/// tree and the size to encode all literal, length and distance symbols and their
/// extra bits.
///
/// litlens: lz77 lit/lengths
/// dists: ll77 distances
/// lstart: start of block
/// lend: end of block (not inclusive)
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn EstimateCost(lz77: *const ZopfliLZ77Store, lstart: size_t, lend: size_t) -> c_double {
    ZopfliCalculateBlockSizeAutoType(lz77, lstart, lend)
}

/// Gets the cost which is the sum of the cost of the left and the right section
/// of the data.
/// type: FindMinimumFun
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn SplitCost(i: size_t, context: *const c_void) -> c_double {
    let c = unsafe {
        assert!(!context.is_null());
        &*(context as *const SplitCostContext)
    };

    EstimateCost(c.lz77, c.start, i) + EstimateCost(c.lz77, i, c.end)
}

#[repr(C)]
pub struct SplitCostContext {
    lz77: *const ZopfliLZ77Store,
    start: size_t,
    end: size_t,
}

/// Finds next block to try to split, the largest of the available ones.
/// The largest is chosen to make sure that if only a limited amount of blocks is
/// requested, their sizes are spread evenly.
/// lz77size: the size of the LL77 data, which is the size of the done array here.
/// done: array indicating which blocks starting at that position are no longer
///     splittable (splitting them increases rather than decreases cost).
/// splitpoints: the splitpoints found so far.
/// npoints: the amount of splitpoints found so far.
/// lstart: output variable, giving start of block.
/// lend: output variable, giving end of block.
/// returns 1 if a block was found, 0 if no block found (all are done).
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn FindLargestSplittableBlock(lz77size: size_t, done: *const c_uchar, splitpoints: *const size_t, npoints: size_t, lstart: *mut size_t, lend: *mut size_t) -> c_int {
    let mut longest = 0;
    let mut found = 0;

    for i in 0..(npoints + 1) {
        let start = if i == 0 { 0 } else { unsafe { *splitpoints.offset(i as isize - 1) } };
        let end = if i == npoints { lz77size - 1 } else { unsafe { *splitpoints.offset(i as isize) } };
        if unsafe { *done.offset(start as isize) } == 0 && end - start > longest {
            unsafe {
                *lstart = start;
                *lend = end;
            }
            found = 1;
            longest = end - start;
        }
    }

    found
}
