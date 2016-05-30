use libc::{size_t, c_void, c_double};

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
