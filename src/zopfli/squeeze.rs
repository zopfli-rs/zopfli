use libc::{c_int, c_void, c_uint, c_double};

use symbols::{ZopfliGetDistExtraBits, ZopfliGetLengthExtraBits, ZopfliGetLengthSymbol};

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
