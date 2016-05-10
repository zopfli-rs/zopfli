use libc::{c_uint};

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
