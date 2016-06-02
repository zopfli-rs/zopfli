use std::slice;
use std::io;

use adler32::adler32;
use libc::{c_uchar, size_t, c_uint, c_double};

use deflate::ZopfliDeflate;
use zopfli::ZopfliOptions;

#[link(name = "zopfli")]
extern {
    fn ZopfliAppendDataUChar(value: c_uchar, data: *const *const c_uchar, size: *const size_t);
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliZlibCompress(options_ptr: *const ZopfliOptions, in_data: *const c_uchar, insize: size_t, out: *const *const c_uchar, outsize: *const size_t) {
    let options = unsafe {
        assert!(!options_ptr.is_null());
        &*options_ptr
    };

    let mut bp = 0;
    let bp_ptr: *mut c_uchar = &mut bp;

    let input_bytes = unsafe { slice::from_raw_parts(in_data, insize) };

    let checksum = adler32(io::Cursor::new(&input_bytes)).expect("Error with adler32");
    let cmf = 120;  /* CM 8, CINFO 7. See zlib spec.*/
    let flevel = 3;
    let fdict = 0;
    let mut cmfflg: c_uint = 256 * cmf + fdict * 32 + flevel * 64;
    let fcheck = 31 - cmfflg % 31;
    cmfflg += fcheck;

    unsafe {
        ZopfliAppendDataUChar((cmfflg / 256) as c_uchar, out, outsize);
        ZopfliAppendDataUChar((cmfflg % 256) as c_uchar, out, outsize);
    }

    ZopfliDeflate(options, 2 /* dynamic block */, 1 /* final */, in_data, insize, bp_ptr, out, outsize);

    unsafe {
        ZopfliAppendDataUChar(((checksum >> 24) % 256) as c_uchar, out, outsize);
        ZopfliAppendDataUChar(((checksum >> 16) % 256) as c_uchar, out, outsize);
        ZopfliAppendDataUChar(((checksum >> 8) % 256) as c_uchar, out, outsize);
        ZopfliAppendDataUChar((checksum % 256) as c_uchar, out, outsize);
    }

    if options.verbose != 0 {
        unsafe {
            println!("Original Size: {}, Zlib: {}, Compression: {}% Removed", insize, *outsize, 100.0 * (insize - *outsize) as c_double / insize as c_double);
        }
    }
}
