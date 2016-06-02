use crc::crc32;
use libc::{c_uchar, c_double, size_t};

use deflate::ZopfliDeflate;
use zopfli::ZopfliOptions;

#[link(name = "zopfli")]
extern {
    fn ZopfliAppendDataUChar(value: c_uchar, data: *const *const c_uchar, size: *const size_t);
}

/// Compresses the data according to the gzip specification, RFC 1952.
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliGzipCompress(options_ptr: *const ZopfliOptions, in_data: &[u8], out: *const *const c_uchar, outsize: *const size_t) {
    let options = unsafe {
        assert!(!options_ptr.is_null());
        &*options_ptr
    };

    unsafe {
        ZopfliAppendDataUChar(31, out, outsize);  /* ID1 */
        ZopfliAppendDataUChar(139, out, outsize);  /* ID2 */
        ZopfliAppendDataUChar(8, out, outsize);  /* CM */
        ZopfliAppendDataUChar(0, out, outsize);  /* FLG */
        /* MTIME */
        ZopfliAppendDataUChar(0, out, outsize);
        ZopfliAppendDataUChar(0, out, outsize);
        ZopfliAppendDataUChar(0, out, outsize);
        ZopfliAppendDataUChar(0, out, outsize);

        ZopfliAppendDataUChar(2, out, outsize);  /* XFL, 2 indicates best compression. */
        ZopfliAppendDataUChar(3, out, outsize);  /* OS follows Unix conventions. */
    }

    let mut bp = 0;
    let bp_ptr: *mut c_uchar = &mut bp;

    ZopfliDeflate(options_ptr, 2 /* Dynamic block */, 1, in_data, bp_ptr, out, outsize);

    let crc = crc32::checksum_ieee(in_data);

    unsafe {
        ZopfliAppendDataUChar((crc >> 0) as u8, out, outsize);
        ZopfliAppendDataUChar((crc >> 8) as u8, out, outsize);
        ZopfliAppendDataUChar((crc >> 16) as u8, out, outsize);
        ZopfliAppendDataUChar((crc >> 24) as u8, out, outsize);
    }

    let input_size = in_data.len() as u32;
    unsafe {
        ZopfliAppendDataUChar((input_size >> 0) as u8, out, outsize);
        ZopfliAppendDataUChar((input_size >> 8) as u8, out, outsize);
        ZopfliAppendDataUChar((input_size >> 16) as u8, out, outsize);
        ZopfliAppendDataUChar((input_size >> 24) as u8, out, outsize);
    }

    if options.verbose != 0 {
        let insize = in_data.len();
        unsafe {
            println!("Original Size: {}, Gzip: {}, Compression: {}% Removed", insize, *outsize, 100.0 * (insize - *outsize) as c_double / insize as c_double);
        }
    }
}
