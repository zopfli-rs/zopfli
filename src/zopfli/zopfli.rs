use libc::{c_int, c_uchar, size_t};

use deflate::ZopfliDeflate;
use gzip::ZopfliGzipCompress;

/// Options used throughout the program.
#[repr(C)]
pub struct ZopfliOptions {
  /* Whether to print output */
  pub verbose: c_int,
  /* Whether to print more detailed output */
  pub verbose_more: c_int,
  /*
  Maximum amount of times to rerun forward and backward pass to optimize LZ77
  compression cost. Good values: 10, 15 for small files, 5 for files over
  several MB in size or it will be too slow.
  */
  pub numiterations: c_int,
  /*
  Maximum amount of blocks to split into (0 for unlimited, but this can give
  extreme results that hurt compression on some files). Default value: 15.
  */
  pub blocksplittingmax: c_int,
}

#[repr(C)]
pub enum ZopfliFormat {
  ZOPFLI_FORMAT_GZIP,
  ZOPFLI_FORMAT_ZLIB,
  ZOPFLI_FORMAT_DEFLATE
}

#[link(name = "zopfli")]
extern {
    fn ZopfliZlibCompress(options_ptr: *const ZopfliOptions, in_data: *const c_uchar, insize: size_t, out: *const *const c_uchar, outsize: *const size_t);
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCompress(options_ptr: *const ZopfliOptions, output_type: ZopfliFormat, in_data: *const c_uchar, insize: size_t, out: *const *const c_uchar, outsize: *const size_t) {
    match output_type {
        ZopfliFormat::ZOPFLI_FORMAT_GZIP => ZopfliGzipCompress(options_ptr, in_data, insize, out, outsize),
        ZopfliFormat::ZOPFLI_FORMAT_ZLIB => unsafe { ZopfliZlibCompress(options_ptr, in_data, insize, out, outsize) },
        ZopfliFormat::ZOPFLI_FORMAT_DEFLATE => {
            let mut bp = 0;
            let bp_ptr: *mut c_uchar = &mut bp;
            ZopfliDeflate(options_ptr, 2 /* Dynamic block */, 1, in_data, insize, bp_ptr, out, outsize);
        }
    }
}
