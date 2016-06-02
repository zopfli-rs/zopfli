use std::io::prelude::*;
use std::fs::File;
use std::ffi::CStr;
use std::{ptr, slice};

use libc::{c_int, c_uchar, size_t, c_char};

use deflate::ZopfliDeflate;
use gzip::ZopfliGzipCompress;
use zlib::ZopfliZlibCompress;

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

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCompress(options_ptr: *const ZopfliOptions, output_type: ZopfliFormat, in_data: &[u8], out: *const *const c_uchar, outsize: *const size_t) {
    match output_type {
        ZopfliFormat::ZOPFLI_FORMAT_GZIP => ZopfliGzipCompress(options_ptr, in_data, out, outsize),
        ZopfliFormat::ZOPFLI_FORMAT_ZLIB => ZopfliZlibCompress(options_ptr, in_data, out, outsize),
        ZopfliFormat::ZOPFLI_FORMAT_DEFLATE => {
            let mut bp = 0;
            let bp_ptr: *mut c_uchar = &mut bp;
            ZopfliDeflate(options_ptr, 2 /* Dynamic block */, 1, in_data, bp_ptr, out, outsize);
        }
    }
}

/// outfilename: filename to write output to, or 0 to write to stdout instead
#[no_mangle]
#[allow(non_snake_case)]
pub extern fn CompressFile(options_ptr: *const ZopfliOptions, output_type: ZopfliFormat, infilename: *const c_char, outfilename: *const c_char) {

    let infilename_str = unsafe {
        assert!(!infilename.is_null());
        CStr::from_ptr(infilename).to_str().expect("Invalid input filename string")
    };

    // TODO: Allow specifying output to STDOUT
    let outfilename_str = unsafe {
        assert!(!outfilename.is_null());
        CStr::from_ptr(outfilename).to_str().expect("Invalid output filename string")
    };

    let mut file = match File::open(infilename_str) {
        Err(why) => panic!("couldn't open {}: {}", infilename_str, why),
        Ok(file) => file,
    };

    let mut in_data = vec![];
    file.read_to_end(&mut in_data).expect("couldn't read the input file");

    let out = ptr::null();
    let mut outsize = 0;
    let outsize_ptr: *mut size_t = &mut outsize;

    ZopfliCompress(options_ptr, output_type, &in_data, &out, outsize_ptr);

    let out_slice = unsafe { slice::from_raw_parts(out, outsize) };

    let mut buffer = File::create(outfilename_str).expect("couldn't create output file");

    buffer.write_all(out_slice).expect("couldn't write to output file");
}
