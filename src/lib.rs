#![deny(trivial_casts, trivial_numeric_casts)]

extern crate libc;
extern crate crc;
extern crate adler32;

mod iter;
mod blocksplitter;
mod cache;
mod deflate;
mod gzip;
mod hash;
mod katajainen;
mod lz77;
mod squeeze;
mod symbols;
mod tree;
mod util;
mod zlib;

use deflate::{deflate, BlockType};
use gzip::gzip_compress;
use zlib::zlib_compress;

/// Options used throughout the program.
pub struct Options {
  /* Whether to print output */
  pub verbose: bool,
  /* Whether to print more detailed output */
  verbose_more: bool,
  /*
  Maximum amount of times to rerun forward and backward pass to optimize LZ77
  compression cost. Good values: 10, 15 for small files, 5 for files over
  several MB in size or it will be too slow.
  */
  numiterations: i32,
  /*
  Maximum amount of blocks to split into (0 for unlimited, but this can give
  extreme results that hurt compression on some files). Default value: 15.
  */
  blocksplittingmax: i32,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            verbose: false,
            verbose_more: false,
            numiterations: 15,
            blocksplittingmax: 15,
        }
    }
}

pub enum Format {
    Gzip,
    Zlib,
    Deflate,
}

pub fn compress(options: &Options, output_type: &Format, in_data: &[u8], out: &mut Vec<u8>) {
    match *output_type {
        Format::Gzip => gzip_compress(options, in_data, out),
        Format::Zlib => zlib_compress(options, in_data, out),
        Format::Deflate => deflate(options, BlockType::Dynamic, true, in_data, out),
    }
}
