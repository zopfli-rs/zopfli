#![feature(iter_arith)]

#![deny(trivial_casts, trivial_numeric_casts)]

extern crate libc;
extern crate crc;
extern crate adler32;

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

use std::io::prelude::*;
use std::fs::File;

use deflate::{deflate, BlockType};
use gzip::gzip_compress;
use zlib::zlib_compress;

/// Options used throughout the program.
pub struct Options {
  /* Whether to print output */
  verbose: bool,
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

fn compress(options: &Options, output_type: &Format, in_data: &[u8], out: &mut Vec<u8>) {
    match *output_type {
        Format::Gzip => gzip_compress(options, in_data, out),
        Format::Zlib => zlib_compress(options, in_data, out),
        Format::Deflate => {
            let mut bp = 0;
            deflate(options, BlockType::Dynamic, true, in_data, &mut bp, out);
        }
    }
}

/// outfilename: filename to write output to, or 0 to write to stdout instead
pub fn compress_file(options: &Options, output_type: &Format, infilename: &str, outfilename: &str) {

    let mut file = match File::open(infilename) {
        Err(why) => panic!("couldn't open {}: {}", infilename, why),
        Ok(file) => file,
    };

    let mut in_data = vec![];
    file.read_to_end(&mut in_data).expect("couldn't read the input file");

    let mut out = vec![];

    compress(options, output_type, &in_data, &mut out);

    let mut buffer = File::create(outfilename).expect("couldn't create output file");

    buffer.write_all(&out).expect("couldn't write to output file");
}
