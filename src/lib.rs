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
        Format::Deflate => deflate(options, BlockType::Dynamic, true, in_data, out),
    }
}

/// outfilename: filename to write output to, or 0 to write to stdout instead
pub fn compress_file(options: &Options, output_type: &Format, infilename: &str, outfilename: &str) {
    let mut file = File::open(infilename)
        .unwrap_or_else(|why| panic!("couldn't open {}: {}", infilename, why));

    let mut in_data = Vec::with_capacity(file.metadata().map(|x| x.len()).unwrap_or(0) as usize);
    let mut out_data = vec![];

    // Read the contents of the file into in_data
    file.read_to_end(&mut in_data).ok()
        // Panic if the file could not be read
        .map_or_else(|| panic!("couldn't read the input file"), |_| {
            // Compress `in_data` and store the result in `out_data`
            compress(options, output_type, &in_data, &mut out_data);
            // Attempt to create the output file
            File::create(outfilename).ok()
                // Panic if the output file could not be opened
                .map_or_else(|| panic!("couldn't create output file"), |mut buffer| {
                    // Write the `out_data` into the newly created file.
                    buffer.write_all(&out_data).expect("couldn't write to output file")
                });
        });
}
