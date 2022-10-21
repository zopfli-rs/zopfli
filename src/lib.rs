#![deny(trivial_casts, trivial_numeric_casts)]

//! A reimplementation of the [Zopfli](https://github.com/google/zopfli) compression library in Rust.
//!
//! Zopfli is a state of the art DEFLATE compressor that heavily prioritizes compression over speed.
//! It usually compresses much better than other DEFLATE compressors, generating standard DEFLATE
//! streams that can be decompressed with any DEFLATE decompressor, at the cost of being
//! significantly slower.

mod blocksplitter;
mod cache;
mod deflate;
mod gzip;
mod hash;
mod iter;
mod katajainen;
mod lz77;
mod squeeze;
mod symbols;
mod tree;
mod util;
mod zlib;

use std::io::{self, Read, Write};
use std::num::NonZeroU8;

use crate::deflate::{deflate, BlockType};
use crate::gzip::gzip_compress;
use crate::zlib::zlib_compress;

/// Options for the Zopfli compression algorithm.
pub struct Options {
    /// Maximum amount of times to rerun forward and backward pass to optimize LZ77
    /// compression cost.
    /// Good values: 10, 15 for small files, 5 for files over several MB in size or
    /// it will be too slow.
    ///
    /// Default value: 15.
    pub iteration_count: NonZeroU8,
    /// Maximum amount of blocks to split into (0 for unlimited, but this can give
    /// extreme results that hurt compression on some files).
    ///
    /// Default value: 15.
    pub maximum_block_splits: u16,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            iteration_count: NonZeroU8::new(15).unwrap(),
            maximum_block_splits: 15,
        }
    }
}

/// The output file format to use to store data compressed with Zopfli.
pub enum Format {
    /// The gzip file format, as defined in
    /// [RFC 1952](https://datatracker.ietf.org/doc/html/rfc1952).
    ///
    /// This file format can be easily decompressed with the gzip
    /// program.
    Gzip,
    /// The zlib file format, as defined in
    /// [RFC 1950](https://datatracker.ietf.org/doc/html/rfc1950).
    ///
    /// The zlib format has less header overhead than gzip, but it
    /// stores less metadata about the compressed data and may not
    /// be as fit for purpose.
    Zlib,
    /// The raw DEFLATE stream format, as defined in
    /// [RFC 1951](https://datatracker.ietf.org/doc/html/rfc1951).
    ///
    /// Raw DEFLATE streams are not meant to be stored raw because
    /// they lack error detection and correction metadata. They
    /// usually are embedded in other file formats, such as gzip
    /// and zlib.
    Deflate,
}

/// Compresses data from a source with the Zopfli algorithm, using the specified
/// options, and writes the result to a sink in the defined output format.
pub fn compress<R, W>(
    options: &Options,
    output_format: &Format,
    in_data: R,
    out: W,
) -> io::Result<()>
where
    R: Read,
    W: Write,
{
    match *output_format {
        Format::Gzip => gzip_compress(options, in_data, out),
        Format::Zlib => zlib_compress(options, in_data, out),
        Format::Deflate => deflate(options, BlockType::Dynamic, in_data, out),
    }
}
