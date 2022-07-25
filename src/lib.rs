#![deny(trivial_casts, trivial_numeric_casts)]

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

use std::io::{self, Read, Seek, SeekFrom, Write};

use crate::deflate::{deflate, BlockType};
use crate::gzip::gzip_compress;
use crate::zlib::zlib_compress;

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

pub fn compress_seekable<R, W>(
    options: &Options,
    output_type: &Format,
    mut in_data: R,
    out: W,
) -> io::Result<()>
where
    R: Read + Seek,
    W: Write,
{
    let insize = {
        let offset = in_data.seek(SeekFrom::Current(0))?;
        let insize = in_data.seek(SeekFrom::End(0))?;
        in_data.seek(SeekFrom::Start(offset))?;
        insize
    };

    compress(options, output_type, in_data, insize, out)
}

pub fn compress<R, W>(
    options: &Options,
    output_type: &Format,
    in_data: R,
    insize: u64,
    out: W,
) -> io::Result<()>
where
    R: Read,
    W: Write,
{
    match *output_type {
        Format::Gzip => gzip_compress(options, in_data, insize, out),
        Format::Zlib => zlib_compress(options, in_data, insize, out),
        Format::Deflate => deflate(options, BlockType::Dynamic, in_data, insize, out),
    }
}
