#![deny(trivial_casts, trivial_numeric_casts, missing_docs)]

//! A reimplementation of the [Zopfli](https://github.com/google/zopfli) compression library in Rust.
//!
//! Zopfli is a state of the art DEFLATE compressor that heavily prioritizes compression over speed.
//! It usually compresses much better than other DEFLATE compressors, generating standard DEFLATE
//! streams that can be decompressed with any DEFLATE decompressor, at the cost of being
//! significantly slower.
//!
//! # Features
//!
//! This crate exposes the following features. You can enable or disable them in your `Cargo.toml`
//! as needed.
//!
//! - `gzip` (enabled by default): enables support for compression in the gzip format. Implies `std`.
//! - `zlib` (enabled by default): enables support for compression in the Zlib format. Implies `std`.
//! - `std` (enabled by default): enables linking against the Rust standard library. When not enabled,
//!                               the crate is built with the `#![no_std]` attribute and can be used
//!                               in any environment where [`alloc`](https://doc.rust-lang.org/alloc/)
//!                               (i.e., a memory allocator) is available. In addition, the crate
//!                               exposes minimalist versions of the `std` I/O traits it needs to
//!                               function, allowing users to implement them. Disabling `std` requires
//!                               enabling `nightly` due to dependencies on unstable language features.
//! - `nightly`: enables performance optimizations that are specific to the nightly Rust toolchain.
//!              Currently, this feature improves rustdoc generation and enables the namesake feature
//!              on `crc32fast` and `simd-adler32`, but this may change in the future.

#![cfg_attr(
    feature = "nightly",
    feature(doc_auto_cfg),
    feature(error_in_core),
    feature(core_intrinsics)
)]

#[macro_use]
extern crate alloc;
extern crate core;

pub use deflate::{BlockType, DeflateEncoder};
#[cfg(test)]
use proptest::prelude::*;

mod blocksplitter;
mod cache;
mod deflate;
#[cfg(feature = "gzip")]
mod gzip;
mod hash;
#[cfg(doc)]
mod io;
mod iter;
mod katajainen;
mod lz77;
mod squeeze;
mod symbols;
mod tree;
mod util;
#[cfg(feature = "zlib")]
mod zlib;

use core::num::NonZeroU64;
#[cfg(not(doc))]
use std::io::{Error, Write};

#[cfg(doc)]
pub use io::{Error, ErrorKind, Write};

/// Options for the Zopfli compression algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(test, derive(proptest_derive::Arbitrary))]
pub struct Options {
    /// Maximum amount of times to rerun forward and backward pass to optimize LZ77
    /// compression cost.
    /// Good values: 10, 15 for small files, 5 for files over several MB in size or
    /// it will be too slow.
    ///
    /// Default value: 15.
    #[cfg_attr(
        test,
        proptest(
            strategy = "(1..=10u64).prop_map(|iteration_count| NonZeroU64::new(iteration_count))"
        )
    )]
    pub iteration_count: Option<NonZeroU64>,
    /// Stop after rerunning forward and backward pass this many times without finding
    /// a smaller representation of the block.
    pub iterations_without_improvement: Option<NonZeroU64>,
    /// Maximum amount of blocks to split into (0 for unlimited, but this can give
    /// extreme results that hurt compression on some files).
    ///
    /// Default value: 15.
    pub maximum_block_splits: u16,
}

impl Default for Options {
    fn default() -> Options {
        Options {
            iteration_count: Some(NonZeroU64::new(15).unwrap()),
            iterations_without_improvement: None,
            maximum_block_splits: 15,
        }
    }
}

/// The output file format to use to store data compressed with Zopfli.
#[derive(Debug, Copy, Clone)]
pub enum Format {
    /// The gzip file format, as defined in
    /// [RFC 1952](https://datatracker.ietf.org/doc/html/rfc1952).
    ///
    /// This file format can be easily decompressed with the gzip
    /// program.
    #[cfg(feature = "gzip")]
    Gzip,
    /// The zlib file format, as defined in
    /// [RFC 1950](https://datatracker.ietf.org/doc/html/rfc1950).
    ///
    /// The zlib format has less header overhead than gzip, but it
    /// stores less metadata.
    #[cfg(feature = "zlib")]
    Zlib,
    /// The raw DEFLATE stream format, as defined in
    /// [RFC 1951](https://datatracker.ietf.org/doc/html/rfc1951).
    ///
    /// Raw DEFLATE streams are not meant to be stored as-is because
    /// they lack error detection and correction metadata. They
    /// usually are embedded in other file formats, such as gzip
    /// and zlib.
    Deflate,
}

/// Compresses data from a source with the Zopfli algorithm, using the specified
/// options, and writes the result to a sink in the defined output format.
pub fn compress<R: std::io::Read, W: Write>(
    options: &Options,
    output_format: &Format,
    in_data: R,
    out: W,
) -> Result<(), Error> {
    match output_format {
        #[cfg(feature = "gzip")]
        Format::Gzip => gzip::gzip_compress(*options, in_data, out),
        #[cfg(feature = "zlib")]
        Format::Zlib => zlib::zlib_compress(*options, in_data, out),
        Format::Deflate => deflate::deflate(*options, BlockType::Dynamic, in_data, out),
    }
}

#[cfg(test)]
mod test {
    use std::io;

    use miniz_oxide::inflate;
    use proptest::proptest;

    use super::*;

    proptest! {
        #[test]
        fn deflating_is_reversible(
            options: Options,
            btype: BlockType,
            data in prop::collection::vec(any::<u8>(), 0..64 * 1024)
        ) {
            let mut compressed_data = Vec::with_capacity(data.len());

            let mut encoder = DeflateEncoder::new(options, btype, &mut compressed_data);
            io::copy(&mut &*data, &mut encoder).unwrap();
            encoder.finish().unwrap();

            let decompressed_data = inflate::decompress_to_vec(&compressed_data).expect("Could not inflate compressed stream");
            prop_assert_eq!(data, decompressed_data, "Decompressed data should match input data");
        }
    }
}
