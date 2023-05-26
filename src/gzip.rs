use std::io::Read;

use crate::{
    deflate::{deflate, BlockType},
    util::HashingAndCountingRead,
    Error, Options, Write,
};

static HEADER: &[u8] = &[
    31,  // ID1
    139, // ID2
    8,   // CM
    0,   // FLG
    0,   // MTIME
    0, 0, 0, 2, // XFL, 2 indicates best compression.
    3, // OS follows Unix conventions.
];

/// Compresses the data according to the gzip specification, RFC 1952.
pub fn gzip_compress<R: Read, W: Write>(
    options: &Options,
    in_data: R,
    mut out: W,
) -> Result<(), Error> {
    let mut crc_hasher = crc32fast::Hasher::new();
    let mut insize = 0;

    let in_data = HashingAndCountingRead::new(in_data, &mut crc_hasher, Some(&mut insize));

    out.write_all(HEADER)?;

    deflate(options, BlockType::Dynamic, in_data, &mut out)?;

    out.write_all(&crc_hasher.finalize().to_le_bytes())?;

    out.write_all(&insize.to_le_bytes())
}
