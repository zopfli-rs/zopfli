use std::io::{self, Read, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::deflate::{deflate, BlockType};
use crate::util::HashingAndCountingRead;
use crate::Options;

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
pub fn gzip_compress<R, W>(options: &Options, in_data: R, mut out: W) -> io::Result<()>
where
    R: Read,
    W: Write,
{
    let mut crc_hasher = crc32fast::Hasher::new();
    let mut insize = 0;

    let in_data = HashingAndCountingRead::new(in_data, &mut crc_hasher, Some(&mut insize));

    out.by_ref().write_all(HEADER)?;

    deflate(options, BlockType::Dynamic, in_data, out.by_ref())?;

    out.by_ref()
        .write_u32::<LittleEndian>(crc_hasher.finalize())?;

    out.write_u32::<LittleEndian>(insize)
}
