use std::io::{self, Write};
use crc::crc32;
use byteorder::{LittleEndian, WriteBytesExt};

use deflate::{deflate, BlockType};
use Options;

static HEADER: &'static [u8] = &[
    31,  // ID1
    139, // ID2
    8,   // CM
    0,   // FLG

    0,   // MTIME
    0,
    0,
    0,

    2,   // XFL, 2 indicates best compression.
    3,   // OS follows Unix conventions.
];

/// Compresses the data according to the gzip specification, RFC 1952.
pub fn gzip_compress<W>(options: &Options, in_data: &[u8], mut out: W) -> io::Result<()>
    where W: Write
{
    try!(out.by_ref().write_all(HEADER));

    try!(deflate(options, BlockType::Dynamic, in_data, out.by_ref()));

    try!(out.by_ref().write_u32::<LittleEndian>(crc32::checksum_ieee(in_data)));
    out.write_u32::<LittleEndian>(in_data.len() as u32)
}
