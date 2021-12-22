use adler32::adler32;
use byteorder::{BigEndian, WriteBytesExt};
use std::io::{self, Write};

use crate::deflate::{deflate, BlockType};
use crate::Options;

pub fn zlib_compress<W>(options: &Options, in_data: &[u8], mut out: W) -> io::Result<()>
where
    W: Write,
{
    let cmf = 120; /* CM 8, CINFO 7. See zlib spec.*/
    let flevel = 3;
    let fdict = 0;
    let mut cmfflg = 256 * cmf + fdict * 32 + flevel * 64;
    let fcheck = 31 - cmfflg % 31;
    cmfflg += fcheck;

    out.by_ref().write_u16::<BigEndian>(cmfflg)?;

    deflate(options, BlockType::Dynamic, in_data, out.by_ref())?;

    let checksum = adler32(io::Cursor::new(&in_data)).expect("Error with adler32");
    out.write_u32::<BigEndian>(checksum)
}
