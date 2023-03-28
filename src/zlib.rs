use std::io::{self, Read, Write};

use crate::deflate::{deflate, BlockType};
use crate::util::HashingAndCountingRead;
use crate::Options;

pub fn zlib_compress<R, W>(options: &Options, in_data: R, mut out: W) -> io::Result<()>
where
    R: Read,
    W: Write,
{
    let cmf = 120; /* CM 8, CINFO 7. See zlib spec.*/
    let flevel = 3;
    let fdict = 0;
    let mut cmfflg: u16 = 256 * cmf + fdict * 32 + flevel * 64;
    let fcheck = 31 - cmfflg % 31;
    cmfflg += fcheck;

    let mut rolling_adler = simd_adler32::Adler32::new();

    let in_data = HashingAndCountingRead::new(in_data, &mut rolling_adler, None);

    out.write_all(&cmfflg.to_be_bytes())?;

    deflate(options, BlockType::Dynamic, in_data, out.by_ref())?;

    out.write_all(&rolling_adler.finish().to_be_bytes())
}
