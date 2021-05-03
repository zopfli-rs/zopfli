use std::io::{self, Read, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use iter_read::IterRead;

use crate::deflate::{deflate, BlockType};
use crate::Options;

const CRC_IEEE: crc::Crc<u32> = crc::Crc::<u32>::new(&crc::CRC_32_ISO_HDLC);

static HEADER: &'static [u8] = &[
    31,  // ID1
    139, // ID2
    8,   // CM
    0,   // FLG
    0,   // MTIME
    0, 0, 0, 2, // XFL, 2 indicates best compression.
    3, // OS follows Unix conventions.
];

/// Compresses the data according to the gzip specification, RFC 1952.
pub fn gzip_compress<R, W>(options: &Options, in_data: R, insize: u64, mut out: W) -> io::Result<()>
where
    R: Read,
    W: Write,
{
    let mut crc_digest = CRC_IEEE.digest();
    let mut read_error_kind = None;

    let in_data = IterRead::new(
        in_data
            .bytes()
            .filter_map(|byte_result| {
                read_error_kind = byte_result.as_ref().map_or_else(
                    |error| Some(error.kind()),
                    |byte| {
                        crc_digest.update(&[*byte]);
                        None
                    },
                );

                byte_result.ok()
            })
            .fuse(),
    );

    out.by_ref().write_all(HEADER)?;

    deflate(options, BlockType::Dynamic, in_data, insize, out.by_ref())?;

    // in_data is fused and stops reading bytes after the first error, so
    // this if is evaluated as soon as an error occurs. The deflate function
    // has received EOF at this point, so the last block has been written.
    if let Some(error_kind) = read_error_kind {
        return Err(error_kind.into());
    }

    out.by_ref()
        .write_u32::<LittleEndian>(crc_digest.finalize())?;

    out.write_u32::<LittleEndian>(insize as u32)
}
