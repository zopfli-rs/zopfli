use crc::crc32;

use deflate::{deflate, BlockType};
use Options;

/// Compresses the data according to the gzip specification, RFC 1952.
pub fn gzip_compress(options: &Options, in_data: &[u8], out: &mut Vec<u8>) {
    out.push(31);  /* ID1 */
    out.push(139);  /* ID2 */
    out.push(8);  /* CM */
    out.push(0);  /* FLG */
    /* MTIME */
    out.push(0);
    out.push(0);
    out.push(0);
    out.push(0);

    out.push(2);  /* XFL, 2 indicates best compression. */
    out.push(3);  /* OS follows Unix conventions. */

    deflate(options, BlockType::Dynamic, true, in_data, out);

    let crc = crc32::checksum_ieee(in_data);

    out.push(crc as u8);
    out.push((crc >> 8) as u8);
    out.push((crc >> 16) as u8);
    out.push((crc >> 24) as u8);

    let input_size = in_data.len() as u32;
    out.push(input_size as u8);
    out.push((input_size >> 8) as u8);
    out.push((input_size >> 16) as u8);
    out.push((input_size >> 24) as u8);
}
