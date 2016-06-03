use std::io;

use adler32::adler32;
use libc::{c_uchar, c_uint, c_double};

use deflate::deflate;
use Options;

pub fn zlib_compress(options: &Options, in_data: &[u8], out: &mut Vec<u8>) {
    let mut bp = 0;
    let bp_ptr: *mut c_uchar = &mut bp;

    let checksum = adler32(io::Cursor::new(&in_data)).expect("Error with adler32");
    let cmf = 120;  /* CM 8, CINFO 7. See zlib spec.*/
    let flevel = 3;
    let fdict = 0;
    let mut cmfflg: c_uint = 256 * cmf + fdict * 32 + flevel * 64;
    let fcheck = 31 - cmfflg % 31;
    cmfflg += fcheck;

    out.push((cmfflg / 256) as c_uchar);
    out.push((cmfflg % 256) as c_uchar);

    deflate(options, 2 /* dynamic block */, 1 /* final */, in_data, bp_ptr, out);

    out.push(((checksum >> 24) % 256) as c_uchar);
    out.push(((checksum >> 16) % 256) as c_uchar);
    out.push(((checksum >> 8) % 256) as c_uchar);
    out.push((checksum % 256) as c_uchar);

    if options.verbose != 0 {
        let insize = in_data.len();
        let outsize = out.len();
        println!("Original Size: {}, Zlib: {}, Compression: {}% Removed", insize, outsize, 100.0 * (insize - outsize) as c_double / insize as c_double);
    }
}
