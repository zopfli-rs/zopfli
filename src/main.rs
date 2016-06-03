use std::env;

extern crate zopfli;

fn main() {
    let options = zopfli::Options::new();
    let output_type = zopfli::Format::Gzip;

    // TODO: CLI arguments
    // TODO: Allow specifying output to STDOUT

    let extension = match output_type {
        zopfli::Format::Gzip => ".gz",
        zopfli::Format::Zlib => ".zlib",
        zopfli::Format::Deflate => ".deflate",
    };

    for filename in env::args().skip(1) {
        zopfli::compress_file(&options, &output_type, &filename, &format!("{}{}", filename, extension));
    }
}
