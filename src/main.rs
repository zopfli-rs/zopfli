use std::env;

extern crate zopfli;

fn main() {
    let options = zopfli::Options::new();
    let output_type = zopfli::ZopfliFormat::ZOPFLI_FORMAT_GZIP;

    // TODO: CLI arguments
    // TODO: Allow specifying output to STDOUT

    let extension = match output_type {
        zopfli::ZopfliFormat::ZOPFLI_FORMAT_GZIP => ".gz",
        zopfli::ZopfliFormat::ZOPFLI_FORMAT_ZLIB => ".zlib",
        zopfli::ZopfliFormat::ZOPFLI_FORMAT_DEFLATE => ".deflate",
    };

    for filename in env::args().skip(1) {
        zopfli::CompressFile(&options, &output_type, &filename, &format!("{}{}", filename, extension));
    }
}
