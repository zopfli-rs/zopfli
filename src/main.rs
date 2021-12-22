use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};

extern crate zopfli;

fn main() {
    let options = zopfli::Options::default();
    let output_type = zopfli::Format::Gzip;

    // TODO: CLI arguments
    // TODO: Allow specifying output to STDOUT

    let extension = match output_type {
        zopfli::Format::Gzip => ".gz",
        zopfli::Format::Zlib => ".zlib",
        zopfli::Format::Deflate => ".deflate",
    };

    for filename in env::args().skip(1) {
        let mut file = File::open(&filename)
            .unwrap_or_else(|why| panic!("couldn't open {}: {}", filename, why));
        let filesize = file.metadata().map(|x| x.len()).unwrap_or(0) as usize;

        let mut data = Vec::with_capacity(filesize);
        // Read the contents of the file into in_data; panic if the file could not be read
        file.read_to_end(&mut data)
            .unwrap_or_else(|why| panic!("couldn't read {}: {}", filename, why));

        let out_filename = format!("{}{}", filename, extension);

        // Attempt to create the output file, panic if the output file could not be opened
        let out_file = File::create(&out_filename)
            .unwrap_or_else(|why| panic!("couldn't create output file {}: {}", out_filename, why));
        let mut out_file = WriteStatistics::new(BufWriter::new(out_file));

        zopfli::compress(&options, &output_type, &data, &mut out_file).unwrap_or_else(|why| {
            panic!("couldn't write to output file {}: {}", out_filename, why)
        });

        if options.verbose {
            let out_size = out_file.count;
            println!(
                "Original Size: {}, Compressed: {}, Compression: {}% Removed",
                filesize,
                out_size,
                100.0 * (filesize - out_size) as f64 / filesize as f64
            );
        }
    }
}

struct WriteStatistics<W> {
    inner: W,
    count: usize,
}

impl<W> WriteStatistics<W> {
    fn new(inner: W) -> Self {
        WriteStatistics {
            inner: inner,
            count: 0,
        }
    }
}

impl<W> Write for WriteStatistics<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let res = self.inner.write(buf);
        if let Ok(size) = res {
            self.count += size;
        }
        res
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}
