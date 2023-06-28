use std::io::{Error, ErrorKind, Read};

/// Number of distinct literal/length symbols in DEFLATE
pub const ZOPFLI_NUM_LL: usize = 288;
/// Number of distinct distance symbols in DEFLATE
pub const ZOPFLI_NUM_D: usize = 32;

/// The window size for deflate. Must be a power of two. This should be 32768, the
/// maximum possible by the deflate spec. Anything less hurts compression more than
/// speed.
pub const ZOPFLI_WINDOW_SIZE: usize = 32768;

/// The window mask used to wrap indices into the window. This is why the
/// window size must be a power of two.
pub const ZOPFLI_WINDOW_MASK: usize = ZOPFLI_WINDOW_SIZE - 1;

/// Maximum length that can be encoded in deflate.
pub const ZOPFLI_MAX_MATCH: usize = 258;
/// Minimum length that can be encoded in deflate.
pub const ZOPFLI_MIN_MATCH: usize = 3;

/// For longest match cache. max 256. Uses huge amounts of memory but makes it
/// faster. Uses this many times three bytes per single byte of the input data.
/// This is so because longest match finding has to find the exact distance
/// that belongs to each length for the best lz77 strategy.
/// Good values: e.g. 5, 8.
pub const ZOPFLI_CACHE_LENGTH: usize = 8;

/// limit the max hash chain hits for this hash value. This has an effect only
/// on files where the hash value is the same very often. On these files, this
/// gives worse compression (the value should ideally be 32768, which is the
/// `ZOPFLI_WINDOW_SIZE`, while zlib uses 4096 even for best level), but makes it
/// faster on some specific files.
/// Good value: e.g. 8192.
pub const ZOPFLI_MAX_CHAIN_HITS: usize = 8192;

/// A hasher that may be used by [`HashingAndCountingRead`].
pub trait Hasher {
    fn update(&mut self, data: &[u8]);
}

#[cfg(feature = "gzip")]
impl Hasher for &mut crc32fast::Hasher {
    fn update(&mut self, data: &[u8]) {
        crc32fast::Hasher::update(self, data)
    }
}

#[cfg(feature = "zlib")]
impl Hasher for &mut simd_adler32::Adler32 {
    fn update(&mut self, data: &[u8]) {
        simd_adler32::Adler32::write(self, data)
    }
}

/// A reader that wraps another reader, a hasher and an optional counter,
/// updating the hasher state and incrementing a counter of bytes read so
/// far for each block of data read.
#[cfg(feature = "std")]
pub struct HashingAndCountingRead<'counter, R: Read, H: Hasher> {
    inner: R,
    hasher: H,
    bytes_read: Option<&'counter mut u32>,
}

impl<'counter, R: Read, H: Hasher> HashingAndCountingRead<'counter, R, H> {
    pub fn new(inner: R, hasher: H, bytes_read: Option<&'counter mut u32>) -> Self {
        Self {
            inner,
            hasher,
            bytes_read,
        }
    }
}

impl<R: Read, H: Hasher> Read for HashingAndCountingRead<'_, R, H> {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, Error> {
        match self.inner.read(buf) {
            Ok(bytes_read) => {
                self.hasher.update(&buf[..bytes_read]);

                if let Some(total_bytes_read) = &mut self.bytes_read {
                    **total_bytes_read = total_bytes_read
                        .checked_add(bytes_read.try_into().map_err(|_| ErrorKind::Other)?)
                        .ok_or(ErrorKind::Other)?;
                }

                Ok(bytes_read)
            }
            Err(err) => Err(err),
        }
    }
}
