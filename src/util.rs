use libc::size_t;

/// Number of distinct literal/length symbols in DEFLATE
pub const ZOPFLI_NUM_LL: size_t = 288;
/// Number of distinct distance symbols in DEFLATE
pub const ZOPFLI_NUM_D: size_t = 32;

/// The window size for deflate. Must be a power of two. This should be 32768, the
/// maximum possible by the deflate spec. Anything less hurts compression more than
/// speed.
pub const ZOPFLI_WINDOW_SIZE: size_t = 32768;

/// The window mask used to wrap indices into the window. This is why the
/// window size must be a power of two.
pub const ZOPFLI_WINDOW_MASK: size_t = 32767; // ZOPFLI_WINDOW_SIZE - 1

/// Maximum length that can be encoded in deflate.
pub const ZOPFLI_MAX_MATCH: size_t = 258;
/// Minimum length that can be encoded in deflate.
pub const ZOPFLI_MIN_MATCH: size_t = 3;

/// For longest match cache. max 256. Uses huge amounts of memory but makes it
/// faster. Uses this many times three bytes per single byte of the input data.
/// This is so because longest match finding has to find the exact distance
/// that belongs to each length for the best lz77 strategy.
/// Good values: e.g. 5, 8.
pub const ZOPFLI_CACHE_LENGTH: size_t = 8;

/// limit the max hash chain hits for this hash value. This has an effect only
/// on files where the hash value is the same very often. On these files, this
/// gives worse compression (the value should ideally be 32768, which is the
/// ZOPFLI_WINDOW_SIZE, while zlib uses 4096 even for best level), but makes it
/// faster on some specific files.
/// Good value: e.g. 8192.
pub const ZOPFLI_MAX_CHAIN_HITS: size_t = 8192;

/// A block structure of huge, non-smart, blocks to divide the input into, to allow
/// operating on huge files without exceeding memory, such as the 1GB wiki9 corpus.
/// The whole compression algorithm, including the smarter block splitting, will
/// be executed independently on each huge block.
/// Dividing into huge blocks hurts compression, but not much relative to the size.
pub const ZOPFLI_MASTER_BLOCK_SIZE: size_t = 1000000;
