use libc::c_int;

/// Options used throughout the program.
#[repr(C)]
pub struct ZopfliOptions {
  /* Whether to print output */
  verbose: c_int,
  /* Whether to print more detailed output */
  verbose_more: c_int,
  /*
  Maximum amount of times to rerun forward and backward pass to optimize LZ77
  compression cost. Good values: 10, 15 for small files, 5 for files over
  several MB in size or it will be too slow.
  */
  numiterations: c_int,
  /*
  No longer used, left for compatibility.
  */
  blocksplittinglast: c_int,
  /*
  Maximum amount of blocks to split into (0 for unlimited, but this can give
  extreme results that hurt compression on some files). Default value: 15.
  */
  blocksplittingmax: c_int,
}
