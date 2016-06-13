# Zopfli in Rust

[![Build Status](https://travis-ci.org/carols10cents/zopfli.svg?branch=master)](https://travis-ci.org/carols10cents/zopfli)

This is a reimplementation of the [Zopfli](https://github.com/google/zopfli) compression tool in Rust.

I have totally ignored zopflipng.

More info about why and how I did this can be found in [the slides for a talk I gave about it](https://github.com/carols10cents/rust-out-your-c-talk).

## How to build

To build the code, run:

```
$ cargo build --release
```

and the executable will be in `target/release/zopfli`.

This should work on stable or beta Rust.

## Running the tests

There are some unit tests, mostly around the boundary package merge algorithm implementation in katajainen.rs, that can be run with:

```
$ cargo test
```

Golden master tests, to check that compressed files are exactly the same as the C implementation would generate, can be run using:

```
$ ./test/run.sh
```

and then checking that git reports no changes to the files in `test/results`.
