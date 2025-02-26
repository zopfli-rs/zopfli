#!/bin/bash
set -eu
set -o pipefail

# Clone Google/zopfli repo
if [ ! -d "./zopfli" ]
then
    git clone https://github.com/google/zopfli
fi

# Build Google/zopfli
CFLAGS='-fno-omit-frame-pointer' CXXFLAGS='-fno-omit-frame-pointer' \
  make -C zopfli

# Build zopfli-rs
cargo build --release