#!/bin/bash
set -eu
set -o pipefail

# Clone Google/zopfli repo
if [ ! -d "./zopfli" ]
then
    git clone https://github.com/google/zopfli
fi

# Build Google/zopfli
cd zopfli && make && cd -

# Build zopfli-rs
cd .. && cargo build --release && cd -