CC = gcc
CXX = g++

CFLAGS = -W -Wall -Wextra -ansi -pedantic -lm -O2 -Wno-unused-function
CXXFLAGS = -W -Wall -Wextra -ansi -pedantic -O2

ZOPFLI_RUST_DEBUG := target/debug/libzopfli.a
ZOPFLI_RUST_RELEASE := target/release/libzopfli.a
ZOPFLILIB_SRC = src/zopfli/deflate.c src/zopfli/hash.c\
                src/zopfli/util.c src/zopfli/zopfli_lib.c
ZOPFLILIB_OBJ := $(patsubst src/zopfli/%.c,%.o,$(ZOPFLILIB_SRC))
ZOPFLIBIN_SRC := src/zopfli/zopfli_bin.c

.PHONY: zopfli

.PHONY: target/debug/libzopfli.a
target/debug/libzopfli.a:
	cargo build --verbose

.PHONY: target/release/libzopfli.a
target/release/libzopfli.a:
	cargo build --verbose --release

# Zopfli binary
zopfli: $(ZOPFLI_RUST_RELEASE)
	$(CC) $(ZOPFLI_RUST_RELEASE) $(ZOPFLILIB_SRC) $(ZOPFLIBIN_SRC) $(CFLAGS) -o zopfli

# Zopfli debug binary
zopflidebug: $(ZOPFLI_RUST_DEBUG)
	$(CC) $(ZOPFLI_RUST_DEBUG) $(ZOPFLILIB_SRC) $(ZOPFLIBIN_SRC) $(CFLAGS) -o zopfli

# Zopfli shared library
libzopfli:
	$(CC) $(ZOPFLILIB_SRC) $(CFLAGS) -fPIC -c
	$(CC) $(ZOPFLILIB_OBJ) $(CFLAGS) -shared -Wl,-soname,libzopfli.so.1 -o libzopfli.so.1.0.1

# Remove all libraries and binaries
clean:
	cargo clean && rm -f zopfli $(ZOPFLILIB_OBJ) libzopfli*
