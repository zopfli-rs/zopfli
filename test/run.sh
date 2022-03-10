#!/bin/bash
set -eu
set -o pipefail

# Clean output from previous runs
rm -f test/data/*.gz

# Test all input cases individually so we can report which failed
for input in test/data/*; do
	printf "Compressing ${input}... "
	env RUST_BACKTRACE=1 ./zopfli ${input} && echo "done"
done

# Move newly compressed data to its own directory
mkdir -p test/temp_compressed/
mv test/data/*.gz test/temp_compressed/

# Compare newly compressed data with expected
diff -r test/results/ test/temp_compressed/

# Verify that all compressed files decompress back to the input sources
mkdir -p test/temp_decompressed/
for gz in test/results/*.gz; do
	BASENAME=$(basename -s .gz $gz)
	printf "Validating ${gz}... "
	cat $gz | gzip -d > test/temp_decompressed/${BASENAME}
	diff -r test/data/${BASENAME} test/temp_decompressed/${BASENAME} && echo "done"
done

# Clean up temporary output
rm -rf test/temp_decompressed/
rm -rf test/temp_compressed/
