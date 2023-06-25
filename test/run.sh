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

# Compare newly compressed data size with expected
cd test/results
find . -type f | while read -r file; do
  reference_size=$(stat -c%s "../../test/temp_compressed/$file")
  current_result_size=$(stat -c%s "$file")
  if [[ $current_result_size > $reference_size ]]; then
    echo "File $file is larger than expected ($current_result_size vs $reference_size bytes)"
    exit 1
  elif [[ $current_result_size < $reference_size ]]; then
    echo "Compression ratio improved for $file! ($current_result_size vs $reference_size bytes)"
  fi
done
cd ../..

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
