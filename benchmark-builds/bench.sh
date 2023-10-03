#!/bin/bash
set -eu
set -o pipefail

#Get Google's version and build rust version
./prepare.sh

# Clean output from previous runs
rm -f  benchmark-builds*.json
rm -f  data/*.gz
rm -rf data_google data_rust

# Duplicate the input data for the bench execution,
# to avoid any side effect due to either having to
# 'delete the output' or 'overwrite the output' during the bench.
cp -r data data_google
cp -r data data_rust

# Bench all input cases individually,
# and store respective results.
for input in $(ls data | grep -v '\.gz$'); do
	printf "Benching with ${input}...\n"
	command_name="Google zopfli with '${input}'"
	command="zopfli/zopfli data_google/${input}" 
	command_name1="Rust zopfli with '${input}'"
	command1="../target/release/zopfli data_rust/${input}"
	hyperfine \
		--ignore-failure \
		-N \
		--export-json="benchmark-builds_${input}.json" \
		--warmup=3 \
		--time-unit=millisecond \
		--command-name="${command_name}" "${command}" \
		--command-name="${command_name1}" "${command1}"
done