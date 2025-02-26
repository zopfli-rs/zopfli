#!/bin/sh -eu

# Helper script to record profiles of Zopfli's CPU usage on a Linux computer with `perf`,
# and display them with `samply`. The profiles are not recorded directly with `samply`
# because doing so does not simplify the script much anyway, `perf`'s output format is
# more interoperable with other visualization frontends, and `perf record` exposes
# much more native features and knobs that are useful.

mkdir -p test/perf/profiles

if [ -d /sys/devices/cpu_core ]; then
  # For hybrid Intel CPUs (Alder Lake+, with performance and efficiency cores),
  # only record cycles for the performance cores Zopfli will be scheduled to,
  # as Samply does not have proper support for combining both efficiency and
  # performance events
  readonly CPU_CYCLES_EVENT='cpu_core/cycles/'
else
  readonly CPU_CYCLES_EVENT='cpu-cycles'
fi

case "$1" in
  'rust')
    readonly BINARY_PATH='target/release/zopfli'
    readonly PROFILE_FILE_PREFIX='rust'
    cargo build --release;;
  'c')
    # Run the prepare.sh script first to get the Zopfli sources
    readonly BINARY_PATH='zopfli/zopfli'
    readonly PROFILE_FILE_PREFIX='c'
    make -C zopfli;;
  *)
    echo 'Invalid Zopfli flavor specified, expected either "rust" or "c"' >&2
    exit 1;;
esac
shift

PERF_DATA_FILE="test/perf/profiles/${PROFILE_FILE_PREFIX}_zopfli_$(date +%s).perf.data"
readonly PERF_DATA_FILE

(set -x; perf record \
  --call-graph=fp -F 5000 \
  --event="$CPU_CYCLES_EVENT" \
  -o "$PERF_DATA_FILE" -- "$BINARY_PATH" "$@")

PROCESSED_PROFILE="${TMPDIR:-/tmp}/profile-$(date +%s).bin"
readonly PROCESSED_PROFILE
trap 'rm -f "$PROCESSED_PROFILE" || true' EXIT INT TERM

(set -x; samply import -o "$PROCESSED_PROFILE" "$PERF_DATA_FILE")
