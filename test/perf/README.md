## Performance testing and bechmarking

The purpose of this folder is to ease comparing the performance of the current rust version vs. the [original version](https://github.com/google/zopfli/blob/master/README).

The measurement is done using [hyperfine](https://github.com/sharkdp/hyperfine).

## How to

To bench, run: (requires `git` and [zopfli requirements](https://github.com/google/zopfli/blob/master/README))

```
$ cd test/perf && ./bench.sh
```
