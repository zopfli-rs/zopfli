## What kind of benchmarks?

The purpose of this, is to bench the current rust version vs. the [original version](https://github.com/google/zopfli/blob/master/README).

The measurement is done using [hyperfine](https://github.com/sharkdp/hyperfine).

## How to

To bench, run: (requires `git` and [zopfli requirements](https://github.com/google/zopfli/blob/master/README))

```
$ cd benchmark-builds && ./bench.sh
```
