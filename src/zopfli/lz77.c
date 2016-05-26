/*
Copyright 2011 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: lode.vandevenne@gmail.com (Lode Vandevenne)
Author: jyrki.alakuijala@gmail.com (Jyrki Alakuijala)
*/

#include "lz77.h"
#include "symbols.h"
#include "util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

extern void ZopfliInitLZ77Store(ZopfliLZ77Store* store);

extern void ZopfliCleanLZ77Store(ZopfliLZ77Store* store);

extern size_t CeilDiv(size_t a, size_t b);

extern void ZopfliCopyLZ77Store(const ZopfliLZ77Store* source, ZopfliLZ77Store* dest);
extern void ZopfliAppendLZ77Store(const ZopfliLZ77Store* store, ZopfliLZ77Store* target);

extern size_t ZopfliLZ77GetByteRange(const ZopfliLZ77Store* lz77, size_t lstart, size_t lend);

extern void ZopfliLZ77GetHistogram(const ZopfliLZ77Store* lz77, size_t lstart, size_t lend, size_t* ll_counts, size_t* d_counts);

void ZopfliInitBlockState(const ZopfliOptions* options,
                          size_t blockstart, size_t blockend, int add_lmc,
                          ZopfliBlockState* s) {
  s->options = options;
  s->blockstart = blockstart;
  s->blockend = blockend;
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  if (add_lmc) {
    s->lmc = ZopfliInitCache(blockend - blockstart);
  } else {
    s->lmc = 0;
  }
#endif
}

void ZopfliCleanBlockState(ZopfliBlockState* s) {
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  if (s->lmc) {
    ZopfliCleanCache(s->lmc);
  }
#endif
}

typedef struct lz77_store_S lz77_store_t;
extern lz77_store_t * lz77_store_from_c(ZopfliLZ77Store *);
extern void lz77_store_result(lz77_store_t *, ZopfliLZ77Store *);

extern void ZopfliLZ77Greedy(ZopfliBlockState* s, const unsigned char* in, size_t instart, size_t inend, ZopfliLZ77Store* store, ZopfliHash* h);
