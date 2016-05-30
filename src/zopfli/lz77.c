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
#include "util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

extern void ZopfliInitLZ77Store(ZopfliLZ77Store* store);

extern void ZopfliCleanLZ77Store(ZopfliLZ77Store* store);

extern void ZopfliCopyLZ77Store(const ZopfliLZ77Store* source, ZopfliLZ77Store* dest);
extern void ZopfliAppendLZ77Store(const ZopfliLZ77Store* store, ZopfliLZ77Store* target);

extern void ZopfliInitBlockState(const ZopfliOptions* options, size_t blockstart, size_t blockend, int add_lmc, ZopfliBlockState* s);

void ZopfliCleanBlockState(ZopfliBlockState* s) {
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  if (s->lmc) {
    ZopfliCleanCache(s->lmc);
  }
#endif
}

typedef struct lz77_store_S lz77_store_t;
extern lz77_store_t * lz77_store_from_c(ZopfliLZ77Store *);

extern void ZopfliLZ77Greedy(ZopfliBlockState* s, const unsigned char* in, size_t instart, size_t inend, ZopfliLZ77Store* store, ZopfliHash* h);
