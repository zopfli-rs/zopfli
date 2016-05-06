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

extern void ZopfliInitLZ77Store(const unsigned char* data, ZopfliLZ77Store* store);

void ZopfliCleanLZ77Store(ZopfliLZ77Store* store) {
  free(store->litlens);
  free(store->dists);
  free(store->pos);
  free(store->ll_symbol);
  free(store->d_symbol);
  free(store->ll_counts);
  free(store->d_counts);
}

static size_t CeilDiv(size_t a, size_t b) {
  return (a + b - 1) / b;
}

void ZopfliCopyLZ77Store(
    const ZopfliLZ77Store* source, ZopfliLZ77Store* dest) {
  size_t i;
  size_t llsize = ZOPFLI_NUM_LL * CeilDiv(source->size, ZOPFLI_NUM_LL);
  size_t dsize = ZOPFLI_NUM_D * CeilDiv(source->size, ZOPFLI_NUM_D);
  ZopfliCleanLZ77Store(dest);
  ZopfliInitLZ77Store(source->data, dest);
  dest->litlens =
      (unsigned short*)malloc(sizeof(*dest->litlens) * source->size);
  dest->dists = (unsigned short*)malloc(sizeof(*dest->dists) * source->size);
  dest->pos = (size_t*)malloc(sizeof(*dest->pos) * source->size);
  dest->ll_symbol =
      (unsigned short*)malloc(sizeof(*dest->ll_symbol) * source->size);
  dest->d_symbol =
      (unsigned short*)malloc(sizeof(*dest->d_symbol) * source->size);
  dest->ll_counts = (size_t*)malloc(sizeof(*dest->ll_counts) * llsize);
  dest->d_counts = (size_t*)malloc(sizeof(*dest->d_counts) * dsize);

  /* Allocation failed. */
  if (!dest->litlens || !dest->dists) exit(-1);
  if (!dest->pos) exit(-1);
  if (!dest->ll_symbol || !dest->d_symbol) exit(-1);
  if (!dest->ll_counts || !dest->d_counts) exit(-1);

  dest->size = source->size;
  for (i = 0; i < source->size; i++) {
    dest->litlens[i] = source->litlens[i];
    dest->dists[i] = source->dists[i];
    dest->pos[i] = source->pos[i];
    dest->ll_symbol[i] = source->ll_symbol[i];
    dest->d_symbol[i] = source->d_symbol[i];
  }
  for (i = 0; i < llsize; i++) {
    dest->ll_counts[i] = source->ll_counts[i];
  }
  for (i = 0; i < dsize; i++) {
    dest->d_counts[i] = source->d_counts[i];
  }
}

/*
Appends the length and distance to the LZ77 arrays of the ZopfliLZ77Store.
context must be a ZopfliLZ77Store*.
*/
void ZopfliStoreLitLenDist(unsigned short length, unsigned short dist,
                           size_t pos, ZopfliLZ77Store* store) {
  size_t i;
  /* Needed for using ZOPFLI_APPEND_DATA multiple times. */
  size_t origsize = store->size;
  size_t llstart = ZOPFLI_NUM_LL * (origsize / ZOPFLI_NUM_LL);
  size_t dstart = ZOPFLI_NUM_D * (origsize / ZOPFLI_NUM_D);

  /* Everytime the index wraps around, a new cumulative histogram is made: we're
  keeping one histogram value per LZ77 symbol rather than a full histogram for
  each to save memory. */
  if (origsize % ZOPFLI_NUM_LL == 0) {
    size_t llsize = origsize;
    for (i = 0; i < ZOPFLI_NUM_LL; i++) {
      ZOPFLI_APPEND_DATA(
          origsize == 0 ? 0 : store->ll_counts[origsize - ZOPFLI_NUM_LL + i],
          &store->ll_counts, &llsize);
    }
  }
  if (origsize % ZOPFLI_NUM_D == 0) {
    size_t dsize = origsize;
    for (i = 0; i < ZOPFLI_NUM_D; i++) {
      ZOPFLI_APPEND_DATA(
          origsize == 0 ? 0 : store->d_counts[origsize - ZOPFLI_NUM_D + i],
          &store->d_counts, &dsize);
    }
  }

  ZOPFLI_APPEND_DATA(length, &store->litlens, &store->size);
  store->size = origsize;
  ZOPFLI_APPEND_DATA(dist, &store->dists, &store->size);
  store->size = origsize;
  ZOPFLI_APPEND_DATA(pos, &store->pos, &store->size);
  assert(length < 259);

  if (dist == 0) {
    store->size = origsize;
    ZOPFLI_APPEND_DATA(length, &store->ll_symbol, &store->size);
    store->size = origsize;
    ZOPFLI_APPEND_DATA(0, &store->d_symbol, &store->size);
    store->ll_counts[llstart + length]++;
  } else {
    store->size = origsize;
    ZOPFLI_APPEND_DATA(ZopfliGetLengthSymbol(length),
                       &store->ll_symbol, &store->size);
    store->size = origsize;
    ZOPFLI_APPEND_DATA(ZopfliGetDistSymbol(dist),
                       &store->d_symbol, &store->size);
    store->ll_counts[llstart + ZopfliGetLengthSymbol(length)]++;
    store->d_counts[dstart + ZopfliGetDistSymbol(dist)]++;
  }
}

void ZopfliAppendLZ77Store(const ZopfliLZ77Store* store,
                           ZopfliLZ77Store* target) {
  size_t i;
  for (i = 0; i < store->size; i++) {
    ZopfliStoreLitLenDist(store->litlens[i], store->dists[i],
                          store->pos[i], target);
  }
}

size_t ZopfliLZ77GetByteRange(const ZopfliLZ77Store* lz77,
                              size_t lstart, size_t lend) {
  size_t l = lend - 1;
  if (lstart == lend) return 0;
  return lz77->pos[l] + ((lz77->dists[l] == 0) ?
      1 : lz77->litlens[l]) - lz77->pos[lstart];
}

static void ZopfliLZ77GetHistogramAt(const ZopfliLZ77Store* lz77, size_t lpos,
                                     size_t* ll_counts, size_t* d_counts) {
  /* The real histogram is created by using the histogram for this chunk, but
  all superfluous values of this chunk subtracted. */
  size_t llpos = ZOPFLI_NUM_LL * (lpos / ZOPFLI_NUM_LL);
  size_t dpos = ZOPFLI_NUM_D * (lpos / ZOPFLI_NUM_D);
  size_t i;
  for (i = 0; i < ZOPFLI_NUM_LL; i++) {
    ll_counts[i] = lz77->ll_counts[llpos + i];
  }
  for (i = lpos + 1; i < llpos + ZOPFLI_NUM_LL && i < lz77->size; i++) {
    ll_counts[lz77->ll_symbol[i]]--;
  }
  for (i = 0; i < ZOPFLI_NUM_D; i++) {
    d_counts[i] = lz77->d_counts[dpos + i];
  }
  for (i = lpos + 1; i < dpos + ZOPFLI_NUM_D && i < lz77->size; i++) {
    if (lz77->dists[i] != 0) d_counts[lz77->d_symbol[i]]--;
  }
}

void ZopfliLZ77GetHistogram(const ZopfliLZ77Store* lz77,
                           size_t lstart, size_t lend,
                           size_t* ll_counts, size_t* d_counts) {
  size_t i;
  if (lstart + ZOPFLI_NUM_LL * 3 > lend) {
    memset(ll_counts, 0, sizeof(*ll_counts) * ZOPFLI_NUM_LL);
    memset(d_counts, 0, sizeof(*d_counts) * ZOPFLI_NUM_D);
    for (i = lstart; i < lend; i++) {
      ll_counts[lz77->ll_symbol[i]]++;
      if (lz77->dists[i] != 0) d_counts[lz77->d_symbol[i]]++;
    }
  } else {
    /* Subtract the cumulative histograms at the end and the start to get the
    histogram for this range. */
    ZopfliLZ77GetHistogramAt(lz77, lend - 1, ll_counts, d_counts);
    if (lstart > 0) {
      size_t ll_counts2[ZOPFLI_NUM_LL];
      size_t d_counts2[ZOPFLI_NUM_D];
      ZopfliLZ77GetHistogramAt(lz77, lstart - 1, ll_counts2, d_counts2);

      for (i = 0; i < ZOPFLI_NUM_LL; i++) {
        ll_counts[i] -= ll_counts2[i];
      }
      for (i = 0; i < ZOPFLI_NUM_D; i++) {
        d_counts[i] -= d_counts2[i];
      }
    }
  }
}

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

extern int GetLengthScore(int length, int distance);

void ZopfliVerifyLenDist(const unsigned char* data, size_t datasize, size_t pos,
                         unsigned short dist, unsigned short length) {

  /* TODO(lode): make this only run in a debug compile, it's for assert only. */
  size_t i;

  assert(pos + length <= datasize);
  for (i = 0; i < length; i++) {
    if (data[pos - dist + i] != data[pos + i]) {
      assert(data[pos - dist + i] == data[pos + i]);
      break;
    }
  }
}

extern unsigned short ZopfliHashSameAt(const ZopfliHash* h, size_t index);

extern LongestMatch ZopfliFindLongestMatch(ZopfliBlockState* s, const ZopfliHash* h, const unsigned char* array, size_t pos, size_t size, size_t limit, unsigned short* sublen);

typedef struct lz77_store_S lz77_store_t;
extern lz77_store_t * lz77_store_from_c(ZopfliLZ77Store *);
extern void lz77_store_lit_len_dist(lz77_store_t *, unsigned short length, unsigned short dist, size_t pos);
extern void lz77_store_result(lz77_store_t *, ZopfliLZ77Store *);

void ZopfliLZ77Greedy(ZopfliBlockState* s, const unsigned char* in,
                      size_t instart, size_t inend,
                      ZopfliLZ77Store* store, ZopfliHash* h) {

  lz77_store_t *rust_store = lz77_store_from_c(store);

  size_t i = 0, j;
  unsigned short leng;
  unsigned short dist;
  int lengthscore;
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;
  unsigned short dummysublen[259];

  LongestMatch longest_match;

#ifdef ZOPFLI_LAZY_MATCHING
  /* Lazy matching. */
  unsigned prev_length = 0;
  unsigned prev_match = 0;
  int prevlengthscore;
  int match_available = 0;
#endif

  if (instart == inend) return;

  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) {
    ZopfliUpdateHash(in, i, inend, h);
  }

  for (i = instart; i < inend; i++) {
    ZopfliUpdateHash(in, i, inend, h);

    longest_match = ZopfliFindLongestMatch(s, h, in, i, inend, ZOPFLI_MAX_MATCH, dummysublen);
    dist = longest_match.distance;
    leng = longest_match.length;
    lengthscore = GetLengthScore(leng, dist);

#ifdef ZOPFLI_LAZY_MATCHING
    /* Lazy matching. */
    prevlengthscore = GetLengthScore(prev_length, prev_match);
    if (match_available) {
      match_available = 0;
      if (lengthscore > prevlengthscore + 1) {
        lz77_store_lit_len_dist(rust_store, in[i - 1], 0, i - 1);
        if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH) {
          match_available = 1;
          prev_length = leng;
          prev_match = dist;
          continue;
        }
      } else {
        /* Add previous to output. */
        leng = prev_length;
        dist = prev_match;
        lengthscore = prevlengthscore;
        /* Add to output. */
        ZopfliVerifyLenDist(in, inend, i - 1, dist, leng);
        lz77_store_lit_len_dist(rust_store, leng, dist, i - 1);
        for (j = 2; j < leng; j++) {
          assert(i < inend);
          i++;
          ZopfliUpdateHash(in, i, inend, h);
        }
        continue;
      }
    }
    else if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH) {
      match_available = 1;
      prev_length = leng;
      prev_match = dist;
      continue;
    }
    /* End of lazy matching. */
#endif

    /* Add to output. */
    if (lengthscore >= ZOPFLI_MIN_MATCH) {
      ZopfliVerifyLenDist(in, inend, i, dist, leng);
      lz77_store_lit_len_dist(rust_store, leng, dist, i);
    } else {
      leng = 1;
      lz77_store_lit_len_dist(rust_store, in[i], 0, i);
    }
    for (j = 1; j < leng; j++) {
      assert(i < inend);
      i++;
      ZopfliUpdateHash(in, i, inend, h);
    }
  }

  lz77_store_result(rust_store, store);
}
