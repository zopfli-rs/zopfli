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

#include "tree.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "katajainen.h"
#include "util.h"

void ZopfliLengthsToSymbols(const unsigned* lengths, size_t n, unsigned maxbits,
                            unsigned* symbols) {
  size_t* bl_count = (size_t*)malloc(sizeof(size_t) * (maxbits + 1));
  size_t* next_code = (size_t*)malloc(sizeof(size_t) * (maxbits + 1));
  unsigned bits, i;
  unsigned code;

  for (i = 0; i < n; i++) {
    symbols[i] = 0;
  }

  /* 1) Count the number of codes for each code length. Let bl_count[N] be the
  number of codes of length N, N >= 1. */
  for (bits = 0; bits <= maxbits; bits++) {
    bl_count[bits] = 0;
  }
  for (i = 0; i < n; i++) {
    assert(lengths[i] <= maxbits);
    bl_count[lengths[i]]++;
  }
  /* 2) Find the numerical value of the smallest code for each code length. */
  code = 0;
  bl_count[0] = 0;
  for (bits = 1; bits <= maxbits; bits++) {
    code = (code + bl_count[bits-1]) << 1;
    next_code[bits] = code;
  }
  /* 3) Assign numerical values to all codes, using consecutive values for all
  codes of the same length with the base values determined at step 2. */
  for (i = 0;  i < n; i++) {
    unsigned len = lengths[i];
    if (len != 0) {
      symbols[i] = next_code[len];
      next_code[len]++;
    }
  }

  free(bl_count);
  free(next_code);
}

extern void ZopfliCalculateEntropy(const size_t* count, size_t n, double* bitlengths);

void ZopfliCalculateBitLengths(const size_t* count, size_t n, int maxbits,
                               unsigned* bitlengths) {
  int error = ZopfliLengthLimitedCodeLengths(count, n, maxbits, bitlengths);
  (void) error;
  assert(!error);
}
