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

#include "blocksplitter.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "deflate.h"
#include "squeeze.h"
#include "util.h"

extern double SplitCost(size_t i, void* context);

/* Actually writes to out, outsize which are splitpoints, npoints */
void AddSorted(size_t value, size_t** out, size_t* outsize) {
  size_t i;
  ZOPFLI_APPEND_DATA(value, out, outsize);
  for (i = 0; i + 1 < *outsize; i++) {
    if ((*out)[i] > value) {
      size_t j;
      for (j = *outsize - 1; j > i; j--) {
        (*out)[j] = (*out)[j - 1];
      }
      (*out)[i] = value;
      break;
    }
  }
}

extern void ZopfliBlockSplitLZ77(const ZopfliOptions* options, const ZopfliLZ77Store* lz77, size_t maxblocks, size_t** splitpoints, size_t* npoints);

void ZopfliAppendDataSizeT(size_t value, size_t** data, size_t* size) {
    ZOPFLI_APPEND_DATA(value, data, size);
}

extern void ZopfliBlockSplit(const ZopfliOptions* options, const unsigned char* in, size_t instart, size_t inend, size_t maxblocks, size_t** splitpoints, size_t* npoints);
