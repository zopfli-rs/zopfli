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

/*
Functions for basic LZ77 compression and utilities for the "squeeze" LZ77
compression.
*/

#ifndef ZOPFLI_LZ77_H_
#define ZOPFLI_LZ77_H_

#include <stdlib.h>

#include "cache.h"
#include "hash.h"
#include "zopfli.h"

/*
Some state information for compressing a block.
This is currently a bit under-used (with mainly only the longest match cache),
but is kept for easy future expansion.
*/
typedef struct ZopfliBlockState {
  const ZopfliOptions* options;

#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  /* Cache for length/distance pairs found so far. */
  ZopfliLongestMatchCache* lmc;
#endif

  /* The start (inclusive) and end (not inclusive) of the current block. */
  size_t blockstart;
  size_t blockend;
} ZopfliBlockState;

#endif  /* ZOPFLI_LZ77_H_ */
