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

#include "cache.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef ZOPFLI_LONGEST_MATCH_CACHE

extern ZopfliLongestMatchCache* ZopfliInitCache(size_t blocksize);

extern void ZopfliCleanCache(ZopfliLongestMatchCache* lmc);

extern void ZopfliSublenToCache(const unsigned short* sublen, size_t pos, size_t length, ZopfliLongestMatchCache* lmc);
extern void ZopfliCacheToSublen(const ZopfliLongestMatchCache* lmc, size_t pos, size_t length, unsigned short* sublen);
extern unsigned ZopfliMaxCachedSublen(const ZopfliLongestMatchCache* lmc, size_t pos, size_t length);

#endif  /* ZOPFLI_LONGEST_MATCH_CACHE */
