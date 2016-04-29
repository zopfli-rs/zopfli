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

#include "hash.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define HASH_SHIFT 5
#define HASH_MASK 32767

void ZopfliAllocHash(size_t window_size, ZopfliHash* h) {
  h->head = (int*)malloc(sizeof(*h->head) * 65536);
  h->prev = (unsigned short*)malloc(sizeof(*h->prev) * window_size);
  h->hashval = (int*)malloc(sizeof(*h->hashval) * window_size);

#ifdef ZOPFLI_HASH_SAME
  h->same = (unsigned short*)malloc(sizeof(*h->same) * window_size);
#endif

#ifdef ZOPFLI_HASH_SAME_HASH
  h->head2 = (int*)malloc(sizeof(*h->head2) * 65536);
  h->prev2 = (unsigned short*)malloc(sizeof(*h->prev2) * window_size);
  h->hashval2 = (int*)malloc(sizeof(*h->hashval2) * window_size);
#endif
}

void ZopfliResetHash(size_t window_size, ZopfliHash* h) {
  size_t i;

  h->val = 0;
  for (i = 0; i < 65536; i++) {
    h->head[i] = -1;  /* -1 indicates no head so far. */
  }
  for (i = 0; i < window_size; i++) {
    h->prev[i] = i;  /* If prev[j] == j, then prev[j] is uninitialized. */
    h->hashval[i] = -1;
  }

#ifdef ZOPFLI_HASH_SAME
  for (i = 0; i < window_size; i++) {
    h->same[i] = 0;
  }
#endif

#ifdef ZOPFLI_HASH_SAME_HASH
  h->val2 = 0;
  for (i = 0; i < 65536; i++) {
    h->head2[i] = -1;
  }
  for (i = 0; i < window_size; i++) {
    h->prev2[i] = i;
    h->hashval2[i] = -1;
  }
#endif
}

void ZopfliCleanHash(ZopfliHash* h) {
  free(h->head);
  free(h->prev);
  free(h->hashval);

#ifdef ZOPFLI_HASH_SAME_HASH
  free(h->head2);
  free(h->prev2);
  free(h->hashval2);
#endif

#ifdef ZOPFLI_HASH_SAME
  free(h->same);
#endif
}

extern void UpdateHashValue(ZopfliHash* h, unsigned char c);
extern void ZopfliUpdateHash(const unsigned char* array, size_t pos, size_t end, ZopfliHash* h);
extern void ZopfliWarmupHash(const unsigned char* array, size_t pos, size_t end, ZopfliHash* h);
