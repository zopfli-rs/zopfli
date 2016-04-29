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
The hash for ZopfliFindLongestMatch of lz77.c.
*/

#ifndef ZOPFLI_HASH_H_
#define ZOPFLI_HASH_H_

#include "util.h"

typedef struct ZopfliHash ZopfliHash;

ZopfliHash* ZopfliInitHash(size_t window_size);

/* Resets all fields of ZopfliHash. */
void ZopfliResetHash(size_t window_size, ZopfliHash* h);

/* Frees ZopfliHash memory. */
void ZopfliCleanHash(ZopfliHash* h);

/*
Updates the hash values based on the current position in the array. All calls
to this must be made for consecutive bytes.
*/
void ZopfliUpdateHash(const unsigned char* array, size_t pos, size_t end,
                      ZopfliHash* h);

/*
Prepopulates hash:
Fills in the initial values in the hash, before ZopfliUpdateHash can be used
correctly.
*/
void ZopfliWarmupHash(const unsigned char* array, size_t pos, size_t end,
                      ZopfliHash* h);

#endif  /* ZOPFLI_HASH_H_ */
