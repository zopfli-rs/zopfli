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


extern ZopfliHash* ZopfliInitHash(size_t window_size);
extern void ZopfliResetHash(size_t window_size, ZopfliHash* h);
extern void ZopfliCleanHash(ZopfliHash* h);
extern void ZopfliUpdateHash(const unsigned char* array, size_t pos, size_t end, ZopfliHash* h);
extern void ZopfliWarmupHash(const unsigned char* array, size_t pos, size_t end, ZopfliHash* h);
