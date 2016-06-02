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

#ifndef ZOPFLI_DEFLATE_H_
#define ZOPFLI_DEFLATE_H_

/*
Functions to compress according to the DEFLATE specification, using the
"squeeze" LZ77 compression backend.
*/

#include "zopfli.h"
#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddBit(int bit, unsigned char* bp, unsigned char** out, size_t* outsize);
void AddBits(unsigned symbol, unsigned length,
                    unsigned char* bp, unsigned char** out, size_t* outsize);

void AddHuffmanBits(unsigned symbol, unsigned length,
                           unsigned char* bp, unsigned char** out,
                           size_t* outsize);

void AddNonCompressedBlock(const ZopfliOptions* options, int final,
                                 const unsigned char* in, size_t instart,
                                 size_t inend,
                                 unsigned char* bp,
                                 unsigned char** out, size_t* outsize);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  /* ZOPFLI_DEFLATE_H_ */
