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
Functions to choose good boundaries for block splitting. Deflate allows encoding
the data in multiple blocks, with a separate Huffman tree for each block. The
Huffman tree itself requires some bytes to encode, so by choosing certain
blocks, you can either hurt, or enhance compression. These functions choose good
ones that enhance it.
*/

#ifndef ZOPFLI_BLOCKSPLITTER_H_
#define ZOPFLI_BLOCKSPLITTER_H_

#include <stdlib.h>

#include "lz77.h"
#include "zopfli.h"

#endif  /* ZOPFLI_BLOCKSPLITTER_H_ */
