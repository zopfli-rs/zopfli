/*
Copyright 2016 Google Inc. All Rights Reserved.

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
Utilities for using the lz77 symbols of the deflate spec.
*/

#ifndef ZOPFLI_SYMBOLS_H_
#define ZOPFLI_SYMBOLS_H_

/* __has_builtin available in clang */
#ifdef __has_builtin
# if __has_builtin(__builtin_clz)
#   define ZOPFLI_HAS_BUILTIN_CLZ
# endif
/* __builtin_clz available beginning with GCC 3.4 */
#elif __GNUC__ * 100 + __GNUC_MINOR__ >= 304
# define ZOPFLI_HAS_BUILTIN_CLZ
#endif

extern int ZopfliGetDistExtraBits(int dist);
extern int ZopfliGetDistExtraBitsValue(int dist);
extern int ZopfliGetDistSymbol(int dist);
extern int ZopfliGetLengthExtraBits(int l);
extern int ZopfliGetLengthExtraBitsValue(int l);
extern int ZopfliGetLengthSymbol(int l);
extern int ZopfliGetLengthSymbolExtraBits(int s);

/* Gets the amount of extra bits for the given distance symbol. */
static int ZopfliGetDistSymbolExtraBits(int s) {
  static const int table[30] = {
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8,
    9, 9, 10, 10, 11, 11, 12, 12, 13, 13
  };
  return table[s];
}

#endif  /* ZOPFLI_SYMBOLS_H_ */
