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

#include "deflate.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "blocksplitter.h"
#include "squeeze.h"
#include "symbols.h"
#include "tree.h"

/*
bp = bitpointer, always in range [0, 7].
The outsize is number of necessary bytes to encode the bits.
Given the value of bp and the amount of bytes, the amount of bits represented
is not simply bytesize * 8 + bp because even representing one bit requires a
whole byte. It is: (bp == 0) ? (bytesize * 8) : ((bytesize - 1) * 8 + bp)
*/
/* Actually writes to bp, out, outsize */
void AddBit(int bit,
                   unsigned char* bp, unsigned char** out, size_t* outsize) {
  if (*bp == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
  (*out)[*outsize - 1] |= bit << *bp;
  *bp = (*bp + 1) & 7;
}

/* Actually writes to bp, out, outsize */
void AddBits(unsigned symbol, unsigned length,
                    unsigned char* bp, unsigned char** out, size_t* outsize) {
  /* TODO(lode): make more efficient (add more bits at once). */
  unsigned i;
  for (i = 0; i < length; i++) {
    unsigned bit = (symbol >> i) & 1;
    if (*bp == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
    (*out)[*outsize - 1] |= bit << *bp;
    *bp = (*bp + 1) & 7;
  }
}

/*
Adds bits, like AddBits, but the order is inverted. The deflate specification
uses both orders in one standard.
*/
/* Actually writes to bp, out, outsize */
void AddHuffmanBits(unsigned symbol, unsigned length,
                           unsigned char* bp, unsigned char** out,
                           size_t* outsize) {
  /* TODO(lode): make more efficient (add more bits at once). */
  unsigned i;
  for (i = 0; i < length; i++) {
    unsigned bit = (symbol >> (length - i - 1)) & 1;
    if (*bp == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
    (*out)[*outsize - 1] |= bit << *bp;
    *bp = (*bp + 1) & 7;
  }
}

extern void PatchDistanceCodesForBuggyDecoders(unsigned* d_lengths);

extern void AddDynamicTree(const unsigned* ll_lengths, const unsigned* d_lengths, unsigned char* bp, unsigned char** out, size_t* outsize);

extern size_t CalculateTreeSize(const unsigned* ll_lengths, const unsigned* d_lengths);

/*
Adds all lit/len and dist codes from the lists as huffman symbols. Does not add
end code 256. expected_data_size is the uncompressed block size, used for
assert, but you can set it to 0 to not do the assertion.
*/
/* Passthrough */
void AddLZ77Data(const ZopfliLZ77Store* lz77,
                        size_t lstart, size_t lend,
                        size_t expected_data_size,
                        const unsigned* ll_symbols, const unsigned* ll_lengths,
                        const unsigned* d_symbols, const unsigned* d_lengths,
                        unsigned char* bp,
                        unsigned char** out, size_t* outsize) {
  size_t testlength = 0;
  size_t i;

  for (i = lstart; i < lend; i++) {
    unsigned dist = lz77->dists[i];
    unsigned litlen = lz77->litlens[i];
    if (dist == 0) {
      assert(litlen < 256);
      assert(ll_lengths[litlen] > 0);
      AddHuffmanBits(ll_symbols[litlen], ll_lengths[litlen], bp, out, outsize);
      testlength++;
    } else {
      unsigned lls = ZopfliGetLengthSymbol(litlen);
      unsigned ds = ZopfliGetDistSymbol(dist);
      assert(litlen >= 3 && litlen <= 288);
      assert(ll_lengths[lls] > 0);
      assert(d_lengths[ds] > 0);
      AddHuffmanBits(ll_symbols[lls], ll_lengths[lls], bp, out, outsize);
      AddBits(ZopfliGetLengthExtraBitsValue(litlen),
              ZopfliGetLengthExtraBits(litlen),
              bp, out, outsize);
      AddHuffmanBits(d_symbols[ds], d_lengths[ds], bp, out, outsize);
      AddBits(ZopfliGetDistExtraBitsValue(dist),
              ZopfliGetDistExtraBits(dist),
              bp, out, outsize);
      testlength += litlen;
    }
  }
  assert(expected_data_size == 0 || testlength == expected_data_size);
}

extern void GetFixedTree(unsigned* ll_lengths, unsigned* d_lengths);
extern size_t CalculateBlockSymbolSizeGivenCounts(const size_t* ll_counts, const size_t* d_counts, const unsigned* ll_lengths, const unsigned* d_lengths, const ZopfliLZ77Store* lz77, size_t lstart, size_t lend);
extern size_t CalculateBlockSymbolSize(const unsigned* ll_lengths, const unsigned* d_lengths, const ZopfliLZ77Store* lz77, size_t lstart, size_t lend);
extern void OptimizeHuffmanForRle(int length, size_t* counts);

extern double TryOptimizeHuffmanForRle(const ZopfliLZ77Store* lz77, size_t lstart, size_t lend, const size_t* ll_counts, const size_t* d_counts, unsigned* ll_lengths, unsigned* d_lengths);

extern double ZopfliCalculateBlockSize(const ZopfliLZ77Store* lz77, size_t lstart, size_t lend, int btype);

double ZopfliCalculateBlockSizeAutoType(const ZopfliLZ77Store* lz77,
                                        size_t lstart, size_t lend) {
  double uncompressedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 0);
  /* Don't do the expensive fixed cost calculation for larger blocks that are
     unlikely to use it. */
  double fixedcost = (lz77->size > 1000) ?
      uncompressedcost : ZopfliCalculateBlockSize(lz77, lstart, lend, 1);
  double dyncost = ZopfliCalculateBlockSize(lz77, lstart, lend, 2);
  return (uncompressedcost < fixedcost && uncompressedcost < dyncost)
      ? uncompressedcost
      : (fixedcost < dyncost ? fixedcost : dyncost);
}

/* Since an uncompressed block can be max 65535 in size, it actually adds
multible blocks if needed. */
/* Actually writes to bp, out, outsize
 AND passthrough */
void AddNonCompressedBlock(const ZopfliOptions* options, int final,
                                  const unsigned char* in, size_t instart,
                                  size_t inend,
                                  unsigned char* bp,
                                  unsigned char** out, size_t* outsize) {
  size_t pos = instart;
  (void)options;
  for (;;) {
    size_t i;
    unsigned short blocksize = 65535;
    unsigned short nlen;
    int currentfinal;

    if (pos + blocksize > inend) blocksize = inend - pos;
    currentfinal = pos + blocksize >= inend;

    nlen = ~blocksize;

    AddBit(final && currentfinal, bp, out, outsize);
    /* BTYPE 00 */
    AddBit(0, bp, out, outsize);
    AddBit(0, bp, out, outsize);

    /* Any bits of input up to the next byte boundary are ignored. */
    *bp = 0;

    ZOPFLI_APPEND_DATA(blocksize % 256, out, outsize);
    ZOPFLI_APPEND_DATA((blocksize / 256) % 256, out, outsize);
    ZOPFLI_APPEND_DATA(nlen % 256, out, outsize);
    ZOPFLI_APPEND_DATA((nlen / 256) % 256, out, outsize);

    for (i = 0; i < blocksize; i++) {
      ZOPFLI_APPEND_DATA(in[pos + i], out, outsize);
    }

    if (currentfinal) break;
    pos += blocksize;
  }
}

extern void AddLZ77Block(const ZopfliOptions* options, int btype, int final, const unsigned char* in, const ZopfliLZ77Store* lz77, size_t lstart, size_t lend, size_t expected_data_size, unsigned char* bp, unsigned char** out, size_t* outsize);

/* Passthrough */
static void AddLZ77BlockAutoType(const ZopfliOptions* options, int final,
                                 const unsigned char* in,
                                 const ZopfliLZ77Store* lz77,
                                 size_t lstart, size_t lend,
                                 size_t expected_data_size,
                                 unsigned char* bp,
                                 unsigned char** out, size_t* outsize) {
  double uncompressedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 0);
  double fixedcost = ZopfliCalculateBlockSize(lz77, lstart, lend, 1);
  double dyncost = ZopfliCalculateBlockSize(lz77, lstart, lend, 2);

  /* Whether to perform the expensive calculation of creating an optimal block
  with fixed huffman tree to check if smaller. Only do this for small blocks or
  blocks which already are pretty good with fixed huffman tree. */
  int expensivefixed = (lz77->size < 1000) || fixedcost <= dyncost * 1.1;

  ZopfliLZ77Store fixedstore;
  if (lstart == lend) {
    /* Smallest empty block is represented by fixed block */
    AddBits(final, 1, bp, out, outsize);
    AddBits(1, 2, bp, out, outsize);  /* btype 01 */
    AddBits(0, 7, bp, out, outsize);  /* end symbol has code 0000000 */
    return;
  }
  ZopfliInitLZ77Store(&fixedstore);
  if (expensivefixed) {
    /* Recalculate the LZ77 with ZopfliLZ77OptimalFixed */
    size_t instart = lz77->pos[lstart];
    size_t inend = instart + ZopfliLZ77GetByteRange(lz77, lstart, lend);

    ZopfliBlockState s;
    ZopfliInitBlockState(options, instart, inend, 1, &s);
    ZopfliLZ77OptimalFixed(&s, in, instart, inend, &fixedstore);
    fixedcost = ZopfliCalculateBlockSize(&fixedstore, 0, fixedstore.size, 1);
    ZopfliCleanBlockState(&s);
  }

  if (uncompressedcost < fixedcost && uncompressedcost < dyncost) {
    AddLZ77Block(options, 0, final, in, lz77, lstart, lend,
                 expected_data_size, bp, out, outsize);
  } else if (fixedcost < dyncost) {
    if (expensivefixed) {
      AddLZ77Block(options, 1, final, in, &fixedstore, 0, fixedstore.size,
                   expected_data_size, bp, out, outsize);
    } else {
      AddLZ77Block(options, 1, final, in, lz77, lstart, lend,
                   expected_data_size, bp, out, outsize);
    }
  } else {
    AddLZ77Block(options, 2, final, in, lz77, lstart, lend,
                 expected_data_size, bp, out, outsize);
  }

  ZopfliCleanLZ77Store(&fixedstore);
}

void AddAllBlocks(size_t npoints, size_t* splitpoints, ZopfliLZ77Store lz77, const ZopfliOptions* options, int final, const unsigned char* in, unsigned char* bp, unsigned char** out, size_t* outsize) {
    size_t i;

    for (i = 0; i <= npoints; i++) {
      size_t start = i == 0 ? 0 : splitpoints[i - 1];
      size_t end = i == npoints ? lz77.size : splitpoints[i];
      AddLZ77BlockAutoType(options, i == npoints && final,
                           in, &lz77, start, end, 0,
                           bp, out, outsize);
    }
}

void BlocksplitAttempt(const ZopfliOptions* options, int final,
                       const unsigned char* in, size_t instart, size_t inend,
                       unsigned char* bp, unsigned char** out,
                       size_t* outsize) {
    size_t i;
    double totalcost = 0;
    ZopfliLZ77Store lz77;
    /* byte coordinates rather than lz77 index */
    size_t* splitpoints_uncompressed = 0;
    size_t npoints = 0;
    size_t* splitpoints = 0;

    ZopfliBlockSplit(options, in, instart, inend,
                   options->blocksplittingmax,
                   &splitpoints_uncompressed, &npoints);
    splitpoints = (size_t*)malloc(sizeof(*splitpoints) * npoints);

    ZopfliInitLZ77Store(&lz77);

    for (i = 0; i <= npoints; i++) {
      size_t start = i == 0 ? instart : splitpoints_uncompressed[i - 1];
      size_t end = i == npoints ? inend : splitpoints_uncompressed[i];
      ZopfliBlockState s;
      ZopfliLZ77Store store;
      ZopfliInitLZ77Store(&store);
      ZopfliInitBlockState(options, start, end, 1, &s);
      ZopfliLZ77Optimal(&s, in, start, end, options->numiterations, &store);
      totalcost += ZopfliCalculateBlockSizeAutoType(&store, 0, store.size);

      ZopfliAppendLZ77Store(&store, &lz77);
      if (i < npoints) splitpoints[i] = lz77.size;

      ZopfliCleanBlockState(&s);
      ZopfliCleanLZ77Store(&store);
    }

    /* Second block splitting attempt */
    if (npoints > 1) {
      size_t* splitpoints2 = 0;
      size_t npoints2 = 0;
      double totalcost2 = 0;

      ZopfliBlockSplitLZ77(options, &lz77,
                           options->blocksplittingmax, &splitpoints2, &npoints2);

      for (i = 0; i <= npoints2; i++) {
        size_t start = i == 0 ? 0 : splitpoints2[i - 1];
        size_t end = i == npoints2 ? lz77.size : splitpoints2[i];
        totalcost2 += ZopfliCalculateBlockSizeAutoType(&lz77, start, end);
      }

      if (totalcost2 < totalcost) {
        free(splitpoints);
        splitpoints = splitpoints2;
        npoints = npoints2;
      } else {
        free(splitpoints2);
      }
    }

    AddAllBlocks(npoints, splitpoints, lz77, options, final, in, bp, out, outsize);

    ZopfliCleanLZ77Store(&lz77);
    free(splitpoints);
    free(splitpoints_uncompressed);
}

/*
Deflate a part, to allow ZopfliDeflate() to use multiple master blocks if
needed.
It is possible to call this function multiple times in a row, shifting
instart and inend to next bytes of the data. If instart is larger than 0, then
previous bytes are used as the initial dictionary for LZ77.
This function will usually output multiple deflate blocks. If final is 1, then
the final bit will be set on the last block.
*/
/* Passthrough of bp/out/outsize
 Allocates npoints/splitpoints/splitpoints_uncompressed */
void ZopfliDeflatePart(const ZopfliOptions* options, int btype, int final,
                       const unsigned char* in, size_t instart, size_t inend,
                       unsigned char* bp, unsigned char** out,
                       size_t* outsize) {

  /* If btype=2 is specified, it tries all block types. If a lesser btype is
  given, then however it forces that one. Neither of the lesser types needs
  block splitting as they have no dynamic huffman trees. */
  if (btype == 0) {
    AddNonCompressedBlock(options, final, in, instart, inend, bp, out, outsize);
    return;
  } else if (btype == 1) {
    ZopfliLZ77Store store;
    ZopfliBlockState s;
    ZopfliInitLZ77Store(&store);
    ZopfliInitBlockState(options, instart, inend, 1, &s);

    ZopfliLZ77OptimalFixed(&s, in, instart, inend, &store);
    AddLZ77Block(options, btype, final, in, &store, 0, store.size, 0,
                 bp, out, outsize);

    ZopfliCleanBlockState(&s);
    ZopfliCleanLZ77Store(&store);
    return;
  }

  BlocksplitAttempt(options, final, in, instart, inend, bp, out, outsize);

}

/* Passthrough */
void ZopfliDeflate(const ZopfliOptions* options, int btype, int final,
                   const unsigned char* in, size_t insize,
                   unsigned char* bp, unsigned char** out, size_t* outsize) {
 size_t offset = *outsize;
#if ZOPFLI_MASTER_BLOCK_SIZE == 0
  ZopfliDeflatePart(options, btype, final, in, 0, insize, bp, out, outsize);
#else
  size_t i = 0;
  do {
    int masterfinal = (i + ZOPFLI_MASTER_BLOCK_SIZE >= insize);
    int final2 = final && masterfinal;
    size_t size = masterfinal ? insize - i : ZOPFLI_MASTER_BLOCK_SIZE;
    ZopfliDeflatePart(options, btype, final2,
                      in, i, i + size, bp, out, outsize);
    i += size;
  } while (i < insize);
#endif
  if (options->verbose) {
    fprintf(stderr,
            "Original Size: %lu, Deflate: %lu, Compression: %f%% Removed\n",
            (unsigned long)insize, (unsigned long)(*outsize - offset),
            100.0 * (double)(insize - (*outsize - offset)) / (double)insize);
  }
}
