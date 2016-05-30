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

#include "squeeze.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "blocksplitter.h"
#include "deflate.h"
#include "util.h"

typedef struct SymbolStats SymbolStats;

extern SymbolStats* symbol_stats_new();

extern void CopyStats(SymbolStats* source, SymbolStats* dest);
extern void AddWeighedStatFreqs(const SymbolStats* stats1, double w1,
                                const SymbolStats* stats2, double w2,
                                SymbolStats* result);

typedef struct RanState RanState;

extern RanState* ran_state_new();
extern void RandomizeStatFreqs(RanState* state, SymbolStats* stats);
extern void ClearStatFreqs(SymbolStats* stats);
/*
Function that calculates a cost based on a model for the given LZ77 symbol.
litlen: means literal symbol if dist is 0, length otherwise.
*/
typedef double CostModelFun(unsigned litlen, unsigned dist, void* context);

extern double GetCostStat(unsigned litlen, unsigned dist, void* context);

extern void CalculateStatistics(SymbolStats* stats);
extern void GetStatistics(const ZopfliLZ77Store* store, SymbolStats* stats);

extern void LZ77OptimalRun(ZopfliBlockState* s, const unsigned char* in, size_t instart, size_t inend, CostModelFun* costmodel, void* costcontext, ZopfliLZ77Store* store, ZopfliHash* h, float* costs);

void ZopfliLZ77Optimal(ZopfliBlockState *s,
                       const unsigned char* in, size_t instart, size_t inend,
                       int numiterations,
                       ZopfliLZ77Store* store) {
  /* Dist to get to here with smallest cost. */
  size_t blocksize = inend - instart;
  ZopfliLZ77Store currentstore;
  ZopfliHash* h = ZopfliInitHash(ZOPFLI_WINDOW_SIZE);
  SymbolStats *stats, *beststats, *laststats;
  int i;
  float* costs = (float*)malloc(sizeof(float) * (blocksize + 1));
  double cost;
  double bestcost = ZOPFLI_LARGE_FLOAT;
  double lastcost = 0;
  /* Try randomizing the costs a bit once the size stabilizes. */
  RanState *ran_state;
  int lastrandomstep = -1;

  if (!costs) exit(-1); /* Allocation failed. */

  ran_state = ran_state_new();
  stats = symbol_stats_new();
  beststats = symbol_stats_new();
  laststats = symbol_stats_new();
  ZopfliInitLZ77Store(&currentstore);

  /* Do regular deflate, then loop multiple shortest path runs, each time using
  the statistics of the previous run. */

  /* Initial run. */
  ZopfliLZ77Greedy(s, in, instart, inend, &currentstore, h);
  GetStatistics(&currentstore, stats);

  /* Repeat statistics with each time the cost model from the previous stat
  run. */
  for (i = 0; i < numiterations; i++) {
    ZopfliCleanLZ77Store(&currentstore);
    ZopfliInitLZ77Store(&currentstore);
    LZ77OptimalRun(s, in, instart, inend,
                   GetCostStat, (void*)stats,
                   &currentstore, h, costs);
    cost = ZopfliCalculateBlockSize(&currentstore, 0, currentstore.size, 2);
    if (s->options->verbose_more || (s->options->verbose && cost < bestcost)) {
      fprintf(stderr, "Iteration %d: %d bit\n", i, (int) cost);
    }
    if (cost < bestcost) {
      /* Copy to the output store. */
      ZopfliCopyLZ77Store(&currentstore, store);
      CopyStats(stats, beststats);
      bestcost = cost;
    }
    CopyStats(stats, laststats);
    ClearStatFreqs(stats);
    GetStatistics(&currentstore, stats);
    if (lastrandomstep != -1) {
      /* This makes it converge slower but better. Do it only once the
      randomness kicks in so that if the user does few iterations, it gives a
      better result sooner. */
      AddWeighedStatFreqs(stats, 1.0, laststats, 0.5, stats);
      CalculateStatistics(stats);
    }
    if (i > 5 && cost == lastcost) {
      CopyStats(beststats, stats);
      RandomizeStatFreqs(ran_state, stats);
      CalculateStatistics(stats);
      lastrandomstep = i;
    }
    lastcost = cost;
  }

  free(costs);
  ZopfliCleanLZ77Store(&currentstore);
  ZopfliCleanHash(h);
}

extern void ZopfliLZ77OptimalFixed(ZopfliBlockState *s, const unsigned char* in, size_t instart, size_t inend, ZopfliLZ77Store* store);
