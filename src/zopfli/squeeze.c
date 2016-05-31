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

extern void ZopfliLZ77Optimal(ZopfliBlockState *s, const unsigned char* in, size_t instart, size_t inend, int numiterations, ZopfliLZ77Store* store);

extern void ZopfliLZ77OptimalFixed(ZopfliBlockState *s, const unsigned char* in, size_t instart, size_t inend, ZopfliLZ77Store* store);
