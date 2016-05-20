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
Bounded package merge algorithm, based on the paper
"A Fast and Space-Economical Algorithm for Length-Limited Coding
Jyrki Katajainen, Alistair Moffat, Andrew Turpin".
*/

#include "katajainen.h"
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

typedef struct Node Node;
typedef struct Leaf Leaf;

/*
Nodes forming chains.
*/
struct Node {
  size_t weight;  /* Total weight (symbol count) of this chain. */
  Node* tail;  /* Previous node(s) of this chain, or 0 if none. */
  int count;  /* Leaf symbol index, or number of leaves before this chain. */
};

struct Leaf {
  size_t weight;
  int count;
};

/*
Memory pool for nodes.
*/
typedef struct NodePool {
  Node* next;  /* Pointer to a free node in the pool. */
} NodePool;

extern void InitNode(size_t weight, int count, Node* tail, Node* node);

/*
Performs a Boundary Package-Merge step. Puts a new chain in the given list. The
new chain is, depending on the weights, a leaf or a combination of two chains
from the previous list.
lists: The lists of chains.
maxbits: Number of lists.
leaves: The leaves, one per symbol.
numsymbols: Number of leaves.
pool: the node memory pool.
index: The index of the list in which a new chain or leaf is required.
*/
static void BoundaryPM(Node* (*lists)[2], Leaf* leaves, int numsymbols,
                       NodePool* pool, int index) {
  Node* newchain;
  Node* oldchain;
  int lastcount = lists[index][1]->count;  /* Count of last chain of list. */

  if (index == 0 && lastcount >= numsymbols) return;

  newchain = pool->next++;
  oldchain = lists[index][1];

  /* These are set up before the recursive calls below, so that there is a list
  pointing to the new node, to let the garbage collection know it's in use. */
  lists[index][0] = oldchain;
  lists[index][1] = newchain;

  if (index == 0) {
    /* New leaf node in list 0. */
    InitNode(leaves[lastcount].weight, lastcount + 1, 0, newchain);
  } else {
    size_t sum = lists[index - 1][0]->weight + lists[index - 1][1]->weight;
    if (lastcount < numsymbols && sum > leaves[lastcount].weight) {
      /* New leaf inserted in list, so count is incremented. */
      InitNode(leaves[lastcount].weight, lastcount + 1, oldchain->tail,
          newchain);
    } else {
      InitNode(sum, lastcount, lists[index - 1][1], newchain);
      /* Two lookahead chains of previous list used up, create new ones. */
      BoundaryPM(lists, leaves, numsymbols, pool, index - 1);
      BoundaryPM(lists, leaves, numsymbols, pool, index - 1);
    }
  }
}

static void BoundaryPMFinal(Node* (*lists)[2],
    Leaf* leaves, int numsymbols, NodePool* pool, int index) {
  int lastcount = lists[index][1]->count;  /* Count of last chain of list. */

  size_t sum = lists[index - 1][0]->weight + lists[index - 1][1]->weight;

  if (lastcount < numsymbols && sum > leaves[lastcount].weight) {
    Node* newchain = pool->next;
    Node* oldchain = lists[index][1]->tail;

    lists[index][1] = newchain;
    newchain->count = lastcount + 1;
    newchain->tail = oldchain;
  } else {
    lists[index][1]->tail = lists[index - 1][1];
  }
}

/*
Initializes each list with as lookahead chains the two leaves with lowest
weights.
*/
static void InitLists(
    NodePool* pool, const Leaf* leaves, int maxbits, Node* (*lists)[2]) {
  int i;
  Node* node0 = pool->next++;
  Node* node1 = pool->next++;
  InitNode(leaves[0].weight, 1, 0, node0);
  InitNode(leaves[1].weight, 2, 0, node1);
  for (i = 0; i < maxbits; i++) {
    lists[i][0] = node0;
    lists[i][1] = node1;
  }
}

extern void ExtractBitLengths(Node* chain, Leaf* leaves, unsigned* bitlengths);

/*
Comparator for sorting the leaves. Has the function signature for qsort.
*/
static int LeafComparator(const void* a, const void* b) {
  return ((const Leaf*)a)->weight - ((const Leaf*)b)->weight;
}

extern int ZopfliLengthLimitedCodeLengths(const size_t* frequencies, int n, int maxbits, unsigned* bitlengths);
