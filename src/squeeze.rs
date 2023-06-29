//! The squeeze functions do enhanced LZ77 compression by optimal parsing with a
//! cost model, rather than greedily choosing the longest length or using a single
//! step of lazy matching like regular implementations.
//!
//! Since the cost model is based on the Huffman tree that can only be calculated
//! after the LZ77 data is generated, there is a chicken and egg problem, and
//! multiple runs are done with updated cost models to converge to a better
//! solution.

use alloc::vec::Vec;
use core::{
    cmp,
    fmt::{Debug, Display, Formatter},
    ops::DerefMut,
};
use std::iter;

use genevo::{
    ga,
    genetic::{Children, Genotype, Parents},
    mutation::value::RandomGenomeMutation,
    operator::{
        prelude::{ElitistReinserter, RandomValueMutator},
        CrossoverOp, GeneticOperator,
    },
    prelude::*,
    selection::truncation::MaximizeSelector,
    simulation::State,
    termination::{StopFlag, Termination},
};
use lockfree_object_pool::LinearObjectPool;
use log::debug;
use once_cell::sync::Lazy;
use ordered_float::OrderedFloat;

use crate::{
    cache::Cache,
    deflate::{calculate_block_size, BlockType},
    hash::ZopfliHash,
    lz77::{find_longest_match, LitLen, Lz77Store, ZopfliBlockState},
    symbols::{get_dist_extra_bits, get_dist_symbol, get_length_extra_bits, get_length_symbol},
    util::{ZOPFLI_MAX_MATCH, ZOPFLI_NUM_D, ZOPFLI_NUM_LL, ZOPFLI_WINDOW_MASK, ZOPFLI_WINDOW_SIZE},
};

const K_INV_LOG2: f64 = core::f64::consts::LOG2_E; // 1.0 / log(2.0)

static LZ77_STORE_POOL: Lazy<LinearObjectPool<Lz77Store>> =
    Lazy::new(|| LinearObjectPool::new(Lz77Store::new, Lz77Store::reset));

/// Cost model which should exactly match fixed tree.
const fn get_cost_fixed(litlen: usize, dist: u16) -> f64 {
    let result = if dist == 0 {
        if litlen <= 143 {
            8
        } else {
            9
        }
    } else {
        let dbits = get_dist_extra_bits(dist);
        let lbits = get_length_extra_bits(litlen);
        let lsym = get_length_symbol(litlen);
        // Every dist symbol has length 5.
        7 + (lsym > 279) as usize + 5 + dbits + lbits
    };
    result as f64
}

/// Cost model based on symbol statistics.
fn get_cost_stat(litlen: usize, dist: u16, stats: &SymbolStats) -> f64 {
    if dist == 0 {
        stats.ll_symbols[litlen]
    } else {
        let lsym = get_length_symbol(litlen);
        let lbits = get_length_extra_bits(litlen) as f64;
        let dsym = get_dist_symbol(dist);
        let dbits = get_dist_extra_bits(dist) as f64;
        lbits + dbits + stats.ll_symbols[lsym] + stats.d_symbols[dsym]
    }
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
struct SymbolTable {
    /* The literal and length symbols. */
    litlens: [usize; ZOPFLI_NUM_LL],
    /* The 32 unique dist symbols, not the 32768 possible dists. */
    dists: [usize; ZOPFLI_NUM_D],
}

impl Default for SymbolTable {
    fn default() -> Self {
        SymbolTable {
            litlens: [0; ZOPFLI_NUM_LL],
            dists: [0; ZOPFLI_NUM_D],
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct SymbolStats {
    table: SymbolTable,
    /* Length of each lit/len symbol in bits. */
    ll_symbols: [f64; ZOPFLI_NUM_LL],
    /* Length of each dist symbol in bits. */
    d_symbols: [f64; ZOPFLI_NUM_D],
}

impl Default for SymbolStats {
    fn default() -> SymbolStats {
        SymbolStats {
            table: Default::default(),
            ll_symbols: [0.0; ZOPFLI_NUM_LL],
            d_symbols: [0.0; ZOPFLI_NUM_D],
        }
    }
}

impl SymbolStats {
    /// Calculates the entropy of each symbol, based on the counts of each symbol. The
    /// result is similar to the result of length_limited_code_lengths, but with the
    /// actual theoritical bit lengths according to the entropy. Since the resulting
    /// values are fractional, they cannot be used to encode the tree specified by
    /// DEFLATE.
    fn calculate_entropy(&mut self) {
        fn calculate_and_store_entropy(count: &[usize], bitlengths: &mut [f64]) {
            let n = count.len();

            let sum = count.iter().sum();

            let log2sum = (if sum == 0 { n } else { sum } as f64).ln() * K_INV_LOG2;

            for i in 0..n {
                // When the count of the symbol is 0, but its cost is requested anyway, it
                // means the symbol will appear at least once anyway, so give it the cost as if
                // its count is 1.
                if count[i] == 0 {
                    bitlengths[i] = log2sum;
                } else {
                    bitlengths[i] = log2sum - (count[i] as f64).ln() * K_INV_LOG2;
                }

                // Depending on compiler and architecture, the above subtraction of two
                // floating point numbers may give a negative result very close to zero
                // instead of zero (e.g. -5.973954e-17 with gcc 4.1.2 on Ubuntu 11.4). Clamp
                // it to zero. These floating point imprecisions do not affect the cost model
                // significantly so this is ok.
                if bitlengths[i] < 0.0 && bitlengths[i] > -1E-5 {
                    bitlengths[i] = 0.0;
                }
                debug_assert!(bitlengths[i] >= 0.0);
            }
        }

        calculate_and_store_entropy(&self.table.litlens, &mut self.ll_symbols);
        calculate_and_store_entropy(&self.table.dists, &mut self.d_symbols);
    }

    /// Appends the symbol statistics from the store.
    fn get_statistics(&mut self, store: &Lz77Store) {
        for &litlen in &store.litlens {
            match litlen {
                LitLen::Literal(lit) => self.table.litlens[lit as usize] += 1,
                LitLen::LengthDist(len, dist) => {
                    self.table.litlens[get_length_symbol(len as usize)] += 1;
                    self.table.dists[get_dist_symbol(dist)] += 1;
                }
            }
        }
        self.table.litlens[256] = 1; /* End symbol. */

        self.calculate_entropy();
    }
}

/// Finds the minimum possible cost this cost model can return for valid length and
/// distance symbols.
fn get_cost_model_min_cost<F: Fn(usize, u16) -> f64>(costmodel: F) -> f64 {
    let mut bestlength = 0; // length that has lowest cost in the cost model
    let mut bestdist = 0; // distance that has lowest cost in the cost model

    // Table of distances that have a different distance symbol in the deflate
    // specification. Each value is the first distance that has a new symbol. Only
    // different symbols affect the cost model so only these need to be checked.
    // See RFC 1951 section 3.2.5. Compressed blocks (length and distance codes).

    const DSYMBOLS: [u16; 30] = [
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
        2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
    ];

    let mut mincost = f64::INFINITY;
    for i in 3..259 {
        let c = costmodel(i, 1);
        if c < mincost {
            bestlength = i;
            mincost = c;
        }
    }

    mincost = f64::INFINITY;
    for dsym in DSYMBOLS {
        let c = costmodel(3, dsym);
        if c < mincost {
            bestdist = dsym;
            mincost = c;
        }
    }
    costmodel(bestlength, bestdist)
}

/// Performs the forward pass for "squeeze". Gets the most optimal length to reach
/// every byte from a previous byte, using cost calculations.
/// `s`: the `ZopfliBlockState`
/// `costmodel`: function to calculate the cost of some lit/len/dist pair.
/// `length_array`: output array of size `(inend - instart)` which will receive the best
///     length to reach this byte from a previous byte.
/// returns the cost that was, according to the `costmodel`, needed to get to the end.
fn get_best_lengths<F: Fn(usize, u16) -> f64, C: Cache>(
    s: &ZopfliBlockState<C>,
    costmodel: F,
    h: &mut ZopfliHash,
) -> (f64, Vec<u16>) {
    let in_data = s.data;
    let instart = s.blockstart;
    let inend = s.blockend;
    // Best cost to get here so far.
    let blocksize = inend - instart;
    let mut length_array = vec![0; blocksize + 1];
    if instart == inend {
        return (0.0, length_array);
    }
    let windowstart = instart.saturating_sub(ZOPFLI_WINDOW_SIZE);

    h.reset();
    let arr = &in_data[..inend];
    h.warmup(arr, windowstart, inend);
    for i in windowstart..instart {
        h.update(arr, i);
    }
    let mut costs: Vec<f32> = iter::repeat(f32::INFINITY).take(blocksize + 1).collect();
    costs[0] = 0.0; /* Because it's the start. */

    let mut i = instart;
    let mut leng;
    let mut longest_match;
    let mut sublen = vec![0; ZOPFLI_MAX_MATCH + 1];
    let mincost = get_cost_model_min_cost(&costmodel);
    while i < inend {
        let mut j = i - instart; // Index in the costs array and length_array.
        h.update(arr, i);

        // If we're in a long repetition of the same character and have more than
        // ZOPFLI_MAX_MATCH characters before and after our position.
        if h.same[i & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH as u16 * 2
            && i > instart + ZOPFLI_MAX_MATCH + 1
            && i + ZOPFLI_MAX_MATCH * 2 + 1 < inend
            && h.same[(i - ZOPFLI_MAX_MATCH) & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH as u16
        {
            let symbolcost = costmodel(ZOPFLI_MAX_MATCH, 1);
            // Set the length to reach each one to ZOPFLI_MAX_MATCH, and the cost to
            // the cost corresponding to that length. Doing this, we skip
            // ZOPFLI_MAX_MATCH values to avoid calling ZopfliFindLongestMatch.

            for _ in 0..ZOPFLI_MAX_MATCH {
                costs[j + ZOPFLI_MAX_MATCH] = costs[j] + symbolcost as f32;
                length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH as u16;
                i += 1;
                j += 1;
                h.update(arr, i);
            }
        }

        longest_match = find_longest_match(
            s,
            h,
            arr,
            i,
            inend,
            ZOPFLI_MAX_MATCH,
            &mut Some(&mut sublen),
        );
        leng = longest_match.length;

        // Literal.
        if i < inend {
            let new_cost = costmodel(arr[i] as usize, 0) + costs[j] as f64;
            debug_assert!(new_cost >= 0.0);
            if new_cost < costs[j + 1] as f64 {
                costs[j + 1] = new_cost as f32;
                length_array[j + 1] = 1;
            }
        }
        // Lengths.
        let kend = cmp::min(leng as usize, inend - i);
        let mincostaddcostj = mincost + costs[j] as f64;

        for (k, &sublength) in sublen.iter().enumerate().take(kend + 1).skip(3) {
            // Calling the cost model is expensive, avoid this if we are already at
            // the minimum possible cost that it can return.
            if costs[j + k] as f64 <= mincostaddcostj {
                continue;
            }

            let new_cost = costmodel(k, sublength) + costs[j] as f64;
            debug_assert!(new_cost >= 0.0);
            if new_cost < costs[j + k] as f64 {
                debug_assert!(k <= ZOPFLI_MAX_MATCH);
                costs[j + k] = new_cost as f32;
                length_array[j + k] = k as u16;
            }
        }
        i += 1;
    }

    debug_assert!(costs[blocksize] >= 0.0);
    (costs[blocksize] as f64, length_array)
}

/// Calculates the optimal path of lz77 lengths to use, from the calculated
/// `length_array`. The `length_array` must contain the optimal length to reach that
/// byte. The path will be filled with the lengths to use, so its data size will be
/// the amount of lz77 symbols.
fn trace(size: usize, length_array: &[u16]) -> Vec<u16> {
    let mut index = size;
    if size == 0 {
        return vec![];
    }
    let mut path = Vec::with_capacity(index);

    while index > 0 {
        let lai = length_array[index];
        let laiu = lai as usize;
        path.push(lai);
        debug_assert!(laiu <= index);
        debug_assert!(laiu <= ZOPFLI_MAX_MATCH);
        debug_assert_ne!(lai, 0);
        index -= laiu;
    }

    path
}

/// Does a single run for `lz77_optimal`. For good compression, repeated runs
/// with updated statistics should be performed.
/// `s`: the block state
/// `in_data`: the input data array
/// `instart`: where to start
/// `inend`: where to stop (not inclusive)
/// `length_array`: array of size `(inend - instart)` used to store lengths
/// `costmodel`: function to use as the cost model for this squeeze run
/// `store`: place to output the LZ77 data
/// returns the cost that was, according to the `costmodel`, needed to get to the end.
///     This is not the actual cost.
#[allow(clippy::too_many_arguments)] // Not feasible to refactor in a more readable way
fn lz77_optimal_run<F: Fn(usize, u16) -> f64, C: Cache>(
    s: &ZopfliBlockState<C>,
    costmodel: F,
    store: &mut Lz77Store,
    h: &mut ZopfliHash,
) {
    let instart = s.blockstart;
    let inend = s.blockend;
    let in_data = s.data;
    let (cost, length_array) = get_best_lengths(s, costmodel, h);
    let path = trace(inend - instart, &length_array);
    store.follow_path(in_data, instart, inend, path, s);
    debug_assert!(cost < f64::INFINITY);
}

/// Does the same as `lz77_optimal`, but optimized for the fixed tree of the
/// deflate standard.
/// The fixed tree never gives the best compression. But this gives the best
/// possible LZ77 encoding possible with the fixed tree.
/// This does not create or output any fixed tree, only LZ77 data optimized for
/// using with a fixed tree.
/// If `instart` is larger than `0`, it uses values before `instart` as starting
/// dictionary.
pub fn lz77_optimal_fixed<C: Cache>(s: &mut ZopfliBlockState<C>, store: &mut Lz77Store) {
    let mut h = ZopfliHash::new();
    lz77_optimal_run(s, get_cost_fixed, store, &mut h);
}

impl Genotype for SymbolTable {
    type Dna = usize;
}

impl RandomGenomeMutation for SymbolTable {
    type Dna = usize;

    fn mutate_genome<R>(
        mut genome: Self,
        mutation_rate: f64,
        min_value: &<Self as Genotype>::Dna,
        max_value: &<Self as Genotype>::Dna,
        rng: &mut R,
    ) -> Self
    where
        R: Rng + Sized,
    {
        genome.litlens.iter_mut().for_each(|litlen| {
            Self::mutate_plus_or_minus_1(mutation_rate, *min_value, *max_value, rng, litlen);
        });
        genome.dists.iter_mut().for_each(|dist| {
            Self::mutate_plus_or_minus_1(mutation_rate, *min_value, *max_value, rng, dist);
        });
        genome
    }
}

impl SymbolTable {
    fn mutate_plus_or_minus_1<R>(
        mutation_rate: f64,
        min_value: usize,
        max_value: usize,
        rng: &mut R,
        litlen: &mut usize,
    ) where
        R: Rng + Sized,
    {
        if rng.gen::<f64>() < mutation_rate {
            if *litlen <= min_value {
                *litlen += 1;
            } else if *litlen >= max_value {
                *litlen -= 1;
            } else if rng.gen_bool(0.5) {
                *litlen += 1;
            } else {
                *litlen -= 1;
            }
        }
    }
}

impl From<SymbolTable> for SymbolStats {
    fn from(value: SymbolTable) -> Self {
        let mut stats = SymbolStats {
            table: value,
            ..Default::default()
        };
        stats.calculate_entropy();
        stats
    }
}

impl Phenotype<SymbolTable> for SymbolStats {
    fn genes(&self) -> SymbolTable {
        self.table
    }

    fn derive(&self, new_genes: SymbolTable) -> Self {
        SymbolStats::from(new_genes)
    }
}

#[repr(transparent)]
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
struct FloatAsFitness(OrderedFloat<f64>);

impl Display for FloatAsFitness {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

impl From<f64> for FloatAsFitness {
    fn from(value: f64) -> Self {
        FloatAsFitness(value.into())
    }
}

impl From<FloatAsFitness> for f64 {
    fn from(value: FloatAsFitness) -> f64 {
        value.0.into()
    }
}

impl Fitness for FloatAsFitness {
    fn zero() -> Self {
        0.0.into()
    }

    fn abs_diff(&self, other: &Self) -> Self {
        f64::abs(self.0 .0 - other.0 .0).into()
    }
}

impl<'a, C> FitnessFunction<SymbolTable, FloatAsFitness> for &'a ZopfliBlockState<'a, C>
where
    C: Cache,
{
    fn fitness_of(&self, a: &SymbolTable) -> FloatAsFitness {
        let stats = SymbolStats::from(*a);
        let pool = &*LZ77_STORE_POOL;
        let mut currentstore = pool.pull();
        let mut h = ZopfliHash::new();
        lz77_optimal_run(
            self,
            |a, b| get_cost_stat(a, b, &stats),
            currentstore.deref_mut(),
            &mut h,
        );
        let cost = calculate_block_size(&currentstore, 0, currentstore.size(), BlockType::Dynamic);
        let mut best = self.best.lock().unwrap();
        let best = best.deref_mut();
        match best {
            None => *best = Some((currentstore.clone(), cost)),
            Some((ref mut best_store, ref mut best_cost)) => {
                if cost < *best_cost {
                    *best_cost = cost;
                    *best_store = currentstore.clone();
                }
            }
        }
        (-cost).into()
    }

    fn average(&self, a: &[FloatAsFitness]) -> FloatAsFitness {
        let mut total = 0.0;
        a.iter().for_each(|value| total += value.0 .0);
        (total / (a.len() as f64)).into()
    }

    fn highest_possible_fitness(&self) -> FloatAsFitness {
        (-1.0).into()
    }

    fn lowest_possible_fitness(&self) -> FloatAsFitness {
        f64::NEG_INFINITY.into()
    }
}

#[derive(Debug)]
pub struct SymbolTableBuilder {
    first_guess: SymbolTable,
    max_litlen_freq: usize,
    max_dist_freq: usize,
}

impl GenomeBuilder<SymbolTable> for SymbolTableBuilder {
    fn build_genome<R>(&self, index: usize, rng: &mut R) -> SymbolTable
    where
        R: Rng + Sized,
    {
        if index == 0 {
            return self.first_guess;
        }
        let mut table = SymbolTable::default();
        for litlen in table.litlens.iter_mut() {
            *litlen = rng.gen_range(0..=self.max_litlen_freq);
        }
        for dist in table.dists.iter_mut() {
            *dist = rng.gen_range(0..=self.max_dist_freq);
        }
        table
    }
}

#[derive(Debug)]
struct GenerationsWithoutImprovementLimiter {
    current_best: f64,
    generations_without_improvement: u64,
    max_generations_without_improvement: Option<u64>,
    generations: u64,
    max_generations: Option<u64>,
}

impl GenerationsWithoutImprovementLimiter {
    fn new(max_generations_without_improvement: Option<u64>, max_generations: Option<u64>) -> Self {
        GenerationsWithoutImprovementLimiter {
            current_best: f64::NEG_INFINITY,
            generations_without_improvement: 0,
            max_generations_without_improvement,
            generations: 0,
            max_generations,
        }
    }
}

impl<A, G> Termination<A> for GenerationsWithoutImprovementLimiter
where
    A: Algorithm<Output = ga::State<G, FloatAsFitness>>,
    G: Genotype,
{
    fn evaluate(&mut self, state: &State<A>) -> StopFlag {
        self.generations += 1;
        if self
            .max_generations
            .is_some_and(|max_gens| self.generations >= max_gens)
        {
            StopFlag::StopNow("Maximum generations reached".into())
        } else {
            let fitness = state.result.best_solution.solution.fitness.0 .0;
            if fitness >= -1.0 {
                StopFlag::StopNow("Already only 1 bit".into())
            } else if fitness > self.current_best {
                self.current_best = fitness;
                self.generations_without_improvement = 0;
                StopFlag::Continue
            } else {
                self.generations_without_improvement += 1;
                if self
                    .max_generations_without_improvement
                    .is_some_and(|max_gens| max_gens <= self.generations_without_improvement)
                {
                    StopFlag::StopNow("Generations-without-improvement limit reached".into())
                } else {
                    StopFlag::Continue
                }
            }
        }
    }
}

#[derive(Copy, Clone, Default, Debug)]
struct SymbolTableCrossBreeder {}

impl GeneticOperator for SymbolTableCrossBreeder {
    fn name() -> String {
        "SymbolTableCrossBreeder".to_string()
    }
}

fn generate_child_chromosomes<T, const N: usize, R>(
    parent0: [T; N],
    parent1: [T; N],
    rng: &mut R,
) -> [[T; N]; 4]
where
    T: Copy,
    R: Rng + Sized,
{
    let cut_point = rng.gen_range(1..N as i32) as usize;
    let mut hybrid_0 = parent0;
    hybrid_0[cut_point..].copy_from_slice(&parent1[cut_point..]);
    let mut hybrid_1 = parent1;
    hybrid_1[cut_point..].copy_from_slice(&parent0[cut_point..]);
    [parent0, parent1, hybrid_0, hybrid_1]
}

impl CrossoverOp<SymbolTable> for SymbolTableCrossBreeder {
    fn crossover<R>(&self, parents: Parents<SymbolTable>, rng: &mut R) -> Children<SymbolTable>
    where
        R: Rng + Sized,
    {
        let litlens = generate_child_chromosomes(parents[0].litlens, parents[1].litlens, rng);
        let dists = generate_child_chromosomes(parents[0].dists, parents[1].dists, rng);
        litlens
            .iter()
            .flat_map(|litlens| {
                dists.iter().map(|dists| SymbolTable {
                    litlens: *litlens,
                    dists: *dists,
                })
            })
            .collect()
    }
}

const POPULATION_SIZE: usize = 10;

/// Calculates lit/len and dist pairs for given data.
/// If `instart` is larger than 0, it uses values before `instart` as starting
/// dictionary.
pub fn lz77_optimal<C: Cache>(
    s: &mut ZopfliBlockState<C>,
    in_data: &[u8],
    max_iterations: Option<u64>,
    max_iterations_without_improvement: Option<u64>,
) -> Lz77Store {
    let instart = s.blockstart;
    let inend = s.blockend;
    /* Dist to get to here with smallest cost. */
    let mut outputstore = Lz77Store::new();

    /* Initial run. */
    outputstore.greedy(s, in_data, instart, inend);
    let mut greedy_stats = SymbolStats::default();
    greedy_stats.get_statistics(&outputstore);
    let max_litlen_freq = greedy_stats.table.litlens.iter().max().unwrap() + 1;
    let max_dist_freq = greedy_stats.table.dists.iter().max().unwrap() + 1;
    let genome_builder = SymbolTableBuilder {
        max_dist_freq,
        max_litlen_freq,
        first_guess: greedy_stats.table,
    };
    let initial_population = build_population()
        .with_genome_builder(genome_builder)
        .of_size(POPULATION_SIZE)
        .uniform_at_random();
    let algorithm = genetic_algorithm()
        .with_evaluation(&*s)
        .with_selection(MaximizeSelector::new(0.85, 12))
        .with_crossover(SymbolTableCrossBreeder::default())
        .with_mutation(RandomValueMutator::new(
            0.2,
            0,
            max_litlen_freq.max(max_dist_freq),
        ))
        .with_reinsertion(ElitistReinserter::new(&*s, false, 0.85))
        .with_initial_population(initial_population)
        .build();
    let mut genetic_algorithm_sim = simulate(algorithm)
        .until(GenerationsWithoutImprovementLimiter::new(
            max_iterations_without_improvement,
            max_iterations,
        ))
        .build();
    loop {
        match genetic_algorithm_sim.step() {
            Ok(SimResult::Intermediate(step)) => {
                let evaluated_population = step.result.evaluated_population;
                let best_solution = step.result.best_solution;
                debug!(
                    "step: generation: {}, average_fitness: {}, \
                     best fitness: {}, duration: {}, processing_time: {}",
                    step.iteration,
                    evaluated_population.average_fitness(),
                    best_solution.solution.fitness,
                    step.duration,
                    step.processing_time,
                );
            }
            Ok(SimResult::Final(step, processing_time, duration, stop_reason)) => {
                debug!(
                    "final result: generation: {},\
                         best fitness: {}, duration: {}, processing_time: {}, stop_reason: {}",
                    step.iteration,
                    step.result.best_solution.solution.fitness,
                    duration,
                    processing_time,
                    stop_reason
                );
                let best = s.best.lock().unwrap();
                return best.clone().unwrap().0;
            }
            Err(e) => panic!("{}", e),
        }
    }
}
