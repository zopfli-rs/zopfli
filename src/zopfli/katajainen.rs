use std::cmp::Ordering;
use std::{mem};

use libc::{size_t, c_int};

// Bounded package merge algorithm, based on the paper
// "A Fast and Space-Economical Algorithm for Length-Limited Coding
// Jyrki Katajainen, Alistair Moffat, Andrew Turpin".

#[derive(Debug)]
struct Node {
    weight: size_t,
    leaf_counts: Vec<c_int>,
}

impl Node {
    pub fn new(weight: size_t, initial_count: c_int, capacity: usize) -> Node {
        let mut n = Node {
            weight: weight,
            leaf_counts: Vec::with_capacity(capacity),
        };
        n.leaf_counts.push(initial_count);
        n
    }
}

#[derive(Debug)]
struct Leaf {
    pub weight: size_t,
    pub index: size_t,
}
impl PartialEq for Leaf {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
    }
}
impl Eq for Leaf { }
impl Ord for Leaf {
    fn cmp(&self, other: &Self) -> Ordering {
        self.weight.cmp(&other.weight)
    }
}
impl PartialOrd for Leaf {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
struct List {
    lookahead1: Node,
    lookahead2: Node,
    next_leaf_index: size_t,
}

pub fn length_limited_code_lengths(frequencies: &[size_t], maxbits: c_int) -> Vec<size_t> {
    let mut leaves = vec![];

    // Count used symbols and place them in the leaves.
    for (i, &freq) in frequencies.iter().enumerate() {
        if freq != 0 {
            leaves.push(Leaf { weight: freq, index: i });
        }
    }

    // Short circuit some special cases
    if leaves.is_empty() {
        // There are no non-zero frequencies.
        return vec![0; frequencies.len()];
    }
    if leaves.len() == 1 {
        let mut result = vec![0; frequencies.len()];
        result[leaves[0].index] = 1;
        return result;
    }
    if leaves.len() == 2 {
        let mut result = vec![0; frequencies.len()];
        result[leaves[0].index] = 1;
        result[leaves[1].index] = 1;
        return result;
    }

    // Sort the leaves from least frequent to most frequent.
    // Add index into the same variable for stable sorting.
    for leaf in leaves.iter_mut() {
        leaf.weight = (leaf.weight << 9) | leaf.index;
    }
    leaves.sort();
    for leaf in leaves.iter_mut() {
        leaf.weight >>= 9;
    }

    let max_num_leaves = 2 * leaves.len() - 2;

    let mut lists = Vec::with_capacity(maxbits as usize);
    for _ in 0..maxbits {
        lists.push(List {
            lookahead1: Node::new(leaves[0].weight, 1, max_num_leaves),
            lookahead2: Node::new(leaves[1].weight, 2, max_num_leaves),
            next_leaf_index: 2,
        });
    }

    // In the last list, 2 * numsymbols - 2 active chains need to be created. Two
    // are already created in the initialization. Each boundary_pm run creates one.
    let num_boundary_pm_runs = max_num_leaves - 2;
    for _ in 0..num_boundary_pm_runs {
        boundary_pm_toplevel(&mut lists[..], &leaves);
    }

    let n = frequencies.len();
    let mut result = vec![0; n];

    let mut a = lists.pop().unwrap().lookahead2.leaf_counts.into_iter().rev().peekable();

    let mut bitlength_value = 1;
    while let Some(leaf_count) = a.next() {
        let next_count = *a.peek().unwrap_or(&0);
        for i in next_count..leaf_count {
            result[leaves[i as usize].index as usize] = bitlength_value;
        }
        bitlength_value += 1;
    }

    result
}

fn lowest_list(lists: &mut [List], leaves: &Vec<Leaf>) {
    // We're in the lowest list, just add another leaf to the lookaheads
    // There will always be more leaves to be added on level 0 so this is safe.
    let mut current_list = lists.get_mut(0).unwrap();
    let ref next_leaf = leaves[current_list.next_leaf_index];
    current_list.lookahead2.weight = next_leaf.weight;

    current_list.lookahead2.leaf_counts[0] = current_list.lookahead1.leaf_counts.last().unwrap() + 1;
    current_list.next_leaf_index += 1;
}

fn next_leaf(lists: &mut [List], leaves: &Vec<Leaf>, current_list_index: usize) {
    let mut current_list = lists.get_mut(current_list_index).unwrap();

    // The next leaf goes next; counting itself makes the leaf_count increase by one.
    current_list.lookahead2.weight = leaves[current_list.next_leaf_index].weight;
    current_list.lookahead2.leaf_counts.clear();
    current_list.lookahead2.leaf_counts.extend(current_list.lookahead1.leaf_counts.iter());
    let last_index = current_list.lookahead2.leaf_counts.len() - 1;
    current_list.lookahead2.leaf_counts[last_index] += 1;
    current_list.next_leaf_index += 1;
}

fn next_tree(weight_sum: size_t, lists: &mut [List], leaves: &Vec<Leaf>, current_list_index: usize) {
    let num_leaf_counts = lists[current_list_index - 1].lookahead2.leaf_counts.len();
    let previous_list_leaf_counts = lists[current_list_index - 1].lookahead2.leaf_counts.as_ptr();
    {
        let ref mut current_list = lists[current_list_index];

        // Make a tree from the lookaheads from the previous list; that goes next.
        // This is not a leaf node, so the leaf count stays the same.
        current_list.lookahead2.weight = weight_sum;
        current_list.lookahead2.leaf_counts.clear();
        for i in 0..num_leaf_counts {
            current_list.lookahead2.leaf_counts.push(unsafe {
                *previous_list_leaf_counts.offset(i as isize)
            });
        }
        current_list.lookahead2.leaf_counts.push(*current_list.lookahead1.leaf_counts.last().unwrap());
    }

    // The previous list needs two new lookahead nodes.
    boundary_pm(lists, leaves, current_list_index - 1);
    boundary_pm(lists, leaves, current_list_index - 1);
}

fn boundary_pm_toplevel(lists: &mut [List], leaves: &Vec<Leaf>) {
    let last_index = lists.len() - 1;
    boundary_pm(lists, leaves, last_index);
}

fn boundary_pm(lists: &mut [List], leaves: &Vec<Leaf>, current_list_index: usize) {
    let next_leaf_index = lists[current_list_index].next_leaf_index;

    if current_list_index == 0 && next_leaf_index == leaves.len() {
        // We've added all the leaves to the lowest list, so we're done here
        return;
    }

    mem::swap(&mut lists[current_list_index].lookahead1, &mut lists[current_list_index].lookahead2);

    if current_list_index == 0 {
        lowest_list(lists, leaves);
    } else {
        // We're at a list other than the lowest list.
        let weight_sum = {
            let previous_list = lists.get(current_list_index - 1).unwrap();
            previous_list.lookahead1.weight + previous_list.lookahead2.weight
        };

        if next_leaf_index < leaves.len() && weight_sum > leaves[next_leaf_index].weight {
            next_leaf(lists, leaves, current_list_index);
        } else {
            next_tree(weight_sum, lists, leaves, current_list_index);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_from_paper_3() {
        let input = [1, 1, 5, 7, 10, 14];
        let output = length_limited_code_lengths(&input, 3);
        let answer = vec![3, 3, 3, 3, 2, 2];
        assert_eq!(output, answer);
    }

    #[test]
    fn test_from_paper_4() {
        let input = [1, 1, 5, 7, 10, 14];
        let output = length_limited_code_lengths(&input, 4);
        let answer = vec![4, 4, 3, 2, 2, 2];
        assert_eq!(output, answer);
    }

    #[test]
    fn maxbits_7() {
        let input = [252, 0, 1, 6, 9, 10, 6, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let output = length_limited_code_lengths(&input, 7);
        let answer = vec![1, 0, 6, 4, 3, 3, 3, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        assert_eq!(output, answer);
    }

    #[test]
    fn maxbits_15() {
        let input = [0, 0, 0, 0, 0, 0, 18, 0, 6, 0, 12, 2, 14, 9, 27, 15, 23, 15, 17, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let output = length_limited_code_lengths(&input, 15);
        let answer = vec! [0, 0, 0, 0, 0, 0, 3, 0, 5, 0, 4, 6, 4, 4, 3, 4, 3, 3, 3, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        assert_eq!(output, answer);
    }

    #[test]
    fn no_frequencies() {
        let input = [0, 0, 0, 0, 0];
        let output = length_limited_code_lengths(&input, 7);
        let answer = vec![0, 0, 0, 0, 0];
        assert_eq!(output, answer);
    }

    #[test]
    fn only_one_frequency() {
        let input = [0, 10, 0];
        let output = length_limited_code_lengths(&input, 7);
        let answer = vec![0, 1, 0];
        assert_eq!(output, answer);
    }

    #[test]
    fn only_two_frequencies() {
        let input = [0, 0, 0, 0, 252, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let output = length_limited_code_lengths(&input, 7);
        let answer = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        assert_eq!(output, answer);
    }
}
