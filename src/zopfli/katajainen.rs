use std::cmp::Ordering;
use std::slice;

use libc::{size_t, c_int, c_uint};

// Bounded package merge algorithm, based on the paper
// "A Fast and Space-Economical Algorithm for Length-Limited Coding
// Jyrki Katajainen, Alistair Moffat, Andrew Turpin".

#[derive(Debug)]
struct Node {
    weight: size_t,
    leaf_counts: Vec<c_int>,
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

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliLengthLimitedCodeLengths(frequencies: *const size_t, n: c_int, maxbits: c_int, bitlengths: *mut c_uint) -> c_int {
    let freqs = unsafe { slice::from_raw_parts(frequencies, n as usize) };
    let result = length_limited_code_lengths(freqs, maxbits);

    for (i, res) in result.into_iter().enumerate() {
        unsafe {
            *bitlengths.offset(i as isize) = res as c_uint;
        }
    }

    return 0;
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

    let mut lists = Vec::with_capacity(maxbits as usize);
    for _ in 0..maxbits {
        lists.push(List {
            lookahead1: Node { weight: leaves[0].weight, leaf_counts: vec![1] },
            lookahead2: Node { weight: leaves[1].weight, leaf_counts: vec![2] },
            next_leaf_index: 2,
        });
    }

    // In the last list, 2 * numsymbols - 2 active chains need to be created. Two
    // are already created in the initialization. Each boundary_pm run creates one.
    let num_boundary_pm_runs = 2 * leaves.len() - 4;
    for _ in 0..num_boundary_pm_runs {
        lists = boundary_pm(lists, &leaves);
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

fn boundary_pm(mut lists: Vec<List>, leaves: &Vec<Leaf>) -> Vec<List> {
    let mut current_list = lists.pop().unwrap();
    if lists.is_empty() && current_list.next_leaf_index == leaves.len() {
        // We've added all the leaves to the lowest list, so we're done here
        lists.push(current_list);
        return lists;
    }

    current_list.lookahead1 = current_list.lookahead2;

    if lists.is_empty() {
        // We're in the lowest list, just add another leaf to the lookaheads
        // There will always be more leaves to be added on level 0 so this is safe.
        let ref next_leaf = leaves[current_list.next_leaf_index];
        current_list.lookahead2 = Node {
            weight: next_leaf.weight,
            leaf_counts: vec![current_list.lookahead1.leaf_counts.last().unwrap() + 1],
        };
        current_list.next_leaf_index += 1;
        lists.push(current_list);
    } else {
        // We're at a list other than the lowest list.
        let previous_list = lists.pop().unwrap();
        let weight_sum = previous_list.lookahead1.weight + previous_list.lookahead2.weight;

        if current_list.next_leaf_index < leaves.len() && weight_sum > leaves[current_list.next_leaf_index].weight {
            // The next leaf goes next; counting itself makes the leaf_count increase by one.
            let mut last_leaf_counts = current_list.lookahead1.leaf_counts.clone();
            let mut last_count = last_leaf_counts.pop().unwrap();
            last_count += 1;
            last_leaf_counts.push(last_count);
            current_list.lookahead2 = Node {
                weight: leaves[current_list.next_leaf_index].weight,
                leaf_counts: last_leaf_counts,
            };
            current_list.next_leaf_index += 1;
            lists.push(previous_list);
            lists.push(current_list);
        } else {
            // Make a tree from the lookaheads from the previous list; that goes next.
            // This is not a leaf node, so the leaf count stays the same.
            let mut last_leaf_counts = previous_list.lookahead2.leaf_counts.clone();
            last_leaf_counts.push(*current_list.lookahead1.leaf_counts.last().unwrap());
            current_list.lookahead2 = Node {
                weight: weight_sum,
                leaf_counts: last_leaf_counts,
            };
            // The previous list needs two new lookahead nodes.
            lists.push(previous_list);
            lists = boundary_pm(lists, leaves);
            lists = boundary_pm(lists, leaves);
            lists.push(current_list);
        }
    }

    lists
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
