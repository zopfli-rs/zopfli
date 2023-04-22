use alloc::vec::Vec;

use crate::util::{ZOPFLI_MIN_MATCH, ZOPFLI_WINDOW_MASK, ZOPFLI_WINDOW_SIZE};

const HASH_SHIFT: i32 = 5;
const HASH_MASK: u16 = 32767;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Which {
    Hash1,
    Hash2,
}

pub struct SmallerHashThing {
    prev: u16,            /* Index to index of prev. occurrence of same hash. */
    hashval: Option<u16>, /* Index to hash value at this index. */
}

pub struct HashThing {
    head: Vec<i32>, /* Hash value to index of its most recent occurrence. */
    prev_and_hashval: Vec<SmallerHashThing>,
    val: u16, /* Current hash value. */
}

impl HashThing {
    fn new() -> HashThing {
        HashThing {
            head: vec![-1; 65536],
            prev_and_hashval: (0..ZOPFLI_WINDOW_SIZE)
                .map(|p| SmallerHashThing {
                    prev: p as u16,
                    hashval: None,
                })
                .collect(),
            val: 0,
        }
    }

    fn reset(&mut self) {
        self.val = 0;

        self.head.fill(-1);

        let mut p = 0;
        self.prev_and_hashval.fill_with(|| {
            let thing = SmallerHashThing {
                prev: p,
                hashval: None,
            };
            p += 1;
            thing
        });
    }

    fn update(&mut self, hpos: usize) {
        let hashval = self.val;
        let index = self.val as usize;
        let head_index = self.head[index];
        let prev = if head_index != -1
            && self.prev_and_hashval[head_index as usize]
                .hashval
                .map_or(false, |hv| hv == self.val)
        {
            head_index as u16
        } else {
            hpos as u16
        };

        self.prev_and_hashval[hpos] = SmallerHashThing {
            prev,
            hashval: Some(hashval),
        };
        self.head[index] = hpos as i32;
    }
}

pub struct ZopfliHash {
    hash1: HashThing,
    hash2: HashThing,
    pub same: [u16; ZOPFLI_WINDOW_SIZE], /* Amount of repetitions of same byte after this .*/
}

impl ZopfliHash {
    pub fn new() -> ZopfliHash {
        ZopfliHash {
            hash1: HashThing::new(),
            hash2: HashThing::new(),
            same: [0; ZOPFLI_WINDOW_SIZE],
        }
    }

    pub fn reset(&mut self) {
        self.hash1.reset();
        self.hash2.reset();
        self.same = [0; ZOPFLI_WINDOW_SIZE];
    }

    pub fn warmup(&mut self, arr: &[u8], pos: usize, end: usize) {
        let c = arr[pos];
        self.update_val(c);

        if pos + 1 < end {
            let c = arr[pos + 1];
            self.update_val(c);
        }
    }

    /// Update the sliding hash value with the given byte. All calls to this function
    /// must be made on consecutive input characters. Since the hash value exists out
    /// of multiple input bytes, a few warmups with this function are needed initially.
    fn update_val(&mut self, c: u8) {
        self.hash1.val = ((self.hash1.val << HASH_SHIFT) ^ c as u16) & HASH_MASK;
    }

    pub fn update(&mut self, array: &[u8], pos: usize) {
        let hash_value = array.get(pos + ZOPFLI_MIN_MATCH - 1).cloned().unwrap_or(0);
        self.update_val(hash_value);

        let hpos = pos & ZOPFLI_WINDOW_MASK;

        self.hash1.update(hpos);

        // Update "same".
        let mut amount = 0;
        let same_index = pos.wrapping_sub(1) & ZOPFLI_WINDOW_MASK;
        let same = self.same[same_index];
        if same > 1 {
            amount = same - 1;
        }

        self.same[hpos] = amount;

        self.hash2.val = (amount.wrapping_sub(ZOPFLI_MIN_MATCH as u16) & 255) ^ self.hash1.val;

        self.hash2.update(hpos);
    }

    pub fn head_at(&self, index: usize, which: Which) -> i32 {
        match which {
            Which::Hash1 => self.hash1.head[index],
            Which::Hash2 => self.hash2.head[index],
        }
    }

    pub fn prev_at(&self, index: usize, which: Which) -> u16 {
        match which {
            Which::Hash1 => self.hash1.prev_and_hashval[index].prev,
            Which::Hash2 => self.hash2.prev_and_hashval[index].prev,
        }
    }

    pub fn hash_val_at(&self, index: usize, which: Which) -> i32 {
        let hashval = match which {
            Which::Hash1 => self.hash1.prev_and_hashval[index].hashval,
            Which::Hash2 => self.hash2.prev_and_hashval[index].hashval,
        };
        hashval.map_or(-1, |hv| hv as i32)
    }

    pub fn val(&self, which: Which) -> u16 {
        match which {
            Which::Hash1 => self.hash1.val,
            Which::Hash2 => self.hash2.val,
        }
    }
}
