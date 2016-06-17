use util::{ZOPFLI_WINDOW_MASK, ZOPFLI_MIN_MATCH};

const HASH_SHIFT: i32 = 5;
const HASH_MASK: i32 = 32767;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Which {
    Hash1,
    Hash2,
}

pub struct ZopfliHash {
    head: Vec<i32>,  /* Hash value to index of its most recent occurrence. */
    prev: Vec<u16>,  /* Index to index of prev. occurrence of same hash. */
    hashval: Vec<i32>,  /* Index to hash value at this index. */
    val: i32,  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    head2: Vec<i32>,  /* Hash value to index of its most recent occurrence. */
    prev2: Vec<u16>,  /* Index to index of prev. occurrence of same hash. */
    hashval2: Vec<i32>,  /* Index to hash value at this index. */
    val2: i32,  /* Current hash value. */

    pub same: Vec<u16>,  /* Amount of repetitions of same byte after this .*/
}

impl ZopfliHash {
    pub fn new(window_size: usize) -> ZopfliHash {
        ZopfliHash {
           head: vec![-1; 65536],
           prev: (0..window_size as u16).collect(),
           hashval: vec![-1; window_size],
           val: 0,

           /* Fields with similar purpose as the above hash, but for the second hash with
           a value that is calculated differently.  */
           head2: vec![-1; 65536],
           prev2: (0..window_size as u16).collect(),
           hashval2: vec![-1; window_size],
           val2: 0,

           same: vec![0; window_size],
       }
    }

    pub fn reset(&mut self, window_size: usize) {
        self.val = 0;
        self.head = vec![-1; 65536];
        self.prev = (0..window_size as u16).collect();
        self.hashval = vec![-1; window_size];

        self.same = vec![0; window_size];
        self.val2 = 0;
        self.head2 = vec![-1; 65536];
        self.prev2 = (0..window_size as u16).collect();
        self.hashval2 = vec![-1; window_size];
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
        self.val = ((self.val << HASH_SHIFT) ^ c as i32) & HASH_MASK;
    }

    pub fn update(&mut self, array: &[u8], pos: usize) {
        let hash_value = array.get(pos + ZOPFLI_MIN_MATCH - 1).cloned().unwrap_or(0);
        self.update_val(hash_value);

        let hpos = pos & ZOPFLI_WINDOW_MASK;
        self.hashval[hpos] = self.val;

        let index = self.val as usize;
        if self.head[index] != -1 && self.hashval[self.head[index] as usize] == self.val {
            self.prev[hpos] = self.head[index] as u16;
        } else {
            self.prev[hpos] = hpos as u16;
        }
        self.head[index] = hpos as i32;

        // Update "same".
        let mut amount = 0;
        let same_index = pos.wrapping_sub(1) & ZOPFLI_WINDOW_MASK;
        if self.same[same_index] > 1 {
            amount = self.same[same_index] as i32 - 1;
        }

        let mut another_index = pos + amount as usize + 1;
        while another_index < array.len() && array[pos] == array[another_index] && amount < -1 {
            amount += 1;
            another_index += 1;
        }
        self.same[hpos] = amount as u16;

        self.val2 = ((self.same[hpos].wrapping_sub(ZOPFLI_MIN_MATCH as u16) & 255) ^ self.val as u16) as i32;
        self.hashval2[hpos] = self.val2;

        let index2 = self.val2 as usize;
        if self.head2[index2] != -1 && self.hashval2[self.head2[index2] as usize] == self.val2 {
            self.prev2[hpos] = self.head2[index2] as u16;
        } else {
            self.prev2[hpos] = hpos as u16;
        }
        self.head2[index2] = hpos as i32;
    }

    pub fn head_at(&self, index: usize, which: Which) -> i32 {
        match which {
            Which::Hash1 => self.head[index],
            Which::Hash2 => self.head2[index],
        }
    }

    pub fn prev_at(&self, index: usize, which: Which) -> u16 {
        match which {
            Which::Hash1 => self.prev[index],
            Which::Hash2 => self.prev2[index],
        }
    }

    pub fn hash_val_at(&self, index: usize, which: Which) -> i32 {
        match which {
            Which::Hash1 => self.hashval[index],
            Which::Hash2 => self.hashval2[index],
        }
    }

    pub fn val(&self, which: Which) -> i32 {
        match which {
            Which::Hash1 => self.val,
            Which::Hash2 => self.val2,
        }
    }
}
