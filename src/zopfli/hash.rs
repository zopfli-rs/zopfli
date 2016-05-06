use std::slice;

use libc::{c_int, c_ushort, c_uchar, size_t};

use symbols::{ZOPFLI_WINDOW_MASK, ZOPFLI_MIN_MATCH};

const HASH_SHIFT: c_int = 5;
const HASH_MASK: c_int = 32767;

#[repr(C)]
pub struct ZopfliHash {
    head: Vec<c_int>,  /* Hash value to index of its most recent occurrence. */
    prev: Vec<c_ushort>,  /* Index to index of prev. occurrence of same hash. */
    hashval: Vec<c_int>,  /* Index to hash value at this index. */
    val: c_int,  /* Current hash value. */

    /* Fields with similar purpose as the above hash, but for the second hash with
    a value that is calculated differently.  */
    head2: Vec<c_int>,  /* Hash value to index of its most recent occurrence. */
    prev2: Vec<c_ushort>,  /* Index to index of prev. occurrence of same hash. */
    hashval2: Vec<c_int>,  /* Index to hash value at this index. */
    val2: c_int,  /* Current hash value. */

    same: Vec<c_ushort>,  /* Amount of repetitions of same byte after this .*/
}

impl ZopfliHash {
    pub fn new(window_size: size_t) -> ZopfliHash {
        ZopfliHash {
           head: vec![-1; 65536],
           prev: (0..window_size as c_ushort).collect::<Vec<_>>(),
           hashval: vec![-1; window_size],
           val: 0,

           /* Fields with similar purpose as the above hash, but for the second hash with
           a value that is calculated differently.  */
           head2: vec![-1; 65536],
           prev2: (0..window_size as c_ushort).collect::<Vec<_>>(),
           hashval2: vec![-1; window_size],
           val2: 0,

           same: vec![0; window_size],
       }
    }

    pub fn warmup(&mut self, arr: &[c_uchar], pos: size_t, end: size_t) {
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
    pub fn update_val(&mut self, c: c_uchar) {
        self.val = ((self.val << HASH_SHIFT) ^ c as c_int) & HASH_MASK;
    }

    pub fn update(&mut self, array: &[c_uchar], pos: size_t) {
        let hpos = (pos & ZOPFLI_WINDOW_MASK) as usize;
        let mut amount: c_int = 0;

        let hash_value = if pos + ZOPFLI_MIN_MATCH as size_t <= array.len() {
            array[(pos + ZOPFLI_MIN_MATCH as size_t - 1) as usize]
        } else {
            0
        };
        self.update_val(hash_value);

        self.hashval[hpos] = self.val;

        let index = self.val as usize;

        if self.head[index] != -1 && self.hashval[self.head[index] as usize] == self.val {
            self.prev[hpos] = self.head[index] as c_ushort;
        } else {
            self.prev[hpos] = hpos as c_ushort;
        }

        self.head[index] = hpos as c_int;

        // Update "same".
        if self.same[(pos.wrapping_sub(1) & ZOPFLI_WINDOW_MASK) as usize] > 1 {
            amount = self.same[(pos.wrapping_sub(1) & ZOPFLI_WINDOW_MASK) as usize] as c_int - 1;
        }

        while pos + amount as size_t + 1 < array.len() && array[pos as usize] == array[(pos + amount as size_t + 1) as usize] && amount < -1 {
            amount += 1;
        }
        self.same[hpos] = amount as c_ushort;

        self.val2 = ((self.same[hpos].wrapping_sub(ZOPFLI_MIN_MATCH as c_ushort) & 255) ^ self.val as c_ushort) as c_int;
        self.hashval2[hpos] = self.val2;

        let index2 = self.val2 as usize;
        if self.head2[index2] != -1 as i32 && self.hashval2[self.head2[index2] as usize] == self.val2 {
            self.prev2[hpos] = self.head2[index2] as c_ushort;
        } else {
            self.prev2[hpos] = hpos as c_ushort;
        }
        self.head2[index2] = hpos as c_int;
    }

    pub fn head_at(&self, index: size_t, which_hash: size_t) -> c_int {
        if which_hash == 1 {
            self.head[index]
        } else {
            self.head2[index]
        }
    }

    pub fn prev_at(&self, index: size_t, which_hash: size_t) -> c_ushort {
        if which_hash == 1 {
            self.prev[index]
        } else {
            self.prev2[index]
        }
    }

    pub fn hash_val_at(&self, index: size_t, which_hash: size_t) -> c_int {
        if which_hash == 1 {
            self.hashval[index]
        } else {
            self.hashval2[index]
        }
    }

    pub fn val(&self, which_hash: size_t) -> c_int {
        if which_hash == 1 {
            self.val
        } else {
            self.val2
        }
    }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliWarmupHash(array: *const c_uchar, pos: size_t, end: size_t, h_ptr: *mut ZopfliHash) {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    let arr = unsafe { slice::from_raw_parts(array, end) };
    h.warmup(arr, pos, end);
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliUpdateHash(array: *const c_uchar, pos: size_t, end: size_t, h_ptr: *mut ZopfliHash) {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    let arr = unsafe { slice::from_raw_parts(array, end) };
    h.update(arr, pos);
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliHashSameAt(h_ptr: *mut ZopfliHash, index: size_t) -> c_ushort {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };
    h.same[index]
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliCleanHash(ptr: *mut ZopfliHash) {
    if ptr.is_null() { return }
    unsafe { Box::from_raw(ptr); }
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliInitHash(window_size: size_t) -> *mut ZopfliHash {
    Box::into_raw(Box::new(ZopfliHash::new(window_size)))
}

#[no_mangle]
#[allow(non_snake_case)]
pub extern fn ZopfliResetHash(window_size: size_t, h_ptr: *mut ZopfliHash) {
    let h = unsafe {
        assert!(!h_ptr.is_null());
        &mut *h_ptr
    };

    h.val = 0;
    h.head = vec![-1; 65536];
    h.prev = (0..window_size as c_ushort).collect::<Vec<_>>();
    h.hashval = vec![-1; window_size];

    h.same = vec![0; window_size];

    h.val2 = 0;
    h.head2 = vec![-1; 65536];
    h.prev2 = (0..window_size as c_ushort).collect::<Vec<_>>();
    h.hashval2 = vec![-1; window_size];
}
