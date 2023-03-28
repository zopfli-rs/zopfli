#[cfg(not(feature = "nightly"))]
compile_error!("zopfli depends on the unstable feature 'core_intrinsics' to function in no_std environments. Please enable its 'nightly' feature");

/// Provides math operations for doubles on `no_std` targets that are not available on `core`.
pub trait F64MathExt {
    /// Computes the absolute value of `self`.
    #[must_use = "method returns a new number and does not mutate the original value"]
    fn abs(self) -> Self;

    /// Returns the natural logarithm of the number.
    ///
    /// This method is allowed to return different values for 0 and negative numbers
    /// than its counterpart on `std` on some exotic platforms.
    #[must_use = "method returns a new number and does not mutate the original value"]
    fn ln(self) -> Self;
}

impl F64MathExt for f64 {
    #[inline]
    fn abs(self) -> Self {
        // SAFETY: trivially safe
        #[cfg(feature = "nightly")]
        unsafe {
            core::intrinsics::fabsf64(self)
        }
        #[cfg(not(feature = "nightly"))]
        // Stub return value to show a single compiler error (see first line)
        0.0
    }

    #[inline]
    fn ln(self) -> Self {
        // SAFETY: trivially safe. See std's implementation for more details about
        //         the different (but safe) behavior on input <= 0
        #[cfg(feature = "nightly")]
        unsafe {
            core::intrinsics::logf64(self)
        }
        #[cfg(not(feature = "nightly"))]
        // Stub return value to show a single compiler error (see first line)
        0.0
    }
}
