use nalgebra::DMatrix;

use crate::constants::{M_E_KEV, PI_OVER_FOUR, ONE_OVER_SQRT_2};
use crate::cdbs::cdbspectrum::Area;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AntiAliasingScheme {
    /// A pixel contributes fully or not at all, depending on whether its center is inside the bin or not
    None,
    /// Rotate the matrix first so that the bins are axis-aligned, then compute overlaps
    Rotate,
    /// A pixel contributes with its fractional coverage which is computed for each bin edge separately, the results getting multiplied together
    EdgewiseScanline,
    /// A pixel contributes with its fractional coverage which is computed for all bin edges simultaneously
    EdgeWalk,
}
