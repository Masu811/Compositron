use std::sync::Arc;

use num_traits::Unsigned;
use thiserror::Error;
use nalgebra::{DMatrix, DVector};

#[derive(Debug, Error)]
pub enum FitError {
    #[error("Data contains only NaN values")]
    AllNANValues,

    #[error("ModelBuilderError")]
    ModelBuilderError {
        #[from]
        inner: varpro::model::builder::error::ModelBuildError,
    },

    #[error("SeparableProblemBuilderError")]
    SeparableProblemBuilderError {
        #[from]
        inner: varpro::problem::SeparableProblemBuilderError,
    },

    #[error("Something went wrong during least squares fit")]
    RuntimeError,
}

#[derive(Debug)]
pub struct DetectorCalibration {
    pub offset: f64,
    pub scale: f64,
}

#[derive(Debug)]
pub struct CoincDetectorCalibration {
    pub first_det: Arc<DetectorCalibration>,
    pub second_det: Arc<DetectorCalibration>,
}

pub enum Spectrum {
    U8(DVector<u8>),
    U16(DVector<u16>),
    U32(DVector<u32>),
    U64(DVector<u64>),
}

#[macro_export]
macro_rules! spectrum_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            $crate::core::utils::Spectrum::U8($v) => $body,
            $crate::core::utils::Spectrum::U16($v) => $body,
            $crate::core::utils::Spectrum::U32($v) => $body,
            $crate::core::utils::Spectrum::U64($v) => $body,
        }
    };
}

pub(crate) use spectrum_match;

impl From<Vec<u8>> for Spectrum {
    fn from(other: Vec<u8>) -> Spectrum {
        Spectrum::U8(DVector::from_vec(other))
    }
}

impl From<Vec<u16>> for Spectrum {
    fn from(other: Vec<u16>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u16::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u8)
            ))
        } else {
            Spectrum::U16(DVector::from_vec(other))
        }
    }
}

impl From<Vec<u32>> for Spectrum {
    fn from(other: Vec<u32>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u32::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u8)
            ))
        } else if max <= &u16::MAX.into() {
            Spectrum::U16(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u16)
            ))
        } else {
            Spectrum::U32(DVector::from_vec(other))
        }
    }
}

impl From<Vec<u64>> for Spectrum {
    fn from(other: Vec<u64>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u64::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u8)
            ))
        } else if max <= &u16::MAX.into() {
            Spectrum::U16(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u16)
            ))
        } else if max <= &u16::MAX.into() {
            Spectrum::U32(DVector::from_iterator(
                other.len(), other.iter().map(|&x| x as u32)
            ))
        } else {
            Spectrum::U64(DVector::from_vec(other))
        }
    }
}

impl Spectrum {
    pub fn len(&self) -> usize {
        spectrum_match!(self, spectrum => spectrum.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub enum Spectrum2D {
    U8(DMatrix<u8>),
    U16(DMatrix<u16>),
    U32(DMatrix<u32>),
    U64(DMatrix<u64>),
}

#[macro_export]
macro_rules! spectrum2d_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            $crate::core::utils::Spectrum2D::U8($v) => $body,
            $crate::core::utils::Spectrum2D::U16($v) => $body,
            $crate::core::utils::Spectrum2D::U32($v) => $body,
            $crate::core::utils::Spectrum2D::U64($v) => $body,
        }
    };
}

pub(crate) use spectrum2d_match;

impl Spectrum2D {
    pub fn nrows(&self) -> usize {
        spectrum2d_match!(self, arr => arr.nrows())
    }

    pub fn ncols(&self) -> usize {
        spectrum2d_match!(self, arr => arr.ncols())
    }
}

pub trait FromColMajor<T: Unsigned> {
    fn from_col_major(data: &Vec<T>, nrows: usize, ncols: usize) -> Self;
}

impl FromColMajor<u8> for Spectrum2D {
    fn from_col_major(data: &Vec<u8>, nrows: usize, ncols: usize) -> Self {
        Spectrum2D::U8(DMatrix::from_column_slice(nrows, ncols, data))
    }
}

impl FromColMajor<u16> for Spectrum2D {
    fn from_col_major(data: &Vec<u16>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u16::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_column_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U16(DMatrix::from_column_slice(nrows, ncols, &data))
        }
    }
}

impl FromColMajor<u32> for Spectrum2D {
    fn from_col_major(data: &Vec<u32>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u32::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_column_slice(nrows, ncols, &data))
        } else if max <= &u16::MAX.into() {
            let data = data.iter().map(|&x| x as u16).collect::<Vec<u16>>();
            Spectrum2D::U16(DMatrix::from_column_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U32(DMatrix::from_column_slice(nrows, ncols, &data))
        }
    }
}

impl FromColMajor<u64> for Spectrum2D {
    fn from_col_major(data: &Vec<u64>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u64::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_column_slice(nrows, ncols, &data))
        } else if max <= &u16::MAX.into() {
            let data = data.iter().map(|&x| x as u16).collect::<Vec<u16>>();
            Spectrum2D::U16(DMatrix::from_column_slice(nrows, ncols, &data))
        } else if max <= &u32::MAX.into() {
            let data = data.iter().map(|&x| x as u32).collect::<Vec<u32>>();
            Spectrum2D::U32(DMatrix::from_column_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U64(DMatrix::from_column_slice(nrows, ncols, &data))
        }
    }
}

pub trait FromRowMajor<T: Unsigned> {
    fn from_row_major(data: &Vec<T>, nrows: usize, ncols: usize) -> Self;
}

impl FromRowMajor<u8> for Spectrum2D {
    fn from_row_major(data: &Vec<u8>, nrows: usize, ncols: usize) -> Self {
        Spectrum2D::U8(DMatrix::from_row_slice(nrows, ncols, data))
    }
}

impl FromRowMajor<u16> for Spectrum2D {
    fn from_row_major(data: &Vec<u16>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u16::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_row_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U16(DMatrix::from_row_slice(nrows, ncols, &data))
        }
    }
}

impl FromRowMajor<u32> for Spectrum2D {
    fn from_row_major(data: &Vec<u32>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u32::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_row_slice(nrows, ncols, &data))
        } else if max <= &u16::MAX.into() {
            let data = data.iter().map(|&x| x as u16).collect::<Vec<u16>>();
            Spectrum2D::U16(DMatrix::from_row_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U32(DMatrix::from_row_slice(nrows, ncols, &data))
        }
    }
}

impl FromRowMajor<u64> for Spectrum2D {
    fn from_row_major(data: &Vec<u64>, ncols: usize, nrows: usize) -> Self {
        let max = data.iter().max().unwrap_or(&u64::MAX);

        if max <= &u8::MAX.into() {
            let data = data.iter().map(|&x| x as u8).collect::<Vec<u8>>();
            Spectrum2D::U8(DMatrix::from_row_slice(nrows, ncols, &data))
        } else if max <= &u16::MAX.into() {
            let data = data.iter().map(|&x| x as u16).collect::<Vec<u16>>();
            Spectrum2D::U16(DMatrix::from_row_slice(nrows, ncols, &data))
        } else if max <= &u32::MAX.into() {
            let data = data.iter().map(|&x| x as u32).collect::<Vec<u32>>();
            Spectrum2D::U32(DMatrix::from_row_slice(nrows, ncols, &data))
        } else {
            Spectrum2D::U64(DMatrix::from_row_slice(nrows, ncols, &data))
        }
    }
}


pub trait LossyIntoF64: Copy { fn lossy_into_f64(x: Self) -> f64; }

impl LossyIntoF64 for u8 {
    fn lossy_into_f64(x: u8) -> f64 { x as f64 }
}

impl LossyIntoF64 for u16 {
    fn lossy_into_f64(x: u16) -> f64 { x as f64 }
}

impl LossyIntoF64 for u32 {
    fn lossy_into_f64(x: u32) -> f64 { x as f64 }
}

impl LossyIntoF64 for u64 {
    fn lossy_into_f64(x: u64) -> f64 { x as f64 }
}

impl LossyIntoF64 for f32 {
    fn lossy_into_f64(x: f32) -> f64 { x as f64 }
}

impl LossyIntoF64 for f64 {
    fn lossy_into_f64(x: f64) -> f64 { x }
}
