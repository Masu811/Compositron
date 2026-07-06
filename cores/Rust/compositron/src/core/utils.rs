use thiserror::Error;

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


pub enum Spectrum {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    U64(Vec<u64>),
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
        Spectrum::U8(other)
    }
}

impl From<Vec<u16>> for Spectrum {
    fn from(other: Vec<u16>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u16::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(other.iter().map(|&x| x as u8).collect())
        } else {
            Spectrum::U16(other)
        }
    }
}

impl From<Vec<u32>> for Spectrum {
    fn from(other: Vec<u32>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u32::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(other.iter().map(|&x| x as u8).collect())
        } else if max <= &u16::MAX.into() {
            Spectrum::U16(other.iter().map(|&x| x as u16).collect())
        } else {
            Spectrum::U32(other)
        }
    }
}

impl From<Vec<u64>> for Spectrum {
    fn from(other: Vec<u64>) -> Spectrum {
        let max = other.iter().max().unwrap_or(&u64::MAX);

        if max <= &u8::MAX.into() {
            Spectrum::U8(other.iter().map(|&x| x as u8).collect())
        } else if max <= &u16::MAX.into() {
            Spectrum::U16(other.iter().map(|&x| x as u16).collect())
        } else if max <= &u16::MAX.into() {
            Spectrum::U32(other.iter().map(|&x| x as u32).collect())
        } else {
            Spectrum::U64(other)
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

pub enum SpectrumSlice<'a> {
    U8(&'a [u8]),
    U16(&'a [u16]),
    U32(&'a [u32]),
    U64(&'a [u64]),
}

macro_rules! spectrum_slice_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            SpectrumSlice::U8($v) => $body,
            SpectrumSlice::U16($v) => $body,
            SpectrumSlice::U32($v) => $body,
            SpectrumSlice::U64($v) => $body,
        }
    };
}

pub(crate) use spectrum_slice_match;

impl<'a> From<&'a [u8]> for SpectrumSlice<'a> {
    fn from(other: &'a [u8]) -> SpectrumSlice<'a> {
        SpectrumSlice::U8(other)
    }
}

impl<'a> From<&'a [u16]> for SpectrumSlice<'a> {
    fn from(other: &'a [u16]) -> SpectrumSlice<'a> {
        SpectrumSlice::U16(other)
    }
}

impl<'a> From<&'a [u32]> for SpectrumSlice<'a> {
    fn from(other: &'a [u32]) -> SpectrumSlice<'a> {
        SpectrumSlice::U32(other)
    }
}

impl<'a> From<&'a [u64]> for SpectrumSlice<'a> {
    fn from(other: &'a [u64]) -> SpectrumSlice<'a> {
        SpectrumSlice::U64(other)
    }
}

impl<'a> SpectrumSlice<'a> {
    pub fn len(&self) -> usize {
        spectrum_slice_match!(self, spectrum => spectrum.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct Spectrum2D {
    pub width: usize,
    pub height: usize,
    pub data: Spectrum,
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
