use std::{fs::File, io::Read};

use num_traits::Unsigned;
use nalgebra::{DMatrix, DVector};

use crate::constants::KEV_PER_M0C;


#[derive(Debug, Clone, Copy)]
pub struct LinearCalibration {
    pub offset: f64,
    pub scale: f64,
}

impl LinearCalibration {
    pub fn from_index(&self, index: usize) -> f64 {
        self.scale * index as f64 + self.offset
    }

    pub fn from_index_f64(&self, index: f64) -> f64 {
        self.scale * index + self.offset
    }

    pub fn to_index_rounded(&self, value: f64) -> usize {
        ((value - self.offset) / self.scale).round() as usize
    }

    pub fn to_index_trunc(&self, value: f64) -> usize {
        ((value - self.offset) / self.scale) as usize
    }

    pub fn to_index_f64(&self, value: f64) -> f64 {
        (value - self.offset) / self.scale
    }
}


#[derive(Debug, Clone)]
pub struct Detector {
    pub name: String,
    pub ecal: Option<LinearCalibration>,
    pub tcal: Option<LinearCalibration>,
    pub eres: Option<f64>,
    pub tres: Option<f64>,
}

impl From<&str> for Detector {
    fn from(value: &str) -> Self {
        Detector {
            name: value.into(), ecal: None, tcal: None, eres: None, tres: None,
        }
    }
}

impl From<&String> for Detector {
    fn from(value: &String) -> Self {
        Detector {
            name: value.into(), ecal: None, tcal: None, eres: None, tres: None,
        }
    }
}

impl From<String> for Detector {
    fn from(value: String) -> Self {
        Detector {
            name: value.into(), ecal: None, tcal: None, eres: None, tres: None,
        }
    }
}

// DBS
#[derive(Debug, Clone)]
pub struct EnergyDetector {
    pub name: String,
    pub ecal: Option<LinearCalibration>,
    pub corrected_ecal: Option<LinearCalibration>,
    pub eres: Option<f64>,
}

// CDBS
#[derive(Debug, Clone)]
pub struct EnergyDetectorPair {
    pub name: String,
    pub first_det: EnergyDetector,
    pub second_det: EnergyDetector,
    pub eres: Option<f64>,
}

// PALS
#[derive(Debug, Clone)]
pub struct TimingDetectorPair {
    pub name: String,
    pub tcal: Option<LinearCalibration>,
    pub corrected_tcal: Option<LinearCalibration>,
    pub tres: Option<f64>,
}

// AMOC
#[derive(Debug, Clone)]
pub struct AMOCDetectorTriplet {
    pub name: String,
    pub energy_det: EnergyDetector,
    pub timing_detpair: TimingDetectorPair,
}


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EcalCorrectionOrder {
    Zeroth,
    First,
    None
}

#[derive(Debug, Clone, Copy)]
pub enum Unit {
    KeV(f64),
    M0C(f64),
    Eres(f64),
}

impl Unit {
    pub fn to_kev(&self, eres: Option<f64>) -> Option<f64> {
        match self {
            Unit::KeV(x) => Some(*x),
            Unit::M0C(x) => Some(*x * KEV_PER_M0C),
            Unit::Eres(x) => eres.map(|e| e * x)
        }
    }
}

#[derive(Debug)]
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


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ByteOrder {
    LittleEndian,
    BigEndian,
}


pub trait FromBytes<const N: usize> {
    fn from_bytes(bytes: [u8; N], order: ByteOrder) -> Self;
}


macro_rules! impl_from_bytes {
    ($target_type:ident, $size:expr) => {
        impl FromBytes<$size> for $target_type {
            fn from_bytes(bytes: [u8; $size], order: ByteOrder) -> Self {
                match order {
                    ByteOrder::LittleEndian => $target_type::from_le_bytes(bytes),
                    ByteOrder::BigEndian => $target_type::from_be_bytes(bytes),
                }
            }
        }
    };
}


impl_from_bytes!(u8, 1);
impl_from_bytes!(u16, 2);
impl_from_bytes!(u32, 4);
impl_from_bytes!(u64, 8);
impl_from_bytes!(u128, 16);
impl_from_bytes!(i8, 1);
impl_from_bytes!(i16, 2);
impl_from_bytes!(i32, 4);
impl_from_bytes!(i64, 8);
impl_from_bytes!(i128, 16);
impl_from_bytes!(f32, 4);
impl_from_bytes!(f64, 8);


pub fn read_num<const N: usize, T: FromBytes<N>>(
    file: &mut impl Read,
    order: ByteOrder,
) -> Result<T, std::io::Error> {
    let mut buffer = [0u8; N];
    file.read_exact(&mut buffer)?;
    Ok(T::from_bytes(buffer, order))
}


pub fn read_num_from_bytes<const N: usize, T: FromBytes<N>>(
    buffer: &[u8], order: ByteOrder,
) -> Result<T, std::io::Error> {
    Ok(T::from_bytes(buffer.try_into().ok().unwrap(), order))
}


pub fn read_str<const N_BYTES: usize>(
    file: &mut impl Read,
) -> Result<String, std::io::Error> {
    let mut buffer = [0u8; N_BYTES];
    file.read_exact(&mut buffer)?;
    Ok(String::from_utf8_lossy(&buffer).into())
}


pub fn read_str_dyn(
    file: &mut impl Read, n_bytes: usize,
) -> Result<String, std::io::Error> {
    let mut buffer = Vec::new();
    buffer.resize(n_bytes, 0u8);
    file.read_exact(&mut buffer)?;
    Ok(String::from_utf8_lossy(&buffer).into())
}
