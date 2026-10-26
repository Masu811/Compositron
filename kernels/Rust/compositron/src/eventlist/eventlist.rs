use std::collections::HashMap;
use std::fmt::Debug;
use std::ops::Sub;
use std::sync::Arc;

use num_traits::AsPrimitive;
use thiserror::Error;

use crate::cdbs::CDBSpectrum;
use crate::core::utils::{
    Detector, EnergyDetector, EnergyDetectorPair, FromColMajor, Spectrum, Spectrum2D, TimingDetectorPair
};
use crate::dbs::DBSpectrum;
use crate::pals::PALSpectrum;
use crate::eventlist::coinc_eval::{self, CoincEvalMode};


#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("No column named '{name}'")]
    ColumnNotFound {
        name: String,
    },

    #[error("No detector named '{name}'")]
    DetectorNotFound {
        name: String,
    },

    #[error("Eventlist columns have different numbers of rows")]
    DifferentColumnLengths,

    #[error("Eventlist column has length 0")]
    ZeroLengthError,

    #[error(
        "Incompatible coincidence column types: column '{col_a_name}' has type \
        {col_a_type:?}, but column '{col_b_name}' has type {col_b_type:?}"
    )]
    IncompatibleColumnTypes {
        col_a_name: String,
        col_a_type: EventListColumnType,
        col_b_name: String,
        col_b_type: EventListColumnType,
    },

    #[error("Number of bins exceeds the supported maximum of 2^16 - 2")]
    TooManyBins,

    #[error("Number of bins is zero")]
    NoBins,

    #[error("Invalid bin edges: min > max")]
    InvalidBinEdges,
}


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EventListColumnType {
    U8,
    U16,
    U32,
    U64,
    I8,
    I16,
    I32,
    I64,
    F32,
    F64,
}


#[derive(Debug, Clone)]
pub enum EventListColumn {
    U8(Arc<[u8]>),
    U16(Arc<[u16]>),
    U32(Arc<[u32]>),
    U64(Arc<[u64]>),
    I8(Arc<[i8]>),
    I16(Arc<[i16]>),
    I32(Arc<[i32]>),
    I64(Arc<[i64]>),
    F32(Arc<[f32]>),
    F64(Arc<[f64]>),
}


#[macro_export]
macro_rules! elistcol_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            $crate::eventlist::eventlist::EventListColumn::U8($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::U16($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::U32($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::U64($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::I8($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::I16($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::I32($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::I64($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::F32($v) => $body,
            $crate::eventlist::eventlist::EventListColumn::F64($v) => $body,
        }
    };
}


#[macro_export]
macro_rules! split_elistcol_match {
    ($s:expr, $v:ident => $body_v:expr, $w:ident => $body_w:expr) => {
        match $s {
            $crate::eventlist::eventlist::EventListColumn::U8($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::U16($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::U32($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::U64($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::I8($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::I16($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::I32($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::I64($v) => $body_v,
            $crate::eventlist::eventlist::EventListColumn::F32($w) => $body_w,
            $crate::eventlist::eventlist::EventListColumn::F64($w) => $body_w,
        }
    };
}


#[derive(Debug, Clone)]
pub enum OwnedEventListColumn {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    U64(Vec<u64>),
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}


#[macro_export]
macro_rules! owned_elistcol_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            $crate::eventlist::eventlist::OwnedEventListColumn::U8($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::U16($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::U32($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::U64($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::I8($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::I16($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::I32($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::I64($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::F32($v) => $body,
            $crate::eventlist::eventlist::OwnedEventListColumn::F64($v) => $body,
        }
    };
}


#[macro_export]
macro_rules! split_owned_elistcol_match {
    ($s:expr, $v:ident => $body_v:expr, $w:ident => $body_w:expr) => {
        match $s {
            $crate::eventlist::eventlist::OwnedEventListColumn::U8($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::U16($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::U32($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::U64($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::I8($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::I16($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::I32($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::I64($v) => $body_v,
            $crate::eventlist::eventlist::OwnedEventListColumn::F32($w) => $body_w,
            $crate::eventlist::eventlist::OwnedEventListColumn::F64($w) => $body_w,
        }
    };
}


#[derive(Debug, Clone, Copy)]
pub struct Bins {
    bin_width: Option<f64>,
    len: Option<u16>,
    min_bin_edge: Option<f64>,
    max_bin_edge: Option<f64>,
}


impl Bins {
    pub fn from_bin_width(
        bin_width: f64, min_bin_edge: Option<f64>, max_bin_edge: Option<f64>,
    ) -> Result<Self, AnalysisError> {
        let len = if let (Some(min), Some(max)) = (min_bin_edge, max_bin_edge) {
            let len_float = (max - min) / bin_width;
            if len_float < 0. {
                return Err(AnalysisError::InvalidBinEdges);
            }
            if len_float >= u16::MAX as f64 {
                return Err(AnalysisError::TooManyBins);
            }
            Some(len_float as u16)
        } else { None };

        Ok(Bins {
            bin_width: Some(bin_width),
            len,
            min_bin_edge,
            max_bin_edge,
        })
    }

    pub fn from_len(
        len: u16, min_bin_edge: Option<f64>, max_bin_edge: Option<f64>,
    ) -> Result<Self, AnalysisError> {
        if len == u16::MAX {
            return Err(AnalysisError::TooManyBins);
        }
        let bin_width = if let (Some(min), Some(max)) = (min_bin_edge, max_bin_edge) {
            let bin_width = (max - min) / len as f64;
            if bin_width < 0. {
                return Err(AnalysisError::InvalidBinEdges);
            }
            Some(bin_width)
        } else { None };

        Ok(Bins {
            bin_width,
            len: Some(len),
            min_bin_edge,
            max_bin_edge,
        })
    }

    pub fn len(&self) -> Option<u16> {
        self.len
    }


    pub fn accomodate_vec<T: AsPrimitive<f64> + PartialOrd>(
        &mut self, arr: &[T], active_indices: &[usize],
    ) -> Result<(), AnalysisError> {
        if self.min_bin_edge.is_none() {
            let Some(min) = active_indices
                .iter()
                .map(|&i| arr[i])
                .min_by(|a, b| {
                    a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                })
            else {
                return Err(AnalysisError::ZeroLengthError);
            };
            self.min_bin_edge = Some(min.as_());
        }

        if self.max_bin_edge.is_none() {
            let Some(max) = active_indices
                .iter()
                .map(|&i| arr[i])
                .max_by(|a, b| {
                    a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                })
            else {
                return Err(AnalysisError::ZeroLengthError);
            };
            self.max_bin_edge = Some(max.as_());
        }

        *self = if let Some(bin_width) = self.bin_width {
            Bins::from_bin_width(
                bin_width, self.min_bin_edge, self.max_bin_edge
            )?
        } else {
            Bins::from_len(
                self.len.unwrap(), self.min_bin_edge, self.max_bin_edge
            )?
        };

        Ok(())
    }


    pub fn accomodate_col(
        &mut self, col: &EventListColumn, active_indices: &[usize],
    ) -> Result<(), AnalysisError> {
        if self.min_bin_edge.is_none() {
            split_elistcol_match!(
                col,
                arr_int => {
                    let Some(min) = active_indices
                        .iter()
                        .map(|&i| arr_int[i])
                        .min()
                    else {
                        return Err(AnalysisError::ZeroLengthError);
                    };
                    self.min_bin_edge = Some(min as f64);
                },
                arr_float => {
                    let Some(min) = active_indices
                        .iter()
                        .map(|&i| arr_float[i])
                        .min_by(|a, b| a.total_cmp(b))
                    else {
                        return Err(AnalysisError::ZeroLengthError);
                    };
                    self.min_bin_edge = Some(min as f64);
                }
            );
        }

        if self.max_bin_edge.is_none() {
            split_elistcol_match!(
                col,
                arr_int => {
                    let Some(max) = active_indices
                        .iter()
                        .map(|&i| arr_int[i])
                        .max()
                    else {
                        return Err(AnalysisError::ZeroLengthError);
                    };
                    self.max_bin_edge = Some(max as f64);
                },
                arr_float => {
                    let Some(max) = active_indices
                        .iter()
                        .map(|&i| arr_float[i])
                        .max_by(|a, b| a.total_cmp(b))
                    else {
                        return Err(AnalysisError::ZeroLengthError);
                    };
                    self.max_bin_edge = Some(max as f64);
                }
            );
        }

        *self = if let Some(bin_width) = self.bin_width {
            Bins::from_bin_width(
                bin_width, self.min_bin_edge, self.max_bin_edge
            )?
        } else {
            Bins::from_len(
                self.len.unwrap(), self.min_bin_edge, self.max_bin_edge
            )?
        };

        Ok(())
    }


    pub fn vec_to_bins<T: AsPrimitive<f64> + PartialOrd>(
        &mut self, arr: &[T], active_indices: &[usize],
    ) -> Result<Vec<u16>, AnalysisError> {
        if self.min_bin_edge.is_none() || self.max_bin_edge.is_none() || self.len.is_none() {
            self.accomodate_vec(arr, active_indices)?;
        }
        let x_min = self.min_bin_edge.unwrap();
        let dx = self.bin_width.unwrap();
        let nbins = self.len.unwrap();
        let dx_inv = 1. / dx;

        if nbins == 0 {
            return Err(AnalysisError::NoBins);
        }

        let max_bin = (nbins - 1).min(u16::MAX) as f64;

        let mut bins = Vec::with_capacity(arr.len());

        for &i in active_indices {
            let bin = ((arr[i]).as_() - x_min) * dx_inv;
            if bin < 0. || bin > max_bin {
                bins.push(u16::MAX);
            } else {
                bins.push(bin as u16);
            }
        }

        Ok(bins)
    }


    pub fn col_to_bins(
        &mut self, col: &EventListColumn, active_indices: &[usize],
    ) -> Result<Vec<u16>, AnalysisError> {
        if self.min_bin_edge.is_none() || self.max_bin_edge.is_none() || self.len.is_none() {
            self.accomodate_col(col, active_indices)?;
        }
        let x_min = self.min_bin_edge.unwrap();
        let dx = self.bin_width.unwrap();
        let nbins = self.len.unwrap();
        let dx_inv = 1. / dx;

        if nbins == 0 {
            return Err(AnalysisError::NoBins);
        }

        let max_bin = nbins.min(u16::MAX) as f64;

        elistcol_match!(
            col,
            arr => {
                let mut bins = Vec::with_capacity(arr.len());

                for &i in active_indices {
                    let bin = (arr[i] as f64 - x_min) * dx_inv;
                    if bin < 0. || bin >= max_bin {
                        bins.push(u16::MAX);
                    } else {
                        bins.push(bin as u16);
                    }
                }

                Ok(bins)
            }
        )
    }
}


impl EventListColumn {
    pub fn len(&self) -> usize {
        elistcol_match!(self, arr => arr.len())
    }

    pub fn column_type(&self) -> EventListColumnType {
        match self {
            EventListColumn::U8(_) => EventListColumnType::U8,
            EventListColumn::U16(_) => EventListColumnType::U16,
            EventListColumn::U32(_) => EventListColumnType::U32,
            EventListColumn::U64(_) => EventListColumnType::U64,
            EventListColumn::I8(_) => EventListColumnType::I8,
            EventListColumn::I16(_) => EventListColumnType::I16,
            EventListColumn::I32(_) => EventListColumnType::I32,
            EventListColumn::I64(_) => EventListColumnType::I64,
            EventListColumn::F32(_) => EventListColumnType::F32,
            EventListColumn::F64(_) => EventListColumnType::F64,
        }
    }
}


#[derive(Debug, Clone)]
pub struct EventList {
    pub detector: Detector,
    pub columns: HashMap<String, EventListColumn>,
    active_indices: Vec<usize>,
    arc_len: usize,
}


#[derive(Debug, Clone)]
pub struct CorrelatedEventLists {
    pub eventlists: HashMap<String, EventList>,
}


#[derive(Debug, Clone)]
pub struct CoincidenceEventList {
    single_columns: HashMap<String, EventList>,
    shared_columns: EventList,
    shared_active_indices: Vec<usize>,
}


impl EventList {
    pub fn try_new_with_active_indices(
        detector: Detector,
        columns: HashMap<String, EventListColumn>,
        active_indices: Vec<usize>,
    ) -> Result<Self, AnalysisError> {
        let nrows = active_indices.len();

        for col in columns.values() {
            if col.len() != nrows {
                return Err(AnalysisError::DifferentColumnLengths);
            }
        }

        Ok(EventList {
            detector,
            columns,
            active_indices,
            arc_len: nrows,
        })
    }


    pub fn try_new(
        detector: Detector, columns: HashMap<String, EventListColumn>,
    ) -> Result<Self, AnalysisError> {
        let ncols = columns.len();

        let nrows = if ncols == 0 {
            0
        } else {
            columns.values().next().unwrap().len()
        };

        EventList::try_new_with_active_indices(
            detector, columns, (0..nrows).collect()
        )
    }


    pub fn ncols(&self) -> usize {
        self.columns.len()
    }


    pub fn nrows(&self) -> usize {
        self.active_indices.len()
    }


    pub fn add_column(
        &mut self, name: &str, column: EventListColumn
    ) -> Result<Option<EventListColumn>, AnalysisError> {
        if self.columns.len() == 0 {
            Ok(self.columns.insert(name.into(), column))
        } else if column.len() == self.arc_len {
            Ok(self.columns.insert(name.into(), column))
        } else {
            Err(AnalysisError::DifferentColumnLengths)
        }
    }


    pub fn to_histogram(
        &self, column: &str, bins: &mut Bins,
    ) -> Result<Spectrum, AnalysisError> {
        let Some(col) = self.columns.get(column) else {
            return Err(AnalysisError::ColumnNotFound { name: column.into() });
        };

        bins.accomodate_col(col, &self.active_indices)?;

        let nbins = bins.len().unwrap() as usize;

        let bins = bins.col_to_bins(col, &self.active_indices)?;

        let mut hist: Vec<u64> = Vec::new();
        hist.resize(nbins, 0);

        for idx in bins.into_iter() {
            if idx < u16::MAX {
                hist[idx as usize] += 1;
            }
        }

        return Ok(Spectrum::from(hist))
    }


    pub fn to_2d_histogram(
        &self,
        col_param_colname: &str,
        row_param_colname: &str,
        col_bins: &mut Bins,
        row_bins: &mut Bins,
    ) -> Result<Spectrum2D, AnalysisError> {
        let Some(col) = self.columns.get(col_param_colname) else {
            return Err(AnalysisError::ColumnNotFound { name: col_param_colname.into() });
        };
        let Some(row) = self.columns.get(row_param_colname) else {
            return Err(AnalysisError::ColumnNotFound { name: row_param_colname.into() });
        };

        col_bins.accomodate_col(col, &self.active_indices)?;
        row_bins.accomodate_col(row, &self.active_indices)?;

        let ncols = col_bins.len().unwrap() as usize;
        let nrows = row_bins.len().unwrap() as usize;

        let col_bins = col_bins.col_to_bins(col, &self.active_indices)?;
        let row_bins = row_bins.col_to_bins(row, &self.active_indices)?;

        let mut hist: Vec<u64> = Vec::new();
        hist.resize(nrows * ncols, 0);

        for (col_idx, row_idx) in col_bins.into_iter().zip(row_bins) {
            if col_idx < u16::MAX && row_idx < u16::MAX {
                hist[col_idx as usize * nrows + row_idx as usize] += 1;
            }
        }

        return Ok(Spectrum2D::from_col_major(&hist, ncols, nrows))
    }


    pub fn to_dbspectrum(
        &self,
        energy_colname: &str,
        bins: &mut Bins,
    ) -> Result<DBSpectrum, AnalysisError> {
        let spectrum = self.to_histogram(energy_colname, bins)?;

        Ok(DBSpectrum::new(
            spectrum,
            EnergyDetector {
                name: self.detector.name.clone(),
                ecal: self.detector.ecal,
                corrected_ecal: None,
                eres: self.detector.eres,
            },
        ))
    }


    pub fn to_cdbspectrum(
        &self,
        col_energy_colname: &str,
        row_energy_colname: &str,
        col_bins: &mut Bins,
        row_bins: &mut Bins,
        detpair: EnergyDetectorPair,
    ) -> Result<CDBSpectrum, AnalysisError> {
        todo!();
    }


    pub fn to_palspectrum(
        &self,
        lifetime_colname: &str,
        bins: &mut Bins,
        detpair: TimingDetectorPair,
    ) -> Result<PALSpectrum, AnalysisError> {
        todo!();
    }


    pub fn to_acarspectrum() {
        todo!();
    }


    pub fn to_amocspectrum() {
        todo!();
    }
}


impl CorrelatedEventLists {
    pub fn new(eventlists: HashMap<String, EventList>) -> Self {
        CorrelatedEventLists { eventlists }
    }


    fn to_hist<T: Copy + Sub<Output = T> + PartialOrd>(
        tarr_c: &[T],
        tarr_r: &[T],
        iarr_c: &[u16],
        iarr_r: &[u16],
        active_indices_c: &[usize],
        active_indices_r: &[usize],
        hist: &mut [u64],
        nrows: usize,
        coinc_time_window: T,
        coinc_eval_mode: CoincEvalMode,
    ) {
        match coinc_eval_mode {
            CoincEvalMode::EachMatch => {
                coinc_eval::to_hist_each_match(
                    tarr_c, tarr_r, iarr_c, iarr_r,
                    active_indices_c, active_indices_r,
                    hist, nrows, coinc_time_window,
                );
            },
            CoincEvalMode::FirstMatch => {
                coinc_eval::to_hist_first_match(
                    tarr_c, tarr_r, iarr_c, iarr_r,
                    active_indices_c, active_indices_r,
                    hist, nrows, coinc_time_window,
                );
            },
            CoincEvalMode::ClosestMatch => {
                coinc_eval::to_hist_closest_match(
                    tarr_c, tarr_r, iarr_c, iarr_r,
                    active_indices_c, active_indices_r,
                    hist, nrows, coinc_time_window,
                );
            },
            CoincEvalMode::UniqueMatch => {
                coinc_eval::to_hist_unique_match(
                    tarr_c, tarr_r, iarr_c, iarr_r,
                    active_indices_c, active_indices_r,
                    hist, nrows, coinc_time_window,
                );
            },
            CoincEvalMode::FirstPivot => {
                coinc_eval::to_hist_first_pivot(
                    tarr_c, tarr_r, iarr_c, iarr_r,
                    active_indices_c, active_indices_r,
                    hist, nrows, coinc_time_window,
                );
            },
        }
    }


    fn to_celist<T: Copy + Sub<Output = T> + PartialOrd>(
        tarr_c: &[T],
        tarr_r: &[T],
        active_indices_c: &[usize],
        active_indices_r: &[usize],
        coinc_time_window: T,
        coinc_eval_mode: CoincEvalMode,
    ) -> (Vec<usize>, Vec<usize>) {
        match coinc_eval_mode {
            CoincEvalMode::EachMatch => {
                coinc_eval::to_celist_each_match(
                    tarr_c, tarr_r, active_indices_c, active_indices_r,
                    coinc_time_window,
                )
            },
            CoincEvalMode::FirstMatch => {
                coinc_eval::to_celist_first_match(
                    tarr_c, tarr_r, active_indices_c, active_indices_r,
                    coinc_time_window,
                )
            },
            CoincEvalMode::ClosestMatch => {
                coinc_eval::to_celist_closest_match(
                    tarr_c, tarr_r, active_indices_c, active_indices_r,
                    coinc_time_window,
                )
            },
            CoincEvalMode::UniqueMatch => {
                coinc_eval::to_celist_unique_match(
                    tarr_c, tarr_r, active_indices_c, active_indices_r,
                    coinc_time_window,
                )
            },
            CoincEvalMode::FirstPivot => {
                coinc_eval::to_celist_first_pivot(
                    tarr_c, tarr_r, active_indices_c, active_indices_r,
                    coinc_time_window,
                )
            },
        }
    }


    pub fn to_2d_coinc_histogram(
        &self,
        col_detname: &str,
        row_detname: &str,
        col_det_timestamp_colname: &str,
        row_det_timestamp_colname: &str,
        col_det_param_colname: &str,
        row_det_param_colname: &str,
        col_det_bins: &mut Bins,
        row_det_bins: &mut Bins,
        coinc_time_window: f64,
        coinc_eval_mode: CoincEvalMode,
    ) -> Result<Spectrum2D, AnalysisError> {
        let Some(col_det_elist) = self.eventlists.get(col_detname) else {
            return Err(AnalysisError::DetectorNotFound {
                name: col_detname.into()
            });
        };
        let Some(row_det_elist) = self.eventlists.get(row_detname) else {
            return Err(AnalysisError::DetectorNotFound {
                name: col_detname.into()
            });
        };

        let Some(timestamps_col_det) = col_det_elist.columns.get(
            col_det_timestamp_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: col_det_timestamp_colname.into()
            });
        };
        let Some(timestamps_row_det) = row_det_elist.columns.get(
            row_det_timestamp_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: row_det_timestamp_colname.into()
            });
        };

        let timestamp_type = timestamps_col_det.column_type();

        if timestamps_row_det.column_type() != timestamp_type {
            return Err(AnalysisError::IncompatibleColumnTypes {
                col_a_name: col_det_timestamp_colname.into(),
                col_a_type: timestamps_col_det.column_type(),
                col_b_name: row_det_timestamp_colname.into(),
                col_b_type: timestamps_row_det.column_type(),
            });
        }

        let Some(params_col_det) = col_det_elist.columns.get(
            col_det_param_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: col_det_param_colname.into()
            });
        };
        let Some(params_row_det) = row_det_elist.columns.get(
            row_det_param_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: row_det_param_colname.into()
            });
        };

        let param_type = params_col_det.column_type();

        if params_row_det.column_type() != param_type {
            return Err(AnalysisError::IncompatibleColumnTypes {
                col_a_name: col_det_param_colname.into(),
                col_a_type: params_col_det.column_type(),
                col_b_name: row_det_param_colname.into(),
                col_b_type: params_row_det.column_type(),
            });
        }

        let act_idcs_c = &col_det_elist.active_indices;
        let act_idcs_r = &row_det_elist.active_indices;

        col_det_bins.accomodate_col(params_col_det, act_idcs_c)?;
        row_det_bins.accomodate_col(params_row_det, act_idcs_r)?;

        let ncols = col_det_bins.len().unwrap() as usize;
        let nrows = row_det_bins.len().unwrap() as usize;

        let mut hist: Vec<u64> = Vec::new();
        hist.resize(ncols * nrows, 0);

        let iarr_c = col_det_bins.col_to_bins(params_col_det, act_idcs_c)?;
        let iarr_r = row_det_bins.col_to_bins(params_row_det, act_idcs_r)?;

        macro_rules! match_timestamp_types {
            ($($ctype:ident),* $(,)?) => {
                match (timestamps_col_det, timestamps_row_det) {
                    $(
                        (EventListColumn::$ctype(tarr_c), EventListColumn::$ctype(tarr_r)) => {
                            CorrelatedEventLists::to_hist(
                                tarr_c, tarr_r, &iarr_c, &iarr_r,
                                act_idcs_c, act_idcs_r, &mut hist, nrows,
                                coinc_time_window.as_(), coinc_eval_mode,
                            );
                        },
                    )*
                    _ => unreachable!(),
                }
            };
        }

        match_timestamp_types!(U8, U16, U32, U64, I8, I16, I32, I64, F32, F64);

        Ok(Spectrum2D::from_col_major(&hist, ncols, nrows))
    }


    pub fn to_coinc_eventlist(
        &self,
        first_detname: &str,
        second_detname: &str,
        detpair: Detector,
        first_det_timestamp_colname: &str,
        second_det_timestamp_colname: &str,
        coinc_time_window: f64,
        coinc_eval_mode: CoincEvalMode,
    ) -> Result<CoincidenceEventList, AnalysisError> {
        let Some(first_det_elist) = self.eventlists.get(first_detname) else {
            return Err(AnalysisError::DetectorNotFound {
                name: first_detname.into()
            });
        };
        let Some(second_det_elist) = self.eventlists.get(second_detname) else {
            return Err(AnalysisError::DetectorNotFound {
                name: second_detname.into()
            });
        };

        let Some(timestamps_first_det) = first_det_elist.columns.get(
            first_det_timestamp_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: first_det_timestamp_colname.into()
            });
        };
        let Some(timestamps_second_det) = second_det_elist.columns.get(
            second_det_timestamp_colname
        ) else {
            return Err(AnalysisError::ColumnNotFound {
                name: second_det_timestamp_colname.into()
            });
        };

        let timestamp_type = timestamps_first_det.column_type();

        if timestamps_second_det.column_type() != timestamp_type {
            return Err(AnalysisError::IncompatibleColumnTypes {
                col_a_name: first_det_timestamp_colname.into(),
                col_a_type: timestamps_first_det.column_type(),
                col_b_name: second_det_timestamp_colname.into(),
                col_b_type: timestamps_second_det.column_type(),
            });
        }

        let act_idcs_1 = &first_det_elist.active_indices;
        let act_idcs_2 = &second_det_elist.active_indices;

        macro_rules! match_timestamp_types {
            ($($ctype:ident),* $(,)?) => {
                match (timestamps_first_det, timestamps_second_det) {
                    $(
                        (EventListColumn::$ctype(tarr_1), EventListColumn::$ctype(tarr_2)) => {
                            CorrelatedEventLists::to_celist(
                                tarr_1, tarr_2, act_idcs_1, act_idcs_2,
                                coinc_time_window.as_(), coinc_eval_mode,
                            )
                        },
                    )*
                    _ => unreachable!(),
                }
            };
        }

        let (first_det_idcs, second_det_idcs) = match_timestamp_types!(
            U8, U16, U32, U64, I8, I16, I32, I64, F32, F64
        );

        let mut single_columns = HashMap::new();

        single_columns.insert(
            first_detname.into(),
            EventList {
                detector: first_det_elist.detector.clone(),
                columns: first_det_elist.columns.clone(),
                active_indices: first_det_idcs,
                arc_len: first_det_elist.arc_len,
            }
        );

        single_columns.insert(
            second_detname.into(),
            EventList {
                detector: second_det_elist.detector.clone(),
                columns: second_det_elist.columns.clone(),
                active_indices: second_det_idcs,
                arc_len: second_det_elist.arc_len,
            }
        );

        Ok(CoincidenceEventList {
            single_columns,
            shared_columns: EventList {
                detector: detpair,
                columns: HashMap::new(),
                active_indices: Vec::new(),
                arc_len: 0,
            },
            shared_active_indices: Vec::new(),
        })
    }


    pub fn to_cdbspectrum(
        &self,
        col_detname: &str,
        row_detname: &str,
        col_det_timestamp_colname: &str,
        row_det_timestamp_colname: &str,
        col_det_energy_colname: &str,
        row_det_energy_colname: &str,
        col_bins: Bins,
        row_bins: Bins,
        coinc_time_window: f64,
        coinc_eval_mode: CoincEvalMode,
    ) -> Result<CDBSpectrum, AnalysisError> {
        todo!();
    }


    pub fn to_palspectrum(
        &self,
        start_detname: &str,
        stop_detname: &str,
        start_det_timestamp_colname: &str,
        stop_det_timestamp_colname: &str,
        start_det_energy_colname: &str,
        stop_det_energy_colname: &str,
        start_energy_window_first_det: Option<(f64, f64)>,
        stop_energy_window_first_det: Option<(f64, f64)>,
        start_energy_window_second_det: Option<(f64, f64)>,
        stop_energy_window_second_det: Option<(f64, f64)>,
        bins: Bins,
        coinc_eval_mode: CoincEvalMode,
    ) -> PALSpectrum {
        todo!();
    }


    pub fn to_acarspectrum(&self) {
        todo!();
    }


    pub fn to_amocspectrum(&self) {
        todo!();
    }
}


impl CoincidenceEventList {
    pub fn ncols(&self) -> usize {
        let mut ncols = self.shared_columns.ncols();

        for elist in self.single_columns.values() {
            ncols += elist.ncols();
        }

        ncols
    }


    pub fn nrows(&self) -> usize {
        self.shared_active_indices.len()
    }


    pub fn to_histogram(
        &self, detname: Option<&str>, colname: &str, bins: &mut Bins,
    ) -> Result<Spectrum, AnalysisError> {
        let elist = match detname {
            Some(name) => self.single_columns.get(name).ok_or(
                AnalysisError::DetectorNotFound { name: name.into() }
            )?,
            None => &self.shared_columns,
        };

        let Some(col) = elist.columns.get(colname) else {
            return Err(AnalysisError::ColumnNotFound { name: colname.into() });
        };

        bins.accomodate_col(col, &elist.active_indices)?;

        let nbins = bins.len().unwrap() as usize;

        let bins = bins.col_to_bins(col, &elist.active_indices)?;

        let mut hist: Vec<u64> = Vec::new();
        hist.resize(nbins, 0);

        for idx in bins.into_iter() {
            if idx < u16::MAX {
                hist[idx as usize] += 1;
            }
        }

        return Ok(Spectrum::from(hist))
    }


    pub fn to_2d_histogram(
        &self,
        col_detname: Option<&str>,
        row_detname: Option<&str>,
        col_param_colname: &str,
        row_param_colname: &str,
        col_bins: &mut Bins,
        row_bins: &mut Bins,
    ) -> Result<Spectrum2D, AnalysisError> {
        let col_elist = match col_detname {
            Some(name) => self.single_columns.get(name).ok_or(
                AnalysisError::DetectorNotFound { name: name.into() }
            )?,
            None => &self.shared_columns,
        };
        let row_elist = match row_detname {
            Some(name) => self.single_columns.get(name).ok_or(
                AnalysisError::DetectorNotFound { name: name.into() }
            )?,
            None => &self.shared_columns,
        };

        let Some(col) = col_elist.columns.get(col_param_colname) else {
            return Err(AnalysisError::ColumnNotFound { name: col_param_colname.into() });
        };
        let Some(row) = row_elist.columns.get(row_param_colname) else {
            return Err(AnalysisError::ColumnNotFound { name: row_param_colname.into() });
        };

        col_bins.accomodate_col(col, &col_elist.active_indices)?;
        row_bins.accomodate_col(row, &row_elist.active_indices)?;

        let ncols = col_bins.len().unwrap() as usize;
        let nrows = row_bins.len().unwrap() as usize;

        let col_bins = col_bins.col_to_bins(col, &col_elist.active_indices)?;
        let row_bins = row_bins.col_to_bins(row, &row_elist.active_indices)?;

        let mut hist: Vec<u64> = Vec::new();
        hist.resize(nrows * ncols, 0);

        for (col_idx, row_idx) in col_bins.into_iter().zip(row_bins) {
            if col_idx < u16::MAX && row_idx < u16::MAX {
                hist[col_idx as usize * nrows + row_idx as usize] += 1;
            }
        }

        return Ok(Spectrum2D::from_col_major(&hist, ncols, nrows))
    }


    pub fn to_dbspectrum(
        &self,
        detname: &str,
        energy_colname: &str,
        bins: &mut Bins,
    ) -> Result<DBSpectrum, AnalysisError> {
        let Some(elist) = self.single_columns.get(detname) else {
            return Err(AnalysisError::DetectorNotFound { name: detname.into() });
        };

        let spectrum = elist.to_histogram(energy_colname, bins)?;

        Ok(DBSpectrum::new(
            spectrum,
            EnergyDetector {
                name: elist.detector.name.clone(),
                ecal: elist.detector.ecal,
                corrected_ecal: None,
                eres: elist.detector.eres,
            },
        ))
    }


    pub fn to_cdbspectrum(
        &self,
        col_energy_colname: &str,
        row_energy_colname: &str,
        col_bins: &mut Bins,
        row_bins: &mut Bins,
        detpair: EnergyDetectorPair,
    ) -> Result<CDBSpectrum, AnalysisError> {
        todo!();
    }


    pub fn to_palspectrum(
        &self,
        lifetime_colname: &str,
        bins: &mut Bins,
        detpair: TimingDetectorPair,
    ) -> Result<PALSpectrum, AnalysisError> {
        todo!();
    }


    pub fn to_acarspectrum() {
        todo!();
    }


    pub fn to_amocspectrum() {
        todo!();
    }
}
