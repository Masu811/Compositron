// NOTE:
// - The spectrum is expected to have its origin in the top left corner,
//   meaning energy increases towards higher column indices and higher row
//   indices.
// - The first detector spans the rows (y-axis),
//   the second detector spans the columns (x-axis)
// - CEL: Constant Energy Line
//   CML: Constant Momentum Line
// - Energy along CEL: u := (x + y - 1022) / sqrt(2)  (= actual Doppler shift)
//   Energy along CML: v := (x - y) / sqrt(2)
// - Energy calibrations and energy resolutions are not assumed to be equal,
//   which means the peak is not assumed to be at 45° in the spectrum, neither
//   in index space (due to ecal assumption) nor energy space (due to eres
//   assumption)

use std::collections::HashMap;

use nalgebra::{DMatrix, DVector};
use thiserror::Error;

use crate::cdbs::anti_aliasing::{
    agg_aa, get_bounding_box, intersect_convex_polygons,
    AntiAliasingError, ConvexToVerticesCounterClockwise,
    Ellipse, Parallelogram, Polygon, Rectangle, Vertex
};
use crate::cdbs::fitting::fit_gauss2d;
use crate::constants::{KEV_PER_M0C, M_E_KEV, ONE_OVER_SQRT_2, TWO_M_E_KEV};
use crate::core::fitting::{self, LMFitError};
use crate::core::utils::{
    spectrum2d_match, EnergyDetectorPair, LinearCalibration, Spectrum2D
};
use crate::dbs::DBSpectrum;

// TODO:
// - Background models
// - Input sanitization
// - Fix: If ecal is very wrong, peak extraction returns empty matrix
// - Fix: Passing inf to areas works only if offset is 0
// - Convenience macro for writing LineshapeParamDefinitions
// - Clipping version of Area and ProjectionBins, that don't error when out of bounds
//
// optimizations for later:
// - coalesce agg calls to reuse buffers
// - use 2D fit results of one spectrum as init for the others
// - couple ecal to DBSpectra, to avoid multiple ecal corrections

#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error(
        "Attempting to integrate an area larger than the spectrum or \
        extracted peak"
    )]
    OutOfBounds,

    #[error(
        "Attempting analysis on background subtracted peak, but no subtraction \
        has been performed yet"
    )]
    NoPeakExtracted,

    #[error("Error during fitting")]
    FitError {
        #[from]
        source: LMFitError,
    },

    #[error("Error during area integration")]
    AntiAliasingError {
        #[from]
        source: AntiAliasingError,
    },

    #[error("Lineshape numerator or denominator areas overlap")]
    OverlappingAreas,

    #[error("Could not convert units of eres to keV due to missing eres")]
    MissingEnergyResolution,

    #[error("Invalid boundaries for area")]
    InvalidAreaBounds,

    #[error("Invalid bins for projection")]
    InvalidBins,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EcalCorrectionOrder {
    Zeroth,
    First,
    None
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BackgroundModel {
    None,
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    FirstDetector,
    SecondDetector,
}

#[derive(Debug, Clone, Copy)]
pub enum Unit {
    KeV(f64),
    M0C(f64),
    Eres(f64),
}

impl Unit {
    pub fn to_kev(&self, eres: Option<f64>) -> Result<f64, AnalysisError> {
        match self {
            Unit::KeV(x) => Ok(*x),
            Unit::M0C(x) => Ok(*x * KEV_PER_M0C),
            Unit::Eres(x) => match eres {
                Some(eres) => Ok(x * eres),
                None => Err(AnalysisError::MissingEnergyResolution),
            },
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Area {
    Diagonal {
        width_cel: Unit,
        width_cml: Unit,
        offset_cel: Unit,
        offset_cml: Unit,
    },
    AxisAligned {
        first_det_bnds: (f64, f64),
        second_det_bnds: (f64, f64),
    },
    Ellipse {
        radius_cel: Unit,
        radius_cml: Unit,
    }
}

#[derive(Debug)]
pub enum ProjectionBins<'a> {
    Linear(Unit),
    CustomZeroCentered(&'a [Unit]),
}

#[derive(Debug, Clone)]
pub struct LineshapeParam {
    pub val: f64,
    pub err: f64,
    pub num: Vec<Area>,
    pub denom: Vec<Area>,
}

#[derive(Debug)]
pub struct LineshapeParamDefinition<'a> {
    pub name: &'a str,
    pub num: &'a [Area],
    pub denom: &'a [Area],
}

pub const STD_LINESHAPE_PARAMS: &[LineshapeParamDefinition; 3] = &[
    LineshapeParamDefinition {
        name: "S",
        num: &[
                Area::Diagonal {
                width_cel: Unit::KeV(2.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        denom: &[
            Area::Diagonal {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
    },
    LineshapeParamDefinition {
        name: "W",
        num: &[
            Area::Diagonal {
                width_cel: Unit::KeV(1.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(3.),
                offset_cml: Unit::KeV(0.),
            },
            Area::Diagonal {
                width_cel: Unit::KeV(1.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(-3.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        denom: &[
            Area::Diagonal {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
    },
    LineshapeParamDefinition {
        name: "P/T",
        num: &[
            Area::Diagonal {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                first_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
                second_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
            },
        ],
    },
];

#[inline]
fn ij_to_xy(i: f64, ecal: LinearCalibration) -> f64 {
    i * ecal.scale + ecal.offset
}

#[inline]
fn xy_to_ij(x: f64, ecal: LinearCalibration) -> f64 {
    (x - ecal.offset) / ecal.scale
}

#[inline]
fn xy_to_uv(x: f64, y: f64) -> (f64, f64) {
    (ONE_OVER_SQRT_2 * (x + y - TWO_M_E_KEV), ONE_OVER_SQRT_2 * (x - y))
}

#[inline]
fn uv_to_xy(u: f64, v: f64) -> (f64, f64) {
    (M_E_KEV + ONE_OVER_SQRT_2 * (u + v), M_E_KEV + ONE_OVER_SQRT_2 * (u - v))
}

#[inline]
fn ij_to_uv(
    i: f64, j: f64, ecal: (LinearCalibration, LinearCalibration)
) -> (f64, f64) {
    let (x, y) = (ij_to_xy(j, ecal.1), ij_to_xy(i, ecal.0));
    xy_to_uv(x, y)
}

#[inline]
fn uv_to_ij(
    u: f64, v: f64, ecal: (LinearCalibration, LinearCalibration)
) -> (f64, f64) {
    let (x, y) = uv_to_xy(u, v);
    (xy_to_ij(y, ecal.0), xy_to_ij(x, ecal.1))
}

fn convert_rectangle(
    first_det_bnds: (f64, f64),
    second_det_bnds: (f64, f64),
    ecal: (LinearCalibration, LinearCalibration),
    peak_bnds: ((usize, usize), (usize, usize)),
) -> Rectangle {
    Rectangle {
        i_min: ecal.0.to_index_f64(first_det_bnds.1) - peak_bnds.0.0 as f64,
        i_max: ecal.0.to_index_f64(first_det_bnds.0) - peak_bnds.0.0 as f64,
        j_min: ecal.1.to_index_f64(second_det_bnds.0) - peak_bnds.1.0 as f64,
        j_max: ecal.1.to_index_f64(second_det_bnds.1) - peak_bnds.1.0 as f64,
    }
}

fn convert_parallelogram(
    width_cel: f64,
    width_cml: f64,
    offset_cel: f64,
    offset_cml: f64,
    ecal: (LinearCalibration, LinearCalibration),
    peak_bnds: ((usize, usize), (usize, usize)),
) -> Parallelogram {
    let (i, j) = uv_to_ij(offset_cel, offset_cml, ecal);
    let w1 = (0.5 * (
        (width_cel / ecal.1.scale).powi(2) + (width_cel / ecal.0.scale).powi(2)
    )).sqrt();
    let w2 = (0.5 * (
        (width_cml / ecal.1.scale).powi(2) + (-width_cml / ecal.0.scale).powi(2)
    )).sqrt();
    Parallelogram {
        center_i: i - peak_bnds.0.0 as f64,
        center_j: j - peak_bnds.1.0 as f64,
        slope_1: -ecal.1.scale / ecal.0.scale,
        slope_2: ecal.1.scale / ecal.0.scale,
        side_length_1: w1,
        side_length_2: w2,
    }
}

fn convert_ellipse(
    radius_cel: f64,
    radius_cml: f64,
    ecal: (LinearCalibration, LinearCalibration),
    peak_bnds: ((usize, usize), (usize, usize)),
) -> Ellipse {
    let w1 = ONE_OVER_SQRT_2 * ecal.1.to_index_f64(radius_cel + radius_cml);
    let w2 = ONE_OVER_SQRT_2 * ecal.0.to_index_f64(radius_cel - radius_cml);
    Ellipse {
        center_i: ecal.0.to_index_f64(M_E_KEV) - peak_bnds.0.0 as f64,
        center_j: ecal.1.to_index_f64(M_E_KEV) - peak_bnds.1.0 as f64,
        radius_i: w1,
        radius_j: w2,
        phi: (ecal.0.scale / ecal.1.scale).atan(),
    }
}

pub struct Projection {
    pub spectrum: DVector<f64>,
    pub parent_detector_name: String,
    pub ecal: Option<LinearCalibration>,
    pub bins: Option<DVector<f64>>,
    pub counts: f64,
    pub dcounts: f64,
}

impl Projection {
    pub fn new(
        spectrum: DVector<f64>,
        parent_detector_name: String,
        ecal: Option<LinearCalibration>,
        bins: Option<DVector<f64>>,
    ) -> Self {
        let counts = spectrum.iter().sum::<f64>();
        let dcounts = counts.sqrt();

        Projection {
            spectrum,
            parent_detector_name,
            ecal,
            bins,
            counts,
            dcounts,
        }
    }

    pub fn fold() -> Result<(DVector<f64>, DVector<f64>), AnalysisError> {
        // Check: are the bins symmetric around 0
        //     yes => fold it and ship it
        //     no => return Err
        todo!()
    }

    pub fn to_dbspectrum() -> Result<DBSpectrum, AnalysisError> {
        // Check: are the bins linear
        //     yes => return DBSpectrum
        //     no => return Err
        todo!()
    }
}

pub struct CDBSpectrum {
    pub spectrum: Spectrum2D,
    pub detpair: EnergyDetectorPair,
    pub counts: u64,
    pub dcounts: f64,
    pub peak: Option<DMatrix<f64>>,
    pub peak_bnds: Option<((usize, usize), (usize, usize))>,
    pub peak_counts: Option<f64>,
    pub dpeak_counts: Option<f64>,
    pub peak_params: HashMap<&'static str, fitting::SimpleFitParam>,
    pub lineshape_params: HashMap<String, LineshapeParam>,
}

impl CDBSpectrum {
    pub fn new(
        spectrum: Spectrum2D, detpair: EnergyDetectorPair,
    ) -> Self {
        let counts = spectrum2d_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );
        CDBSpectrum {
            spectrum,
            detpair,
            counts,
            dcounts: (counts as f64).sqrt(),
            peak: None,
            peak_bnds: None,
            peak_counts: None,
            dpeak_counts: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn new_default_analyzed(
        spectrum: Spectrum2D, detector: EnergyDetectorPair,
    ) -> Result<Self, AnalysisError> {
        let mut s = CDBSpectrum::new(spectrum, detector);

        s.correct_ecal(EcalCorrectionOrder::First)?;

        s.extract_peak((10., 10.), BackgroundModel::None);

        for param in STD_LINESHAPE_PARAMS {
            s.calc_lineshape_param(param, true)?;
        }

        Ok(s)
    }

    pub fn project_axes(&self, onto_axis: Axis) -> Projection {
        match onto_axis {
            Axis::FirstDetector => {
                let spectrum = spectrum2d_match!(
                    &self.spectrum,
                    arr => DVector::from_iterator(arr.nrows(), arr
                        .row_iter()
                        .map(|row| row.iter().map(|&x| x as f64).sum())
                    )
                );

                Projection::new(
                    spectrum,
                    self.detpair.first_det.name.clone(),
                    Some(self.detpair.first_det.ecal),
                    None,
                )
            },
            Axis::SecondDetector => {
                let spectrum = spectrum2d_match!(
                    &self.spectrum,
                    arr => DVector::from_iterator(arr.nrows(), arr
                        .row_iter()
                        .map(|row| row.iter().map(|&x| x as f64).sum())
                    )
                );

                Projection::new(
                    spectrum,
                    self.detpair.second_det.name.clone(),
                    Some(self.detpair.second_det.ecal),
                    None,
                )
            },
        }
    }

    fn fit_2d_peak(
        &mut self, bg_model: BackgroundModel
    ) -> Result<(), AnalysisError> {
        let Some(peak) = &self.peak else {
            return Err(AnalysisError::NoPeakExtracted);
        };

        let Some(peak_bnds) = self.peak_bnds else {
            return Err(AnalysisError::NoPeakExtracted);
        };

        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );

        let nrows = peak_bnds.0.1 - peak_bnds.0.0;
        let ncols = peak_bnds.1.1 - peak_bnds.1.0;

        let x = DVector::from_iterator(
            ncols, (peak_bnds.1.0..peak_bnds.1.1).map(|i| ecal_2.from_index(i))
        );

        let y = DVector::from_iterator(
            nrows, (peak_bnds.0.0..peak_bnds.0.1).map(|i| ecal_1.from_index(i))
        );

        self.peak_params = match bg_model {
            BackgroundModel::None => {
                let (idx, &max) = peak
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap();

                let x0 = idx % peak.ncols() as usize + peak_bnds.1.0;
                let y0 = idx / peak.ncols() as usize + peak_bnds.0.0;

                let y0 = ecal_1.from_index(y0);
                let x0 = ecal_2.from_index(x0);

                fit_gauss2d(&x, &y, peak, &[max, x0, y0, 1., 2., -0.78])?
            },
        };

        Ok(())
    }

    fn subtract_bg(&mut self, bg_model: BackgroundModel) {
        if bg_model == BackgroundModel::None {
            return;
        }

        match bg_model {
            BackgroundModel::None => unreachable!(),
        }
    }

    pub fn extract_peak(
        &mut self,
        peak_window_kev: (f64, f64),
        bg_model: BackgroundModel,
    ) {
        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );

        let first_row_e = ecal_1.to_index_f64(M_E_KEV - peak_window_kev.0 / 2.);
        let last_row_e = ecal_1.to_index_f64(M_E_KEV + peak_window_kev.0 / 2.);

        let first_col_e = ecal_2.to_index_f64(M_E_KEV - peak_window_kev.1 / 2.);
        let last_col_e = ecal_2.to_index_f64(M_E_KEV + peak_window_kev.1 / 2.);

        let first_row = (first_row_e as usize).min(self.spectrum.nrows() - 1);
        let last_row = ((last_row_e + 1.) as usize).min(self.spectrum.nrows() - 1);

        let first_col = (first_col_e as usize).min(self.spectrum.nrows() - 1);
        let last_col = ((last_col_e + 1.) as usize).min(self.spectrum.ncols() - 1);

        let nrows = last_row - first_row;
        let ncols = last_col - first_col;

        self.peak = Some(spectrum2d_match!(&self.spectrum, arr => {
            let view = arr.view((first_row, first_col), (nrows, ncols));
            DMatrix::<f64>::from_fn(nrows, ncols, |i, j| {view[(i, j)] as f64})
        }));

        self.peak_bnds = Some(((first_row, last_row), (first_col, last_col)));

        self.subtract_bg(bg_model);

        let peak_counts = self.peak.as_ref().unwrap().sum();

        self.peak_counts = Some(peak_counts);
        self.dpeak_counts = Some(peak_counts.sqrt());
    }

    pub fn correct_ecal(
        &mut self, order: EcalCorrectionOrder,
    ) -> Result<(), AnalysisError> {
        if order == EcalCorrectionOrder::None {
            return Ok(());
        }

        if self.peak == None {
            self.extract_peak((4., 4.), BackgroundModel::None);
        }

        if !self.peak_params.contains_key("x0")
            || !self.peak_params.contains_key("y0")
        {
            self.fit_2d_peak(BackgroundModel::None)?;
        }

        let x0 = self.peak_params.get("x0").unwrap().val;
        let y0 = self.peak_params.get("y0").unwrap().val;

        let mut ecal_1 = self.detpair.first_det.ecal;
        let mut ecal_2 = self.detpair.second_det.ecal;

        match order {
            EcalCorrectionOrder::Zeroth => {
                ecal_1.offset += M_E_KEV - y0;
                ecal_2.offset += M_E_KEV - x0;
            },
            EcalCorrectionOrder::First => {
                ecal_1.scale *= (M_E_KEV - ecal_1.offset) / (y0 - ecal_1.offset);
                ecal_2.scale *= (M_E_KEV - ecal_2.offset) / (x0 - ecal_2.offset);
            },
            EcalCorrectionOrder::None => unreachable!(),
        }

        self.detpair.first_det.corrected_ecal = Some(ecal_1);
        self.detpair.second_det.corrected_ecal = Some(ecal_2);

        Ok(())
    }

    fn get_max_projection_length(
        &self,
        u: f64,
        ecal: (LinearCalibration, LinearCalibration),
        peak_bnds: ((usize, usize), (usize, usize)),
    ) -> Option<f64> {
        let (i1, j1) = uv_to_ij(u, 0., ecal);
        let (i2, j2) = uv_to_ij(-u, 0., ecal);
        let m = -ecal.1.scale / ecal.0.scale;
        let t1 = i1 - m * j1;
        let t2 = i2 - m * j2;

        let ((i_min, i_max), (j_min, j_max)) = peak_bnds;

        let i_min = i_min as f64 - 0.49;
        let i_max = i_max as f64 + 0.49;
        let j_min = j_min as f64 - 0.49;
        let j_max = j_max as f64 + 0.49;

        let x1 = ij_to_uv(i_min, (i_min - t1) / m, ecal);
        let x2 = ij_to_uv(i_min, (i_min - t2) / m, ecal);
        let x3 = ij_to_uv(m * j_min + t1, j_min, ecal);
        let x4 = ij_to_uv(m * j_min + t2, j_min, ecal);
        let x5 = ij_to_uv(i_max, (i_max - t1) / m, ecal);
        let x6 = ij_to_uv(i_max, (i_max - t2) / m, ecal);
        let x7 = ij_to_uv(m * j_max + t1, j_max, ecal);
        let x8 = ij_to_uv(m * j_max + t2, j_max, ecal);

        [x1, x2, x3, x4, x5, x6, x7, x8]
            .iter()
            .map(|v| v.1.abs())
            .min_by(|x, y| x.total_cmp(y))
    }

    fn calculate_polygon_boundaries(
        &self,
        area: Area,
        ecal: (LinearCalibration, LinearCalibration),
        peak_bnds: ((usize, usize), (usize, usize)),
    ) -> Result<Polygon, AnalysisError> {
        match area {
            Area::Diagonal { width_cel, width_cml, offset_cel, offset_cml } => {
                let mut width_cel = width_cel.to_kev(self.detpair.eres)?;
                let width_cml = width_cml.to_kev(self.detpair.eres)?;
                let offset_cel = offset_cel.to_kev(self.detpair.eres)?;
                let offset_cml = offset_cml.to_kev(self.detpair.eres)?;

                if width_cel.is_infinite() {
                    width_cel = match self.get_max_projection_length(
                        0.5 * width_cml, ecal, peak_bnds
                    ) {
                        Some(x) => 2. * x,
                        None => return Err(AnalysisError::InvalidAreaBounds),
                    };
                }

                if !width_cel.is_finite() || !width_cml.is_finite() ||
                    !offset_cel.is_finite() || !offset_cml.is_finite()
                {
                    return Err(AnalysisError::InvalidAreaBounds);
                }

                Ok(Polygon { vertices: convert_parallelogram(
                    width_cel, width_cml, offset_cel, offset_cml, ecal, peak_bnds,
                ).to_vertices() })
            },
            Area::AxisAligned { first_det_bnds, second_det_bnds } => {
                let (mut y_min, mut y_max) = first_det_bnds;
                let (mut x_min, mut x_max) = second_det_bnds;

                if y_min.is_infinite() {
                    y_min = ecal.0.from_index_f64(peak_bnds.0.0 as f64 - 0.49);
                }
                if y_max.is_infinite() {
                    y_max = ecal.0.from_index_f64(peak_bnds.0.1 as f64 + 0.49);
                }
                if x_min.is_infinite() {
                    x_min = ecal.1.from_index_f64(peak_bnds.1.0 as f64 - 0.49);
                }
                if x_max.is_infinite() {
                    x_max = ecal.1.from_index_f64(peak_bnds.1.1 as f64 + 0.49);
                }

                if !y_min.is_finite() || !y_max.is_finite() ||
                    !x_min.is_finite() || !x_max.is_finite() {
                    return Err(AnalysisError::InvalidAreaBounds);
                }

                Ok(Polygon { vertices: convert_rectangle(
                    (y_min, y_max), (x_min, x_max), ecal, peak_bnds,
                ).to_vertices() })
            },
            Area::Ellipse { radius_cel, radius_cml } => {
                let radius_cel = radius_cel.to_kev(self.detpair.eres)?;
                let radius_cml = radius_cml.to_kev(self.detpair.eres)?;

                if !radius_cel.is_finite() || !radius_cml.is_finite() {
                    return Err(AnalysisError::InvalidAreaBounds);
                }

                Ok(Polygon { vertices: convert_ellipse(
                    radius_cel, radius_cml, ecal, peak_bnds,
                ).to_vertices() })
            },
        }
    }

    fn blend_and_sum(
        &mut self,
        weights: DMatrix<u8>,
        upper_row: usize,
        left_col: usize,
        nrows: usize,
        ncols: usize,
        in_corrected_peak: bool,
    ) -> Result<f64, AnalysisError> {
        let integral = if in_corrected_peak {
            let Some(peak) = &self.peak else {
                return Err(AnalysisError::NoPeakExtracted);
            };

            let view = peak.view((upper_row, left_col), (nrows, ncols));

            view
                .iter()
                .zip(weights.iter())
                .map(|(x, w)| *x as u64 * *w as u64)
                .sum::<u64>()
        } else {
            spectrum2d_match!(
                &self.spectrum,
                arr => {
                    let view = arr.view((upper_row, left_col), (nrows, ncols));

                    view
                        .iter()
                        .zip(weights.iter())
                        .map(|(x, w)| *x as u64 * *w as u64)
                        .sum()
                }
            )
        } as f64 / 255.;

        Ok(integral)
    }

    pub fn integrate(
        &mut self,
        area: Area,
        in_corrected_peak: bool,
    ) -> Result<f64, AnalysisError> {
        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );
        let ecal = (ecal_1, ecal_2);

        let integral_bnds = if in_corrected_peak {
            let Some(bnds) = self.peak_bnds else {
                return Err(AnalysisError::NoPeakExtracted);
            };
            bnds
        } else {
            ((0, self.spectrum.nrows() - 1), (0, self.spectrum.ncols() - 1))
        };

        let nrows = integral_bnds.0.1 - integral_bnds.0.0 + 1;
        let ncols = integral_bnds.1.1 - integral_bnds.1.0 + 1;

        let mut polygon = self.calculate_polygon_boundaries(
            area, ecal, integral_bnds,
        )?;

        let (
            left_col, right_col, lower_row, upper_row
        ) = get_bounding_box(&polygon.to_vertices())?;

        if right_col > ncols - 1 || lower_row > nrows - 1 {
            return Err(AnalysisError::AntiAliasingError {
                source: AntiAliasingError::OutOfBounds
            });
        }

        polygon.vertices
            .iter_mut()
            .for_each(|v| *v = Vertex(
                v.0 - left_col as f64, v.1 - upper_row as f64
            ));

        let nrows_view = lower_row - upper_row + 1;
        let ncols_view = right_col - left_col + 1;

        let weights = agg_aa(nrows_view, ncols_view, &polygon)?;

        self.blend_and_sum(
            weights,
            upper_row,
            left_col,
            nrows_view,
            ncols_view,
            in_corrected_peak,
        )
    }

    pub fn calc_lineshape_param(
        &mut self,
        definition: &LineshapeParamDefinition,
        in_corrected_peak: bool,
    ) -> Result<&LineshapeParam, AnalysisError> {
        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );
        let ecal = (ecal_1, ecal_2);

        let integral_bnds = if in_corrected_peak {
            let Some(bnds) = self.peak_bnds else {
                return Err(AnalysisError::NoPeakExtracted);
            };
            bnds
        } else {
            ((0, self.spectrum.nrows() - 1), (0, self.spectrum.ncols() - 1))
        };

        let nrows = integral_bnds.0.1 - integral_bnds.0.0 + 1;
        let ncols = integral_bnds.1.1 - integral_bnds.1.0 + 1;

        let num_polygons = definition.num
            .iter()
            .map(|&area| self.calculate_polygon_boundaries(
                area, ecal, integral_bnds,
            ))
            .collect::<Result<Vec<_>, _>>()?;

        let denom_polygons = definition.denom
            .iter()
            .map(|&area| self.calculate_polygon_boundaries(
                area, ecal, integral_bnds,
            ))
            .collect::<Result<Vec<_>, _>>()?;

        for field in [&num_polygons, &denom_polygons] {
            for i in 0..field.len() - 1 {
                for j in i + 1..field.len() {
                    let intersection = intersect_convex_polygons(
                        &field[i], &field[j]
                    )?;

                    if intersection.vertices.len() > 0 {
                        return Err(AnalysisError::OverlappingAreas);
                    }
                }
            }
        }

        let num_areas = definition.num
            .iter()
            .map(|&area| self.integrate(area, in_corrected_peak))
            .collect::<Result<Vec<f64>, AnalysisError>>()?;

        let denom_areas = definition.denom
            .iter()
            .map(|&area| self.integrate(area, in_corrected_peak))
            .collect::<Result<Vec<f64>, AnalysisError>>()?;

        let mut common_areas = Vec::new();

        for num_area in &num_polygons {
            for denom_area in &denom_polygons {
                let mut intersection = intersect_convex_polygons(
                    num_area, denom_area
                )?;

                if intersection.vertices.len() > 0 {
                    let (
                        left_col, right_col, lower_row, upper_row
                    ) = get_bounding_box(&intersection.to_vertices())?;

                    if right_col > ncols - 1 || lower_row > nrows - 1 {
                        return Err(AnalysisError::AntiAliasingError {
                            source: AntiAliasingError::OutOfBounds
                        });
                    }

                    intersection.vertices
                        .iter_mut()
                        .for_each(|v| *v = Vertex(
                            v.0 - left_col as f64, v.1 - upper_row as f64
                        ));

                    let nrows_view = lower_row - upper_row + 1;
                    let ncols_view = right_col - left_col + 1;

                    let weights = agg_aa(nrows_view, ncols_view, &intersection)?;

                    let integral = self.blend_and_sum(
                        weights,
                        upper_row,
                        left_col,
                        nrows_view,
                        ncols_view,
                        in_corrected_peak,
                    )?;

                    common_areas.push(integral);
                }
            }
        }

        let num = num_areas.iter().sum::<f64>();
        let denom = denom_areas.iter().sum::<f64>();
        let common = common_areas.iter().sum::<f64>();

        let val = num / denom;
        let err = val * (
            (num - common) / num.powi(2) +
            (denom - common) / denom.powi(2) +
            common * ((denom - num) / (denom * num)).powi(2)
        ).sqrt();

        let param = LineshapeParam {
            val,
            err,
            num: definition.num.to_owned(),
            denom: definition.denom.to_owned(),
        };

        self.lineshape_params.insert(definition.name.into(), param);

        Ok(self.lineshape_params.get(definition.name).unwrap())
    }

    pub fn project_digonal(
        &mut self,
        bins: ProjectionBins,
        width: Unit,
        max_length: Unit,
        in_corrected_peak: bool,
    ) -> Result<Projection, AnalysisError> {
        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );
        let ecal = (ecal_1, ecal_2);

        let integral_bnds = if in_corrected_peak {
            let Some(bnds) = self.peak_bnds else {
                return Err(AnalysisError::NoPeakExtracted);
            };
            bnds
        } else {
            ((0, self.spectrum.nrows() - 1), (0, self.spectrum.ncols() - 1))
        };

        match bins {
            ProjectionBins::Linear(bins) => {
                let v_bin = bins.to_kev(self.detpair.eres)?;
                let Some(v_max) = self.get_max_projection_length(
                    width.to_kev(self.detpair.eres)? / 2., ecal, integral_bnds
                ) else {
                    return Err(AnalysisError::InvalidBins);
                };

                let v_max = v_max.min(max_length.to_kev(self.detpair.eres)?);

                let num_bins = 2 * (v_max / v_bin) as usize;

                if num_bins < 2 {
                    return Err(AnalysisError::InvalidBins);
                }

                let max_offset = ((num_bins / 2) as f64 - 0.5) * v_bin;

                let ecal = LinearCalibration {
                    offset: -max_offset, scale: v_bin
                };

                let spectrum = DVector::from_vec((0..num_bins)
                    .map(|i| self.integrate(
                        Area::Diagonal {
                            width_cel: Unit::KeV(v_bin),
                            width_cml: width,
                            offset_cel: Unit::KeV(ecal.from_index(i)),
                            offset_cml: Unit::KeV(0.),
                        },
                        in_corrected_peak
                    ))
                    .collect::<Result<Vec<_>, _>>()?
                );

                Ok(Projection::new(
                    spectrum,
                    self.detpair.name.clone(),
                    Some(ecal),
                    None,
                ))
            },
            ProjectionBins::CustomZeroCentered(bins) => {
                if bins.len() < 2 {
                    return Err(AnalysisError::InvalidBins);
                }

                let bins_kev = bins
                    .iter()
                    .map(|x| x.to_kev(self.detpair.eres))
                    .collect::<Result<Vec<_>, _>>()?;

                let mut spectrum = DVector::zeros(bins.len() - 1);

                for i in 0..bins_kev.len() - 1 {
                    let left_edge = bins_kev[i];
                    let right_edge = bins_kev[i + 1];

                    let width_cel = right_edge - left_edge;
                    let offset = 0.5 * (left_edge + right_edge);

                    spectrum[i] = self.integrate(
                        Area::Diagonal {
                            width_cel: Unit::KeV(width_cel),
                            width_cml: Unit::KeV(0.),
                            offset_cel: Unit::KeV(offset),
                            offset_cml: Unit::KeV(0.),
                        },
                        in_corrected_peak
                    )?;
                }

                Ok(Projection::new(
                    spectrum,
                    self.detpair.name.clone(),
                    None,
                    Some(DVector::from_column_slice(&bins_kev))
                ))
            },
        }
    }
}
