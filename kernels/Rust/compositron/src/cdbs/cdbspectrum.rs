// NOTE:
// - The spectrum is expected to have its origin in the top left corner,
//   meaning energy increases towards higher column indices and higher row
//   indices.
// - The first detector spans the rows (y-axis),
//   the second detector spans the columns (x-axis)
// - CEL: Constant Energy Line
//   CML: Constant Momentum Line
// - Energy along CEL: u := (x - y) / 2  (= actual Doppler shift)
//   Energy along CML: v := (x + y - 1022) / 2
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
use crate::constants::M_E_KEV;
use crate::core::fitting::{self, LMFitError};
use crate::core::utils::{
    spectrum2d_match, EcalCorrectionOrder, EnergyDetector, EnergyDetectorPair,
    LinearCalibration, Spectrum2D, Spectrum, Unit
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

    #[error("Could not fold projection as bins are not symmetric around 0")]
    FoldError,

    #[error(
        "Could not convert projection to DBSpectrum due to missing energy calibration"
    )]
    ProjectionConversionError,
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
pub enum Area {
    Diagonal {
        width_cel: Unit,
        width_cml: Unit,
    },
    DiagonalWithOffset {
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
    pub in_corrected_peak: bool,
}


#[derive(Debug)]
pub struct LineshapeParamDefinition<'a> {
    pub name: &'a str,
    pub num: &'a [Area],
    pub denom: &'a [Area],
    pub in_corrected_peak: bool,
}


pub const STD_LINESHAPE_PARAMS: &[LineshapeParamDefinition; 3] = &[
    LineshapeParamDefinition {
        name: "S",
        num: &[
            Area::Diagonal {
                width_cel: Unit::KeV(2.),
                width_cml: Unit::Eres(1.),
            },
        ],
        denom: &[
            Area::Diagonal {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
            },
        ],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "W",
        num: &[
            Area::DiagonalWithOffset {
                width_cel: Unit::KeV(1.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(3.),
                offset_cml: Unit::KeV(0.),
            },
            Area::DiagonalWithOffset {
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
            },
        ],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "P/T",
        num: &[
            Area::Diagonal {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                first_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
                second_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
            },
        ],
        in_corrected_peak: false,
    },
];


#[inline]
fn ij_to_xy(i: f64, ecal: LinearCalibration) -> f64 {
    ecal.from_index_f64(i)
}


#[inline]
fn xy_to_ij(x: f64, ecal: LinearCalibration) -> f64 {
    ecal.to_index_f64(x)
}


#[inline]
fn xy_to_uv(x: f64, y: f64) -> (f64, f64) {
    (0.5 * (x - y), 0.5 * (x + y) - M_E_KEV)
}


#[inline]
fn uv_to_xy(u: f64, v: f64) -> (f64, f64) {
    (M_E_KEV + u + v, M_E_KEV - u + v)
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
    let w1 = (0.5 * (
        (radius_cel / ecal.1.scale).powi(2) + (radius_cel / ecal.0.scale).powi(2)
    )).sqrt();
    let w2 = (0.5 * (
        (radius_cml / ecal.1.scale).powi(2) + (-radius_cml / ecal.0.scale).powi(2)
    )).sqrt();
    Ellipse {
        center_i: ecal.0.to_index_f64(M_E_KEV) - peak_bnds.0.0 as f64,
        center_j: ecal.1.to_index_f64(M_E_KEV) - peak_bnds.1.0 as f64,
        radius_i: w1,
        radius_j: w2,
        phi: -(ecal.1.scale / ecal.0.scale).atan(),
    }
}


pub struct FoldedProjection {
    pub energies: DVector<f64>,
    pub spectrum: DVector<f64>,
}


pub struct Projection {
    pub spectrum: DVector<f64>,
    pub parent_detector_name: String,
    pub ecal: Option<LinearCalibration>,
    pub bins: Option<DVector<f64>>,
    pub counts: f64,
}


impl Projection {
    pub fn new(
        spectrum: DVector<f64>,
        parent_detector_name: String,
        ecal: Option<LinearCalibration>,
        bins: Option<DVector<f64>>,
    ) -> Self {
        let counts = spectrum.iter().sum::<f64>();

        Projection {
            spectrum,
            parent_detector_name,
            ecal,
            bins,
            counts,
        }
    }


    pub fn get_energies(&self) -> Result<DVector<f64>, AnalysisError> {
        if let Some(ecal) = self.ecal {
            let len = self.spectrum.len();
            Ok(DVector::from_iterator(len, (0..len)
                .map(|i| ecal.from_index(i))
            ))
        } else if let Some(bins) = &self.bins {
            let len = bins.len() - 1;
            Ok(DVector::from_iterator(len, (0..len)
                .map(|i| 0.5 * (bins[i] + bins[i + 1]))
            ))
        } else {
            Err(AnalysisError::FoldError)
        }
    }


    pub fn fold(&self) -> Result<FoldedProjection, AnalysisError> {
        if let Some(ecal) = self.ecal
            && (ecal.to_index_f64(0.) % 1. - 0.5).abs() < 1e-2
        {
            let last_left_idx = ecal.to_index_rounded(0.);
            let first_right_idx = last_left_idx + 1;

            let len_left = last_left_idx + 1;
            let len_right = self.spectrum.len() - len_left;

            let len_new = len_left.min(len_right);

            let spectrum = DVector::from_iterator(
                len_new,
                (0..len_new)
                    .map(|i| {
                        self.spectrum[last_left_idx - i]
                        + self.spectrum[first_right_idx + i]
                    })
            );

            let energies = DVector::from_iterator(
                len_new,
                (first_right_idx..self.spectrum.len())
                    .map(|i| ecal.from_index(i))
            );

            return Ok(FoldedProjection { energies, spectrum });
        } else if let Some(bins) = &self.bins
            && bins.iter().map(|x| x).any(|x| x.abs() < 1e-2)
        {
            let Some(center_idx) = bins
                .iter()
                .enumerate()
                .min_by(|(_, x), (_, y)| x.abs().total_cmp(&y.abs()))
                .map(|(i, _)| i)
            else {
                return Err(AnalysisError::FoldError);
            };

            let len_left = center_idx;
            let len_right = bins.len() - len_left - 1;

            let len_new = len_left.min(len_right);

            let mut spectrum = Vec::new();
            let mut energies = Vec::new();

            for i in 0..len_new - 1 {
                let left_bin_left_edge = bins[center_idx - i - 1];
                let left_bin_right_edge = bins[center_idx - i];

                let right_bin_left_edge = bins[center_idx + i];
                let right_bin_right_edge = bins[center_idx + i + 1];

                if left_bin_left_edge + right_bin_right_edge > 1e-2
                    || left_bin_right_edge + right_bin_left_edge > 1e-2
                {
                    return Err(AnalysisError::FoldError);
                }

                spectrum.push(
                    self.spectrum[center_idx - i - 1]
                    + self.spectrum[center_idx + i]
                );
                energies.push(
                    0.5 * (right_bin_left_edge + right_bin_right_edge) as f64
                );
            }

            return Ok(FoldedProjection {
                energies: DVector::from_vec(energies),
                spectrum: DVector::from_vec(spectrum),
            })
        }

        return Err(AnalysisError::FoldError);
    }


    pub fn to_dbspectrum(&self) -> Result<DBSpectrum, AnalysisError> {
        let Some(ecal) = self.ecal else {
            return Err(AnalysisError::ProjectionConversionError);
        };

        Ok(DBSpectrum::new(
            Spectrum::from(
                self.spectrum.iter().map(|&x| x as u64).collect::<Vec<u64>>()
            ),
            EnergyDetector {
                name: self.parent_detector_name.clone(),
                ecal,
                corrected_ecal: None,
                eres: None,
            }
        ))
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

                let i0 = idx % peak.nrows() as usize + peak_bnds.0.0;
                let j0 = idx / peak.nrows() as usize + peak_bnds.1.0;

                let y0 = ecal_1.from_index(i0);
                let x0 = ecal_2.from_index(j0);

                fit_gauss2d(&x, &y, peak, &[max, x0, y0, 0.7, 1.8, -0.78])?
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

        let first_col = (first_col_e as usize).min(self.spectrum.ncols() - 1);
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
        v: f64,
        ecal: (LinearCalibration, LinearCalibration),
        peak_bnds: ((usize, usize), (usize, usize)),
    ) -> Option<f64> {
        let (i1, j1) = uv_to_ij(0., v, ecal);
        let (i2, j2) = uv_to_ij(0., -v, ecal);
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
            .map(|v| v.0.abs())
            .min_by(|x, y| x.total_cmp(y))
    }


    fn calculate_polygon_boundaries(
        &self,
        area: Area,
        ecal: (LinearCalibration, LinearCalibration),
        peak_bnds: ((usize, usize), (usize, usize)),
    ) -> Result<Polygon, AnalysisError> {
        match area {
            Area::Diagonal { width_cel, width_cml } => {
                let (
                    Some(mut width_cel), Some(width_cml)
                ) = (
                    width_cel.to_kev(self.detpair.eres),
                    width_cml.to_kev(self.detpair.eres),
                ) else {
                    return Err(AnalysisError::MissingEnergyResolution);
                };

                if width_cel.is_infinite() {
                    width_cel = match self.get_max_projection_length(
                        0.5 * width_cml, ecal, peak_bnds
                    ) {
                        Some(x) => 2. * x,
                        None => return Err(AnalysisError::InvalidAreaBounds),
                    };
                }

                if !width_cel.is_finite() || !width_cml.is_finite() {
                    return Err(AnalysisError::InvalidAreaBounds);
                }

                Ok(Polygon { vertices: convert_parallelogram(
                    width_cel, width_cml, 0., 0., ecal, peak_bnds,
                ).to_vertices() })
            },
            Area::DiagonalWithOffset { width_cel, width_cml, offset_cel, offset_cml } => {
                let (
                    Some(mut width_cel), Some(width_cml),
                    Some(offset_cel), Some(offset_cml)
                ) = (
                    width_cel.to_kev(self.detpair.eres),
                    width_cml.to_kev(self.detpair.eres),
                    offset_cel.to_kev(self.detpair.eres),
                    offset_cml.to_kev(self.detpair.eres),
                ) else {
                    return Err(AnalysisError::MissingEnergyResolution);
                };

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
                let radius_cel = radius_cel.to_kev(self.detpair.eres)
                    .ok_or(AnalysisError::MissingEnergyResolution)?;
                let radius_cml = radius_cml.to_kev(self.detpair.eres)
                    .ok_or(AnalysisError::MissingEnergyResolution)?;

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
                .map(|(x, w)| *x * *w as f64)
                .sum::<f64>()
        } else {
            spectrum2d_match!(
                &self.spectrum,
                arr => {
                    let view = arr.view((upper_row, left_col), (nrows, ncols));

                    view
                        .iter()
                        .zip(weights.iter())
                        .map(|(x, w)| *x as u64 * *w as u64)
                        .sum::<u64>() as f64
                }
            )
        } / 255.;

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
            .for_each(|v| *v = Vertex {
                x: v.x - left_col as f64, y: v.y - upper_row as f64
            });

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
    ) -> Result<&LineshapeParam, AnalysisError> {
        let ecal_1 = self.detpair.first_det.corrected_ecal.unwrap_or(
            self.detpair.first_det.ecal
        );
        let ecal_2 = self.detpair.second_det.corrected_ecal.unwrap_or(
            self.detpair.second_det.ecal
        );
        let ecal = (ecal_1, ecal_2);

        let integral_bnds = if definition.in_corrected_peak {
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
            .map(|&area| self.integrate(area, definition.in_corrected_peak))
            .collect::<Result<Vec<f64>, AnalysisError>>()?;

        let denom_areas = definition.denom
            .iter()
            .map(|&area| self.integrate(area, definition.in_corrected_peak))
            .collect::<Result<Vec<f64>, AnalysisError>>()?;

        let mut common_areas = Vec::new();

        for num_area in &num_polygons {
            for denom_area in &denom_polygons {
                let mut intersection = intersect_convex_polygons(
                    num_area, denom_area
                )?;

                if intersection.vertices.len() == 0 { continue; };

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
                    .for_each(|v| *v = Vertex {
                        x: v.x - left_col as f64, y: v.y - upper_row as f64
                    });

                let nrows_view = lower_row - upper_row + 1;
                let ncols_view = right_col - left_col + 1;

                let weights = agg_aa(nrows_view, ncols_view, &intersection)?;

                let integral = self.blend_and_sum(
                    weights,
                    upper_row,
                    left_col,
                    nrows_view,
                    ncols_view,
                    definition.in_corrected_peak,
                )?;

                common_areas.push(integral);
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
            in_corrected_peak: definition.in_corrected_peak,
        };

        self.lineshape_params.insert(definition.name.into(), param);

        Ok(self.lineshape_params.get(definition.name).unwrap())
    }


    pub fn default_analyze(&mut self) -> Result<(), AnalysisError> {
        self.correct_ecal(EcalCorrectionOrder::First)?;

        self.extract_peak((10., 10.), BackgroundModel::None);

        for param in STD_LINESHAPE_PARAMS {
            self.calc_lineshape_param(param)?;
        }

        Ok(())
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

        let width_kev = width.to_kev(self.detpair.eres)
            .ok_or(AnalysisError::MissingEnergyResolution)?;

        let max_length_kev = max_length.to_kev(self.detpair.eres)
            .ok_or(AnalysisError::MissingEnergyResolution)?;

        match bins {
            ProjectionBins::Linear(bins) => {
                let u_bin = bins.to_kev(self.detpair.eres)
                    .ok_or(AnalysisError::MissingEnergyResolution)?;
                let Some(u_max) = self.get_max_projection_length(
                    width_kev / 2., ecal, integral_bnds
                ) else {
                    return Err(AnalysisError::InvalidBins);
                };

                let u_max = u_max.min(max_length_kev);

                let num_bins = 2 * (u_max / u_bin) as usize;

                if num_bins < 2 {
                    return Err(AnalysisError::InvalidBins);
                }

                let max_offset = ((num_bins / 2) as f64 - 0.5) * u_bin;

                let ecal = LinearCalibration {
                    offset: -max_offset, scale: u_bin
                };

                let spectrum = DVector::from_vec((0..num_bins)
                    .map(|i| self.integrate(
                        Area::DiagonalWithOffset {
                            width_cel: Unit::KeV(u_bin),
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
                    .map(|x| x
                        .to_kev(self.detpair.eres)
                        .ok_or(AnalysisError::MissingEnergyResolution)
                    )
                    .collect::<Result<Vec<_>, _>>()?;

                let mut spectrum = DVector::zeros(bins.len() - 1);

                for i in 0..bins_kev.len() - 1 {
                    let left_edge = bins_kev[i];
                    let right_edge = bins_kev[i + 1];

                    let width_cel = right_edge - left_edge;
                    let offset = 0.5 * (left_edge + right_edge);

                    spectrum[i] = self.integrate(
                        Area::DiagonalWithOffset {
                            width_cel: Unit::KeV(width_cel),
                            width_cml: width,
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
