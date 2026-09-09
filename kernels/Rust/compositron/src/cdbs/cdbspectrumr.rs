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
    CEL,
    CML,
}


#[derive(Debug, Clone, Copy)]
pub enum Area {
    AxisAligned {
        width_cel: Unit,
        width_cml: Unit,
    },
    AxisAlignedWithOffset {
        width_cel: Unit,
        width_cml: Unit,
        offset_cel: Unit,
        offset_cml: Unit,
    },
    Diagonal {
        width_first_det: Unit,
        width_second_det: Unit,
        offset_first_det: Unit,
        offset_second_det: Unit,
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
            Area::AxisAligned {
                width_cel: Unit::KeV(2.),
                width_cml: Unit::Eres(1.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
            },
        ],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "W",
        num: &[
            Area::AxisAlignedWithOffset {
                width_cel: Unit::KeV(1.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(3.),
                offset_cml: Unit::KeV(0.),
            },
            Area::AxisAlignedWithOffset {
                width_cel: Unit::KeV(1.),
                width_cml: Unit::Eres(1.),
                offset_cel: Unit::KeV(-3.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
            },
        ],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "P/T",
        num: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::Eres(1.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(f64::INFINITY),
                width_cml: Unit::KeV(f64::INFINITY),
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
    let w1 = width_cel * (
        (1. / ecal.1.scale).powi(2) + (1. / ecal.0.scale).powi(2)
    ).sqrt();
    let w2 = width_cml * (
        (1. / ecal.1.scale).powi(2) + (1. / ecal.0.scale).powi(2)
    ).sqrt();
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
    let w1 = radius_cel * (
        (1. / ecal.1.scale).powi(2) + (1. / ecal.0.scale).powi(2)
    ).sqrt();
    let w2 = radius_cml * (
        (1. / ecal.1.scale).powi(2) + (1. / ecal.0.scale).powi(2)
    ).sqrt();
    Ellipse {
        center_i: ecal.0.to_index_f64(0.) - peak_bnds.0.0 as f64,
        center_j: ecal.1.to_index_f64(0.) - peak_bnds.1.0 as f64,
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


pub struct CDBSpectrumR {
    pub spectrum: Spectrum2D,
    pub detpair: EnergyDetectorPair,
    pub counts: u64,
    pub peak: Option<DMatrix<f64>>,
    pub peak_bnds: Option<((usize, usize), (usize, usize))>,
    pub peak_counts: Option<f64>,
    pub peak_params: HashMap<&'static str, fitting::SimpleFitParam>,
    pub lineshape_params: HashMap<String, LineshapeParam>,
}


impl CDBSpectrumR {
    pub fn new(spectrum: Spectrum2D, detpair: EnergyDetectorPair) -> Self {
        let counts = spectrum2d_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );
        CDBSpectrumR {
            spectrum,
            detpair,
            counts,
            peak: None,
            peak_bnds: None,
            peak_counts: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn project_axes(&self, onto_axis: Axis) -> Projection {
        match onto_axis {
            Axis::CEL => {
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
            Axis::CML => {
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

                fit_gauss2d(&x, &y, peak, &[max, x0, y0, 0.7, 1.8, 0.])?
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

        let first_row_e = ecal_1.to_index_f64(0. - peak_window_kev.0 / 2.);
        let last_row_e = ecal_1.to_index_f64(0. + peak_window_kev.0 / 2.);

        let first_col_e = ecal_2.to_index_f64(0. - peak_window_kev.1 / 2.);
        let last_col_e = ecal_2.to_index_f64(0. + peak_window_kev.1 / 2.);

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

        self.peak_counts = Some(self.peak.as_ref().unwrap().sum());
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
                ecal_1.offset -= y0;
                ecal_2.offset -= x0;
            },
            EcalCorrectionOrder::First => {
                ecal_1.scale *= -ecal_1.offset / (y0 - ecal_1.offset);
                ecal_2.scale *= -ecal_2.offset / (x0 - ecal_2.offset);
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
            Area::AxisAligned { width_cel, width_cml } => {
                let width_cel_kev = width_cel.to_kev(self.detpair.eres)
                    .ok_or(AnalysisError::MissingEnergyResolution)?;
                let width_cml_kev = width_cml.to_kev(self.detpair.eres)
                    .ok_or(AnalysisError::MissingEnergyResolution)?;

                Ok(Polygon { vertices: convert_rectangle(
                    (-width_cml_kev / 2., width_cml_kev / 2.),
                    (-width_cel_kev / 2., width_cel_kev / 2.),
                    ecal,
                    peak_bnds,
                ).to_vertices() })
            },
            Area::AxisAlignedWithOffset { width_cel, width_cml, offset_cel, offset_cml } => todo!(),
            Area::Diagonal { width_first_det, width_second_det, offset_first_det, offset_second_det } => todo!(),
            Area::Ellipse { radius_cel, radius_cml } => todo!(),
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

                    use std::io::Write;

                    let mut f = std::fs::File::create("area.txt").unwrap();

                    for row in view.column_iter() {
                        for x in row.iter() {
                            write!(f, "{x} ").unwrap();
                        }
                        writeln!(f, "").unwrap();
                    }

                    let mut f = std::fs::File::create("weights.txt").unwrap();

                    for row in weights.column_iter() {
                        for x in row.iter() {
                            write!(f, "{x} ").unwrap();
                        }
                        writeln!(f, "").unwrap();
                    }

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
}
