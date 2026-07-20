use std::collections::HashMap;
use std::f64;

use nalgebra::{DMatrix, DVector};
use thiserror::Error;

use crate::cdbs::anti_aliasing::AntiAliasingScheme;
use crate::cdbs::fitting::fit_gauss2d;
use crate::constants::M_E_KEV;
use crate::core::fitting::{self, LMFitError};
use crate::core::utils::{spectrum2d_match, Spectrum, Spectrum2D};
use crate::dbs::DBSpectrum;

// TODO:
// - Input sanitization
// - Initial guesses for fits
// - Detector names for projection
// - If ecal is very wrong, peak extraction returns empty matrix

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
    // Erf,
    // ErfGauss,
    // ErfLinear1Gauss,
    // ErfLinear2Gauss,
    // ErfGauss1Exp,
    // ErfGauss2Exp,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BinIntegrationScheme {
    Const,
    Bilinear,
}

#[derive(Debug, Clone, Copy)]
pub enum Axis {
    FirstDetector,
    SecondDetector,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PeakLineshapeParamEvalSide {
    Both,
    Left,
    Right,
}

#[derive(Debug, Clone, Copy)]
pub enum Area {
    Diagonal {
        width_cel: f64,
        width_cml: f64,
        offset_cel: f64,
        offset_cml: f64,
        rotation: f64,
    },
    AxisAligned {
        first_det_bnds: (f64, f64),
        second_det_bnds: (f64, f64),
    },
}

#[derive(Debug)]
pub enum LineshapeParam {
    PeakSym {
        val: f64,
        err: f64,
        num: Area,
        denom: Area,
    },
    PeakAsym {
        val: f64,
        err: f64,
        num: Area,
        denom: Area,
        num_side: PeakLineshapeParamEvalSide,
        denom_side: PeakLineshapeParamEvalSide,
    },
    Spectrum {
        val: f64,
        err: f64,
        num: Area,
        denom: Area,
    },
}

#[derive(Debug)]
pub enum LineshapeParamDefinition<'a> {
    PeakSym {
        name: &'a str,
        num: Area,
        denom: Area,
    },
    PeakAsym {
        name: &'a str,
        num: Area,
        denom: Area,
        num_side: PeakLineshapeParamEvalSide,
        denom_side: PeakLineshapeParamEvalSide,
    },
    Spectrum {
        name: &'a str,
        num: Area,
        denom: Area,
    },
}

pub const STD_LINESHAPE_PARAMS: &[LineshapeParamDefinition; 3] = &[
    LineshapeParamDefinition::PeakSym {
        name: "S",
        num: Area::Diagonal {
            width_cel: 1.,
            width_cml: f64::NAN,
            offset_cel: 0.,
            offset_cml: 0.,
            rotation: f64::NAN,
        },
        denom: Area::Diagonal {
            width_cel: f64::INFINITY,
            width_cml: f64::NAN,
            offset_cel: 0.,
            offset_cml: 0.,
            rotation: f64::NAN,
        },
    },
    LineshapeParamDefinition::PeakAsym {
        name: "W",
        num: Area::Diagonal {
            width_cel: 1.,
            width_cml: f64::NAN,
            offset_cel: 3.,
            offset_cml: 0.,
            rotation: f64::NAN,
        },
        denom: Area::Diagonal {
            width_cel: f64::INFINITY,
            width_cml: f64::NAN,
            offset_cel: 0.,
            offset_cml: 0.,
            rotation: f64::NAN,
        },
        num_side: PeakLineshapeParamEvalSide::Both,
        denom_side: PeakLineshapeParamEvalSide::Both,
    },
    LineshapeParamDefinition::PeakSym {
        name: "P/T",
        denom: Area::Diagonal {
            width_cel: f64::INFINITY,
            width_cml: f64::NAN,
            offset_cel: 0.,
            offset_cml: 0.,
            rotation: f64::NAN,
        },
        num: Area::AxisAligned {
            first_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
            second_det_bnds: (f64::NEG_INFINITY, f64::INFINITY),
        },
    },
];

pub struct CDBSpectrum {
    pub spectrum: Spectrum2D,
    pub detpair: String,
    pub ecal: ((f64, f64), (f64, f64)),
    pub corrected_ecal: Option<((f64, f64), (f64, f64))>,
    pub eres: Option<f64>,
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
        spectrum: Spectrum2D, detpair: String, ecal: ((f64, f64), (f64, f64))
    ) -> Self {
        let counts = spectrum2d_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );
        CDBSpectrum {
            spectrum: spectrum,
            detpair: detpair,
            ecal: ecal,
            corrected_ecal: None,
            eres: None,
            counts: counts,
            dcounts: (counts as f64).sqrt(),
            peak: None,
            peak_bnds: None,
            peak_counts: None,
            dpeak_counts: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn project_axes(&self, onto_axis: Axis) -> DBSpectrum {
        let spectrum = match onto_axis {
            Axis::FirstDetector => spectrum2d_match!(
                &self.spectrum,
                arr => arr
                    .row_iter()
                    .map(|row| row.iter().map(|&x| x as u64).sum())
                    .collect::<Vec<u64>>()
            ),
            Axis::SecondDetector => spectrum2d_match!(
                &self.spectrum,
                arr => arr
                    .column_iter()
                    .map(|col| col.iter().map(|&x| x as u64).sum())
                    .collect::<Vec<u64>>()
            ),
        };

        DBSpectrum::new(Spectrum::from(spectrum), "A".into(), (0., 1.), None)
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

        let ecal = self.corrected_ecal.unwrap_or(self.ecal);

        let nrows = peak_bnds.0.1 - peak_bnds.0.0;
        let ncols = peak_bnds.1.1 - peak_bnds.1.0;

        let x = DVector::from_iterator(
            ncols, (peak_bnds.1.0..peak_bnds.1.1).map(|i| {
                ecal.1.0 + ecal.1.1 * i as f64
            })
        );

        let y = DVector::from_iterator(
            nrows, (peak_bnds.0.0..peak_bnds.0.1).map(|i| {
                ecal.0.0 + ecal.0.1 * i as f64
            })
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

                let y0 = ecal.0.0 + ecal.0.1 * y0 as f64;
                let x0 = ecal.1.0 + ecal.1.1 * x0 as f64;

                fit_gauss2d(&x, &y, peak, &[max, x0, y0, 1., 2., -0.78])?
            },
            // BackgroundModel::Erf => {},
            // BackgroundModel::ErfGauss => {},
            // BackgroundModel::ErfLinear1Gauss => {},
            // BackgroundModel::ErfLinear2Gauss => {},
            // BackgroundModel::ErfGauss1Exp => {},
            // BackgroundModel::ErfGauss2Exp => {},
        };

        Ok(())
    }

    fn subtract_bg(&mut self, bg_model: BackgroundModel) {
        if bg_model == BackgroundModel::None {
            return;
        }

        todo!();
    }

    pub fn extract_peak(
        &mut self,
        peak_window_kev: (f64, f64),
        bg_corr: bool,
        bg_model: BackgroundModel
    ) {
        let ecal = self.corrected_ecal.unwrap_or(self.ecal);

        let first_row_e = (M_E_KEV - peak_window_kev.0 / 2. - ecal.0.0) / ecal.0.1;
        let last_row_e = (M_E_KEV + peak_window_kev.0 / 2. - ecal.0.0) / ecal.0.1;

        let first_col_e = (M_E_KEV - peak_window_kev.1 / 2. - ecal.1.0) / ecal.1.1;
        let last_col_e = (M_E_KEV + peak_window_kev.1 / 2. - ecal.1.0) / ecal.1.1;

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

        if bg_corr {
            self.subtract_bg(bg_model);
        }

        let peak_counts = self.peak.as_ref().unwrap().sum();

        self.peak_counts = Some(peak_counts);
        self.dpeak_counts = Some(peak_counts.sqrt());
    }

    pub fn correct_ecal(&mut self, order: EcalCorrectionOrder) -> Result<(), AnalysisError> {
        if order == EcalCorrectionOrder::None {
            return Ok(());
        }

        if self.peak == None {
            self.extract_peak((4., 4.), false, BackgroundModel::None);
        }

        if !self.peak_params.contains_key("x0") || !self.peak_params.contains_key("y0") {
            self.fit_2d_peak(BackgroundModel::None)?;
        }

        let x0 = self.peak_params.get("x0").unwrap().val;
        let y0 = self.peak_params.get("y0").unwrap().val;

        let mut ecal = self.ecal;

        match order {
            EcalCorrectionOrder::Zeroth => {
                ecal.0.0 += M_E_KEV - y0;
                ecal.1.0 += M_E_KEV - x0;
            },
            EcalCorrectionOrder::First => {
                ecal.0.1 *= (M_E_KEV - ecal.0.0) / (y0 - ecal.0.0);
                ecal.1.1 *= (M_E_KEV - ecal.1.0) / (x0 - ecal.1.0);
            },
            EcalCorrectionOrder::None => unreachable!(),
        }

        self.corrected_ecal = Some(ecal);

        Ok(())
    }

    fn integrate(
        &self,
        area: Area,
        in_corrected_peak: bool,
        bin_integration_scheme: BinIntegrationScheme,
        anti_aliasing_scheme: AntiAliasingScheme
    ) {
        todo!();
    }

    pub fn project_digonal(
        &self,
        width: f64,
        bins: &[f64],
    ) -> DBSpectrum {
        todo!();
    }

    pub fn calc_lineshape_param(
        &mut self,
        in_corrected_peak: bool,
        bin_integration_scheme: BinIntegrationScheme,
        anti_aliasing_scheme: AntiAliasingScheme,
    ) {
        todo!();
    }
}
