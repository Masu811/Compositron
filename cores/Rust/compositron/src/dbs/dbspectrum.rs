use std::collections::HashMap;

use thiserror::Error;

use crate::constants::M_E_KEV;
use crate::core::utils::{spectrum_match, LossyIntoF64, Spectrum};
use crate::core::fitting;
use crate::dbs::fitting::*;

// NOTE:
// - spectrum counts represent the center of their bins
// - integration happens in channel space, not in energy space

#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("Input parameter violates constraint {violated_constraint}")]
    InputError {
        violated_constraint: String,
    },

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

    #[error(
        "Attempting to integrate areas on both sides of the peak which are \
        overlapping (area_width > area_dist)"
    )]
    OverlappingAreas,

    #[error("Error during fitting")]
    FitError {
        #[from]
        inner: fitting::FitError,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EcalCorrectionOrder {
    Zeroth,
    First,
    None
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PeakModel {
    Gauss,
    ErfLinear1Gauss,
    ErfLinear2Gauss,
    ErfLinear3Gauss,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BinIntegrationScheme {
    Const,
    Linear,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PeakLineshapeParamEvalSide {
    Both,
    Left,
    Right,
}

#[derive(Debug)]
pub enum LineshapeParam {
    PeakSym {
        val: f64,
        err: f64,
        num_width: f64,
        denom_width: f64,
    },
    PeakAsym {
        val: f64,
        err: f64,
        num_width: f64,
        denom_width: f64,
        num_dist: f64,
        denom_dist: f64,
        side: PeakLineshapeParamEvalSide,
    },
    Spectrum {
        val: f64,
        err: f64,
        num_bnds: (f64, f64),
        denom_bnds: (f64, f64),
    },
}

#[derive(Debug)]
pub enum LineshapeParamDefinition<'a> {
    PeakSym {
        name: &'a str,
        num_width: f64,
        denom_width: f64,
    },
    PeakAsym {
        name: &'a str,
        num_width: f64,
        denom_width: f64,
        num_dist: f64,
        denom_dist: f64,
        num_side: PeakLineshapeParamEvalSide,
        denom_side: PeakLineshapeParamEvalSide,
    },
    Spectrum {
        name: &'a str,
        num_bnds: (f64, f64),
        denom_bnds: (f64, f64),
    },
}

pub const STD_LINESHAPE_PARAMS: &[LineshapeParamDefinition; 4] = &[
    LineshapeParamDefinition::PeakSym {
        name: "S",
        num_width: 1.,
        denom_width: 30.,
    },
    LineshapeParamDefinition::PeakAsym {
        name: "W",
        num_width: 1.,
        denom_width: 30.,
        num_dist: 3.,
        denom_dist: 0.,
        num_side: PeakLineshapeParamEvalSide::Both,
        denom_side: PeakLineshapeParamEvalSide::Both,
    },
    LineshapeParamDefinition::Spectrum {
        name: "V/P",
        num_bnds: (510., 512.),
        denom_bnds: (481., 541.),
    },
    LineshapeParamDefinition::Spectrum {
        name: "P/T",
        num_bnds: (510., 512.),
        denom_bnds: (f64::NEG_INFINITY, f64::INFINITY),
    },
];

enum InputParam {
    NumWidth(f64),
    DenomWidth(f64),
    NumDist(f64),
    DenomDist(f64),
    NumBounds((f64, f64)),
    DenomBounds((f64, f64)),
}

pub struct DBSpectrum {
    pub spectrum: Spectrum,
    pub detname: String,
    pub ecal: (f64, f64),
    pub corrected_ecal: Option<(f64, f64)>,
    pub eres: Option<f64>,
    pub counts: u64,
    pub dcounts: f64,
    pub peak: Option<Vec<f64>>,
    pub peak_bnds: Option<(usize, usize)>,
    pub peak_counts: f64,
    pub dpeak_counts: f64,
    pub peak_params: HashMap<&'static str, fitting::SimpleFitParam>,
    pub lineshape_params: HashMap<String, LineshapeParam>,
}

impl DBSpectrum {
    pub fn new(
        spectrum: Spectrum,
        detname: String,
        ecal: (f64, f64),
        eres: Option<f64>,
    ) -> Self {
        let counts: u64 = spectrum_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum()
        );

        let dcounts = (counts as f64).sqrt();

        DBSpectrum {
            spectrum: spectrum,
            detname: detname,
            ecal: ecal,
            corrected_ecal: None,
            eres: eres,
            counts: counts,
            dcounts: dcounts,
            peak_counts: f64::NAN,
            dpeak_counts: f64::NAN,
            peak: None,
            peak_bnds: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn analyze(
        &mut self,
        ecal_corr_order: EcalCorrectionOrder,
        peak_width: f64,
        subtract_bg: bool,
        peak_model: PeakModel,
        param_definitions: &[LineshapeParamDefinition],
        bin_integration_scheme: BinIntegrationScheme,
        save_extracted_peak: bool,
    ) -> Result<(), AnalysisError> {
        self.correct_ecal(ecal_corr_order)?;

        self.extract_peak(peak_width, subtract_bg, peak_model)?;

        for param in param_definitions {
            self.calc_lineshape_param(param, bin_integration_scheme)?;
        }

        if !save_extracted_peak {
            self.peak = None;
            self.peak_bnds = None;
        }

        Ok(())
    }

    fn get_peak_energies(&self) -> Vec<f64> {
        let ecal = self.corrected_ecal.unwrap_or(self.ecal);
        let bnds = self.peak_bnds.unwrap();

        (bnds.0..=bnds.1)
            .map(|i| ecal.1 * i as f64 + ecal.0)
            .collect::<Vec<f64>>()
    }

    fn fit_peak(
        &mut self, peak_model: PeakModel,
    ) -> Result<(), AnalysisError> {
        let x = self.get_peak_energies();
        let y = self.peak.as_ref().unwrap();

        match peak_model {
            PeakModel::Gauss => {
                self.peak_params = fit_gaussian(&x, &y, vec![511., 1.])?;
            },
            PeakModel::ErfLinear1Gauss => {
                self.peak_params = fit_erf_linear_1_gauss(&x, &y)?;
            },
            PeakModel::ErfLinear2Gauss => {
                self.peak_params = fit_erf_linear_2_gauss(&x, &y)?;
            },
            PeakModel::ErfLinear3Gauss => {
                self.peak_params = fit_erf_linear_3_gauss(&x, &y)?;
            },
        }

        Ok(())
    }

    fn subtract_bg(
        &mut self, peak_model: PeakModel
    ) -> Result<(), AnalysisError> {
        self.fit_peak(peak_model)?;

        let x = self.get_peak_energies();
        let y = self.peak.as_mut().unwrap();

        let corr = background(
            &x,
            self.peak_params.get("erf_amp").unwrap().val,
            self.peak_params.get("x0_1").unwrap().val,
            self.peak_params.get("sig_1").unwrap().val,
            self.peak_params.get("lin").unwrap().val,
            self.peak_params.get("const").unwrap().val,
        );

        for (s, bg) in y.iter_mut().zip(corr) {
            *s = *s - bg;
        }

        Ok(())
    }

    pub fn extract_peak(
        &mut self, peak_width: f64, bg_corr: bool, peak_model: PeakModel,
    ) -> Result<&[f64], AnalysisError> {
        let ecal = self.corrected_ecal.unwrap_or(self.ecal);

        let left_peak_idx = (
            (M_E_KEV - peak_width / 2. - ecal.0) / ecal.1
        ).round() as usize;

        let right_peak_idx = (
            (M_E_KEV + peak_width / 2. - ecal.0) / ecal.1
        ).round() as usize;

        self.peak_bnds = Some((left_peak_idx, right_peak_idx));

        self.peak = Some(spectrum_match!(
            &self.spectrum,
            arr => arr.rows(left_peak_idx, right_peak_idx - left_peak_idx)
                .iter()
                .map(|&x| x as f64)
                .collect()
        ));

        if bg_corr == true {
            self.subtract_bg(peak_model)?;
        }

        Ok(self.peak.as_ref().unwrap())
    }

    pub fn correct_ecal(
        &mut self, order: EcalCorrectionOrder
    ) -> Result<(f64, f64), AnalysisError> {
        if order == EcalCorrectionOrder::None {
            return Ok(self.ecal);
        }

        if self.peak == None {
            self.extract_peak(20., false, PeakModel::Gauss)?;
        }

        self.fit_peak(PeakModel::Gauss)?;

        let c = self.peak_params.get("x0_1").unwrap().val;

        match order {
            EcalCorrectionOrder::Zeroth => {
                self.ecal.0 += M_E_KEV - c;
            },
            EcalCorrectionOrder::First => {
                self.ecal.1 *= (M_E_KEV - self.ecal.0) / (c - self.ecal.0);
            },
            EcalCorrectionOrder::None => unreachable!(),
        }

        Ok(self.ecal)
    }

    fn integrate<T: LossyIntoF64>(
        spectrum: &[T],
        lower_ch_bnd: f64,
        upper_ch_bnd: f64,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        // TODO: We could allow for a lower_ch_bdn of >= -0.5
        // but I'm too lazy to deal with this special edge case now
        if lower_ch_bnd <= 0. && !lower_ch_bnd.is_infinite() {
            return Err(AnalysisError::OutOfBounds);
        }
        if upper_ch_bnd >= spectrum.len() as f64 - 1. && !upper_ch_bnd.is_infinite() {
            return Err(AnalysisError::OutOfBounds);
        }

        let mut counts: f64 = 0.;

        match bin_integration_scheme {
            BinIntegrationScheme::Const => {
                // index of lowest channel still inside ROI
                let lowest_ch = lower_ch_bnd.round() as usize;
                // index of highest channel still inside ROI
                let highest_ch = upper_ch_bnd.round() as usize;

                counts += spectrum[lowest_ch..=highest_ch]
                    .iter()
                    .map(|&x| T::lossy_into_f64(x))
                    .sum::<f64>();

                // Left edge counts

                let lowest_ch_counts = T::lossy_into_f64(spectrum[lowest_ch]);
                let x1 = lowest_ch as f64 - 0.5;
                counts -= lowest_ch_counts * (lower_ch_bnd - x1);

                // Right edge counts

                let highest_ch_counts = T::lossy_into_f64(spectrum[highest_ch]);
                let x2 = highest_ch as f64 + 0.5;
                counts -= highest_ch_counts * (x2 - upper_ch_bnd);
            },
            BinIntegrationScheme::Linear => {
                // index of lowest channel still inside ROI
                let lowest_ch = lower_ch_bnd.floor() as usize;
                // index of highest channel still inside ROI
                let highest_ch = upper_ch_bnd.ceil() as usize;

                counts += spectrum[lowest_ch+1..highest_ch]
                    .iter()
                    .map(|&x| T::lossy_into_f64(x))
                    .sum::<f64>();

                counts += 0.5 * T::lossy_into_f64(spectrum[lowest_ch]);
                counts += 0.5 * T::lossy_into_f64(spectrum[highest_ch]);

                // Left edge counts

                let x1 = lowest_ch as f64;

                let y1 = T::lossy_into_f64(spectrum[lowest_ch]);
                let y2 = T::lossy_into_f64(spectrum[lowest_ch + 1]);

                let y_l = (y2 - y1) * (lower_ch_bnd - x1) + y1;

                counts -= 0.5 * (y1 + y_l) * (lower_ch_bnd - x1);

                // Right edge counts

                let x1 = (highest_ch - 1) as f64;
                let x2 = highest_ch as f64;

                let y1 = T::lossy_into_f64(spectrum[highest_ch - 1]);
                let y2 = T::lossy_into_f64(spectrum[highest_ch]);

                let y_u = (y2 - y1) * (upper_ch_bnd - x1) + y1;

                counts -= 0.5 * (y_u + y2) * (x2 - upper_ch_bnd);
            },
        }

        Ok(counts)
    }

    pub fn integrate_spectrum(
        &self,
        lower_bnd: f64,
        upper_bnd: f64,
        in_corrected_peak: bool,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        if in_corrected_peak && self.peak == None {
            return Err(AnalysisError::NoPeakExtracted);
        }
        if lower_bnd >= upper_bnd {
            return Err(AnalysisError::InputError {
                violated_constraint: format!(
                    "Lower bound < upper bound (got {lower_bnd} < {upper_bnd})"
                )
            });
        }

        let ecal = self.corrected_ecal.unwrap_or(self.ecal);

        let mut lower_ch_bnd = (lower_bnd - ecal.0) / ecal.1;
        let mut upper_ch_bnd = (upper_bnd - ecal.0) / ecal.1;

        if in_corrected_peak {
            let peak_bnds = self.peak_bnds.unwrap();
            lower_ch_bnd -= peak_bnds.0 as f64;
            upper_ch_bnd -= peak_bnds.0 as f64;
            let peak = self.peak.as_ref().unwrap();
            DBSpectrum::integrate(
                peak, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme,
            )
        } else {
            spectrum_match!(
                &self.spectrum, spectrum => DBSpectrum::integrate(
                    spectrum.as_slice(),
                    lower_ch_bnd,
                    upper_ch_bnd,
                    bin_integration_scheme,
                )
            )
        }
    }

    fn check_input(input: InputParam) -> Result<(), AnalysisError> {
        let violation = match input {
            InputParam::NumWidth(x) => x <= 0.,
            InputParam::DenomWidth(x) => x <= 0.,
            InputParam::NumDist(x) => x < 0.,
            InputParam::DenomDist(x) => x < 0.,
            InputParam::NumBounds((x, y)) => y <= x,
            InputParam::DenomBounds((x, y)) => y <= x,
        };

        if violation {
            let msg = match input {
                InputParam::NumWidth(x) => {
                    format!("num_width > 0 (got {x})")
                },
                InputParam::DenomWidth(x) => {
                    format!("denom_width > 0 (got {x})")
                },
                InputParam::NumDist(x) => {
                    format!("num_dist >= 0 (got {x})")
                },
                InputParam::DenomDist(x) => {
                    format!("denom_dist >= 0 (got {x})")
                },
                InputParam::NumBounds((x, y)) => {
                    format!("num_bnds.0 < num_bnds.1 (got {x} < {y})")
                },
                InputParam::DenomBounds((x, y)) => {
                    format!("denom_bnds.0 < denom_bnds.1 (got {x} < {y})")
                },
            };

            Err(AnalysisError::InputError {
                violated_constraint: msg.into()
            })
        } else {
            Ok(())
        }
    }

    fn err_divide(a_n: f64, a_d: f64, a_c: f64) -> (f64, f64) {
        let val = (a_n + a_c) / (a_d + a_c);

        let err = val * (
            a_n / (a_n + a_c).powi(2) +
            a_d / (a_d + a_c).powi(2) +
            a_c * ((a_d - a_n) / ((a_d + a_c) * (a_n + a_c))).powi(2)
        ).sqrt();

        (val, err)
    }

    fn get_area_bnds(
        width: f64,
        dist: f64,
        side: PeakLineshapeParamEvalSide,
    ) -> Vec<f64> {
        let half_width = width / 2.;
        match side {
            PeakLineshapeParamEvalSide::Left => {
                let mid = M_E_KEV - dist - half_width;
                vec![mid - half_width, mid + half_width]
            },
            PeakLineshapeParamEvalSide::Right => {
                let mid = M_E_KEV + dist + half_width;
                vec![mid - half_width, mid + half_width]
            },
            PeakLineshapeParamEvalSide::Both => {
                if dist == 0. {
                    let mid = M_E_KEV;
                    vec![mid - width, mid + width]
                } else {
                    let mid_l = M_E_KEV - dist - half_width;
                    let mid_r = M_E_KEV + dist + half_width;
                    vec![
                        mid_l - half_width, mid_l + half_width,
                        mid_r - half_width, mid_r + half_width,
                    ]
                }
            },
        }
    }

    fn get_overlap(
        &self,
        bnds_a: (f64, f64),
        bnds_b: (f64, f64),
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        let common_bnds = [bnds_a.0.max(bnds_b.0), bnds_a.1.min(bnds_b.1)];

        if common_bnds[0] < common_bnds[1] {
            Ok(self.integrate_multiple_areas(
                &common_bnds, bin_integration_scheme
            )?)
        } else {
            Ok(0.)
        }
    }

    fn integrate_multiple_areas(
        &self,
        bnds: &[f64],
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        let mut area = self.integrate_spectrum(
            bnds[0], bnds[1], true, bin_integration_scheme
        )?;

        if bnds.len() == 4 {
            area += self.integrate_spectrum(
                bnds[2], bnds[3], true, bin_integration_scheme
            )?;
        }

        Ok(area)
    }

    fn calc_peak_param(
        &mut self,
        name: &str,
        num_width: f64,
        denom_width: f64,
        num_dist: f64,
        denom_dist: f64,
        num_side: PeakLineshapeParamEvalSide,
        denom_side: PeakLineshapeParamEvalSide,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<&LineshapeParam, AnalysisError> {
        if self.peak == None {
            return Err(AnalysisError::NoPeakExtracted);
        }

        DBSpectrum::check_input(InputParam::NumWidth(num_width))?;
        DBSpectrum::check_input(InputParam::DenomWidth(denom_width))?;
        DBSpectrum::check_input(InputParam::NumDist(num_dist))?;
        DBSpectrum::check_input(InputParam::DenomDist(denom_dist))?;

        let num_bnds = DBSpectrum::get_area_bnds(
            num_width, num_dist, num_side
        );
        let denom_bnds = DBSpectrum::get_area_bnds(
            denom_width, denom_dist, denom_side
        );

        let mut num = self.integrate_multiple_areas(
            &num_bnds, bin_integration_scheme
        )?;
        let mut denom = self.integrate_multiple_areas(
            &denom_bnds, bin_integration_scheme
        )?;

        let common: f64 = match (num_bnds.len(), denom_bnds.len()) {
            (2, 2) => {
                self.get_overlap(
                    (num_bnds[0], num_bnds[1]),
                    (denom_bnds[0], denom_bnds[1]),
                    bin_integration_scheme,
                )?
            },
            (4, 2) => {
                self.get_overlap(
                    (num_bnds[0], num_bnds[1]),
                    (denom_bnds[0], denom_bnds[1]),
                    bin_integration_scheme,
                )? + self.get_overlap(
                    (num_bnds[2], num_bnds[3]),
                    (denom_bnds[0], denom_bnds[1]),
                    bin_integration_scheme,
                )?
            },
            (2, 4) => {
                self.get_overlap(
                    (num_bnds[0], num_bnds[1]),
                    (denom_bnds[0], denom_bnds[1]),
                    bin_integration_scheme,
                )? + self.get_overlap(
                    (num_bnds[0], num_bnds[1]),
                    (denom_bnds[2], denom_bnds[3]),
                    bin_integration_scheme,
                )?
            },
            (4, 4) => {
                2. * self.get_overlap(
                    (num_bnds[0], num_bnds[1]),
                    (denom_bnds[0], denom_bnds[1]),
                    bin_integration_scheme,
                )?
            },
            _ => unreachable!(),
        };

        num -= common;
        denom -= common;

        let (val, err) = DBSpectrum::err_divide(num, denom, common);

        let p = LineshapeParam::PeakSym {val, err, num_width, denom_width};

        self.lineshape_params.insert(name.into(), p);

        Ok(self.lineshape_params.get(name).unwrap())
    }

    fn calc_spectrum_param(
        &mut self,
        name: &str,
        num_bnds: (f64, f64),
        denom_bnds: (f64, f64),
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<&LineshapeParam, AnalysisError> {
        DBSpectrum::check_input(InputParam::NumBounds(num_bnds))?;
        DBSpectrum::check_input(InputParam::DenomBounds(denom_bnds))?;

        let mut num = self.integrate_spectrum(
            num_bnds.0, num_bnds.1, false, bin_integration_scheme
        )?;
        let mut denom = self.integrate_spectrum(
            num_bnds.0, num_bnds.1, false, bin_integration_scheme
        )?;

        let common = self.get_overlap(
            num_bnds, denom_bnds, bin_integration_scheme,
        )?;

        num -= common;
        denom -= common;

        let (val, err) = DBSpectrum::err_divide(num, denom, common);

        let p = LineshapeParam::Spectrum { val, err, num_bnds, denom_bnds };

        self.lineshape_params.insert(name.into(), p);

        Ok(self.lineshape_params.get(name).unwrap())
    }

    pub fn calc_lineshape_param(
        &mut self,
        definition: &LineshapeParamDefinition,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<&LineshapeParam, AnalysisError> {
        match definition {
            &LineshapeParamDefinition::PeakSym {
                name, num_width, denom_width
            } => {
                self.calc_peak_param(
                    name,
                    num_width,
                    denom_width,
                    0.,
                    0.,
                    PeakLineshapeParamEvalSide::Both,
                    PeakLineshapeParamEvalSide::Both,
                    bin_integration_scheme,
                )
            },
            &LineshapeParamDefinition::PeakAsym {
                name, num_width, denom_width, num_dist, denom_dist, num_side, denom_side
            } => {
                self.calc_peak_param(
                    name,
                    num_width,
                    denom_width,
                    num_dist,
                    denom_dist,
                    num_side,
                    denom_side,
                    bin_integration_scheme,
                )
            },
            &LineshapeParamDefinition::Spectrum {
                name, num_bnds, denom_bnds
            } => {
                self.calc_spectrum_param(
                    name, num_bnds, denom_bnds, bin_integration_scheme
                )
            },
        }
    }
}
