// NOTE:
// - spectrum counts represent the center of their bins
// - integration happens in channel space, not in energy space

use std::collections::HashMap;

use thiserror::Error;

use crate::constants::M_E_KEV;
use crate::core::utils::{
    spectrum_match, EnergyDetector, LinearCalibration, LossyIntoF64, Spectrum,
    EcalCorrectionOrder, Unit
};
use crate::core::fitting::{self, VarproFitError};
use crate::dbs::fitting::*;

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

    #[error("Lineshape numerator or denominator areas overlap")]
    OverlappingAreas,

    #[error("Error during fitting")]
    FitError {
        #[from]
        inner: VarproFitError,
    },

    #[error("Could not convert units of eres to keV due to missing eres")]
    MissingEnergyResolution,
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

#[derive(Debug, Clone, Copy)]
pub enum Area {
    Peak {
        width: Unit,
    },
    PeakWithOffset {
        width: Unit,
        offset: Unit,
    },
    Spectrum {
        left_bnd: f64,
        right_bnd: f64,
    },
}

impl Area {
    pub fn to_bnds_kev(&self, eres: Option<f64>) -> Option<(f64, f64)> {
        match self {
            Area::Peak { width } => {
                let width_kev = width.to_kev(eres)?;
                Some((M_E_KEV - width_kev / 2., M_E_KEV + width_kev / 2.))
            },
            Area::PeakWithOffset { width, offset } => {
                let width_kev = width.to_kev(eres)?;
                let offset_kev = offset.to_kev(eres)?;
                Some((
                    M_E_KEV - width_kev / 2. + offset_kev,
                    M_E_KEV + width_kev / 2. + offset_kev
                ))
            },
            Area::Spectrum { left_bnd, right_bnd } => {
                Some((*left_bnd, *right_bnd))
            },
        }
    }
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

pub const STD_LINESHAPE_PARAMS: &[LineshapeParamDefinition; 4] = &[
    LineshapeParamDefinition {
        name: "S",
        num: &[Area::Peak { width: Unit::KeV(1.) }],
        denom: &[Area::Peak { width: Unit::KeV(20.) }],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "W",
        num: &[
            Area::PeakWithOffset { width: Unit::KeV(1.), offset: Unit::KeV(3.) },
            Area::PeakWithOffset { width: Unit::KeV(1.), offset: Unit::KeV(-3.) },
        ],
        denom: &[Area::Peak { width: Unit::KeV(20.) }],
        in_corrected_peak: true,
    },
    LineshapeParamDefinition {
        name: "V/P",
        num: &[Area::Spectrum { left_bnd: 400., right_bnd: 500. }],
        denom: &[Area::Spectrum { left_bnd: 501., right_bnd: 521. }],
        in_corrected_peak: false,
    },
    LineshapeParamDefinition {
        name: "P/T",
        num: &[Area::Spectrum { left_bnd: 501., right_bnd: 521. }],
        denom: &[Area::Spectrum {
            left_bnd: f64::NEG_INFINITY, right_bnd: f64::INFINITY
        }],
        in_corrected_peak: false,
    },
];

pub struct DBSpectrum {
    pub spectrum: Spectrum,
    pub detector: EnergyDetector,
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
        detector: EnergyDetector,
    ) -> Self {
        let counts: u64 = spectrum_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum()
        );

        let dcounts = (counts as f64).sqrt();

        DBSpectrum {
            spectrum,
            detector,
            counts,
            dcounts,
            peak_counts: f64::NAN,
            dpeak_counts: f64::NAN,
            peak: None,
            peak_bnds: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn default_analyze(&mut self) -> Result<(), AnalysisError> {
        self.correct_ecal(EcalCorrectionOrder::First)?;

        self.extract_peak(30., true, PeakModel::ErfLinear1Gauss)?;

        for param in STD_LINESHAPE_PARAMS {
            self.calc_lineshape_param(param, BinIntegrationScheme::Const)?;
        }

        self.peak = None;
        self.peak_bnds = None;

        Ok(())
    }

    fn get_peak_energies(&self) -> Vec<f64> {
        let ecal = self.detector.corrected_ecal.unwrap_or(self.detector.ecal);

        let bnds = self.peak_bnds.unwrap();

        (bnds.0..=bnds.1).map(|i| ecal.from_index(i)).collect::<Vec<f64>>()
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
        let ecal = self.detector.corrected_ecal.unwrap_or(self.detector.ecal);

        let left_peak_idx = ecal.to_index_rounded(M_E_KEV - peak_width / 2.);

        let right_peak_idx = ecal.to_index_rounded(M_E_KEV + peak_width / 2.);

        self.peak_bnds = Some((left_peak_idx, right_peak_idx));

        self.peak = Some(spectrum_match!(
            &self.spectrum,
            arr => arr.rows(left_peak_idx, right_peak_idx - left_peak_idx + 1)
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
    ) -> Result<LinearCalibration, AnalysisError> {
        if order == EcalCorrectionOrder::None {
            return Ok(self.detector.ecal);
        }

        if self.peak == None {
            self.extract_peak(20., false, PeakModel::Gauss)?;
        }

        self.fit_peak(PeakModel::Gauss)?;

        let c = self.peak_params.get("x0_1").unwrap().val;

        let mut ecal = self.detector.ecal;

        match order {
            EcalCorrectionOrder::Zeroth => {
                ecal.offset += M_E_KEV - c;
            },
            EcalCorrectionOrder::First => {
                ecal.scale *= (M_E_KEV - ecal.offset) / (c - ecal.offset);
            },
            EcalCorrectionOrder::None => unreachable!(),
        }

        self.detector.corrected_ecal = Some(ecal);

        Ok(ecal)
    }

    fn integrate_generic<T: LossyIntoF64>(
        spectrum: &[T],
        mut lower_ch_bnd: f64,
        mut upper_ch_bnd: f64,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        // TODO: We could allow for a lower_ch_bdn of >= -0.5
        // but I'm too lazy to deal with this special edge case now
        if lower_ch_bnd.is_infinite() {
            lower_ch_bnd = 0.;
        }
        if upper_ch_bnd.is_infinite() {
            upper_ch_bnd = spectrum.len() as f64 - 1.;
        }

        if lower_ch_bnd < 0. || lower_ch_bnd > spectrum.len() as f64 - 1. {
            return Err(AnalysisError::OutOfBounds);
        }
        if upper_ch_bnd < 0. || upper_ch_bnd > spectrum.len() as f64 - 1. {
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

    pub fn integrate(
        &self,
        area: Area,
        in_corrected_peak: bool,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<f64, AnalysisError> {
        let (lower_bnd, upper_bnd) = area.to_bnds_kev(self.detector.eres)
            .ok_or(AnalysisError::MissingEnergyResolution)?;

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

        let ecal = self.detector.corrected_ecal.unwrap_or(self.detector.ecal);

        let mut lower_ch_bnd = ecal.to_index_f64(lower_bnd);
        let mut upper_ch_bnd = ecal.to_index_f64(upper_bnd);

        if in_corrected_peak {
            let peak_bnds = self.peak_bnds.unwrap();
            lower_ch_bnd -= peak_bnds.0 as f64;
            upper_ch_bnd -= peak_bnds.0 as f64;
            let peak = self.peak.as_ref().unwrap();
            DBSpectrum::integrate_generic(
                peak, lower_ch_bnd, upper_ch_bnd, bin_integration_scheme,
            )
        } else {
            spectrum_match!(
                &self.spectrum, spectrum => DBSpectrum::integrate_generic(
                    spectrum.as_slice(),
                    lower_ch_bnd,
                    upper_ch_bnd,
                    bin_integration_scheme,
                )
            )
        }
    }

    fn get_area_bnds(&self, areas: &[Area]) -> Result<Vec<(f64, f64)>, AnalysisError> {
        let bnds = areas
            .iter()
            .map(|area| area.to_bnds_kev(self.detector.eres))
            .collect::<Option<Vec<_>>>()
            .ok_or(AnalysisError::MissingEnergyResolution)?;

        for i in 0..areas.len() - 1 {
            for j in i + 1..areas.len() {
                if !(bnds[i].1 < bnds[j].0 || bnds[j].1 < bnds[i].0) {
                    return Err(AnalysisError::OverlappingAreas);
                }
            }
        }

        Ok(bnds)
    }

    pub fn calc_lineshape_param(
        &mut self,
        definition: &LineshapeParamDefinition,
        bin_integration_scheme: BinIntegrationScheme,
    ) -> Result<&LineshapeParam, AnalysisError> {
        let num_bnds = self.get_area_bnds(definition.num)?;
        let denom_bnds = self.get_area_bnds(definition.denom)?;

        let num = definition.num
            .iter()
            .map(|&area| self.integrate(
                area, definition.in_corrected_peak, bin_integration_scheme
            ))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum::<f64>();

        let denom = definition.denom
            .iter()
            .map(|&area| self.integrate(
                area, definition.in_corrected_peak, bin_integration_scheme
            ))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum::<f64>();

        let mut common = 0.;

        for num_bnd in num_bnds.iter() {
            for denom_bnd in denom_bnds.iter() {
                let left_bnd = num_bnd.0.max(denom_bnd.0);
                let right_bnd = num_bnd.1.min(denom_bnd.1);
                if left_bnd < right_bnd {
                    common += self.integrate(
                        Area::Spectrum { left_bnd, right_bnd },
                        definition.in_corrected_peak,
                        bin_integration_scheme,
                    )?;
                }
            }
        }

        let val = num / denom;
        let err = val * (
            (num - common) / num.powi(2) +
            (denom - common) / denom.powi(2) +
            common * ((denom - num) / (denom * num)).powi(2)
        ).sqrt();

        let ls_param = LineshapeParam {
            val,
            err,
            num: definition.num.to_owned(),
            denom: definition.denom.to_owned(),
            in_corrected_peak: definition.in_corrected_peak,
        };

        self.lineshape_params.insert(definition.name.into(), ls_param);

        Ok(self.lineshape_params.get(definition.name).unwrap())
    }
}
