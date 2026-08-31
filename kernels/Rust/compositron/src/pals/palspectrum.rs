use std::collections::HashMap;

use thiserror::Error;

use crate::constants::FWHM_OVER_SIGMA;
use crate::core::utils::{spectrum_match, TimingDetectorPair, Spectrum};
use crate::core::fitting::{FitParam, LMFitError, SimpleFitParam, VarproFitError};
use crate::dbs::fitting::fit_gauss;
use crate::pals::fitting::{fit_lifetime_spectrum, FitResult};
use crate::pals::model::*;


#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("Fit failed")]
    LMFitError {
        #[from]
        source: LMFitError
    },

    #[error("Fit failed")]
    VarproFitError {
        #[from]
        source: VarproFitError
    },
}


#[derive(Debug, Error)]
pub enum ReportError {
    #[error("No fit has been performed yet")]
    NoFitPerformed,

    #[error("I/O Error")]
    IOError {
        #[from]
        source: std::fmt::Error,
    }
}


pub struct PALSpectrum {
    pub spectrum: Spectrum,
    pub detpair: TimingDetectorPair,
    pub counts: u64,
    pub peak_params: HashMap<&'static str, SimpleFitParam>,
    pub fit_result: Option<FitResult>,
}


impl PALSpectrum {
    pub fn new(
        spectrum: Spectrum,
        detpair: TimingDetectorPair,
    ) -> Self {
        let counts = spectrum_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );

        PALSpectrum {
            spectrum,
            detpair,
            counts,
            peak_params: HashMap::new(),
            fit_result: None,
        }
    }


    fn get_peak_center(&mut self) -> Result<f64, AnalysisError> {
        let (argmax, _) = spectrum_match!(
            &self.spectrum, arr => arr
                .iter()
                .enumerate()
                .map(|(i, &y)| (i, y as u64))
                .max_by_key(|&(_, y)| y)
                .unwrap_or((0, 0))
        );

        let peak_center_window = 15;

        let roi = argmax - peak_center_window .. argmax + peak_center_window;

        let y = spectrum_match!(
            &self.spectrum, arr => arr.rows(roi.start, roi.len())
                .iter()
                .map(|&i| i as f64)
                .collect::<Vec<f64>>()
        );

        let x = roi.map(|i| i as f64).collect::<Vec<f64>>();

        let mut peak_params = fit_gauss(&x, &y, vec![argmax as f64, 100.])?;

        let peak_height = peak_params.remove("amp_1").unwrap();
        let peak_center = peak_params.remove("x0_1").unwrap();
        let peak_width = peak_params.remove("sig_1").unwrap();

        let x0 = peak_center.val;

        self.peak_params.insert("center_idx", peak_center);
        self.peak_params.insert("height", peak_height);
        self.peak_params.insert("fwhm", SimpleFitParam {
            val: peak_width.val * FWHM_OVER_SIGMA,
            err: peak_width.err * FWHM_OVER_SIGMA
        });

        Ok(x0)
    }


    pub fn fit(
        &mut self,
        model: LifetimeModel,
        fit_start_idx: usize,
        fit_end_idx: usize,
    ) -> Result<(), AnalysisError> {
        let peak_center = self.peak_params
            .get("center_idx")
            .map(|p| p.val)
            .unwrap_or(self.get_peak_center()?);

        let roi = fit_start_idx..fit_end_idx;

        let y = spectrum_match!(
            &self.spectrum, arr => arr.rows(roi.start, roi.len())
                .iter()
                .map(|&i| i as f64)
                .collect::<Vec<f64>>()
        );

        let tcal = self.detpair.corrected_tcal.unwrap_or(self.detpair.tcal);

        let x = roi.map(|i| tcal.from_index(i)).collect::<Vec<f64>>();

        let fit_result = fit_lifetime_spectrum(
            &x,
            &y,
            model,
            self.counts as f64,
            tcal.from_index_f64(peak_center),
        )?;

        self.fit_result = Some(fit_result);

        Ok(())
    }


    fn print_params(params: &[&FitParam], names: Vec<String>) -> String {
        let mut text = String::new();
        let n_rows = names.len();

        let values = params
            .iter()
            .zip(&names)
            .map(|(param, name)| {
                if name.starts_with("int") {
                    format!("{:.2}", 100. * param.val)
                } else if param.val.abs() >= 1e4 {
                    format!("{:.3e}", param.val)
                } else {
                    format!("{:.2}", param.val)
                }
            })
            .collect::<Vec<String>>();

        let abs_errors = params
            .iter()
            .zip(&names)
            .map(|(param, name)| {
                if name.starts_with("int") {
                    format!("{:.1}", 100. * param.err)
                } else if param.err.abs() >= 1e4 {
                    format!("{:.3e}", param.err)
                } else {
                    format!("{:.2}", param.err)
                }
            })
            .collect::<Vec<String>>();

        let rel_errors = params
            .iter()
            .map(|param| format!("{:.2} %", (param.err / param.val * 100.).abs()))
            .collect::<Vec<String>>();

        let mins = params
            .iter()
            .zip(&names)
            .map(|(param, name)| {
                if name.starts_with("int") {
                    format!("{:.1}", 100. * param.min)
                } else {
                    format!("{:.1}", param.min)
                }
            })
            .collect::<Vec<String>>();


        let maxs = params
            .iter()
            .zip(&names)
            .map(|(param, name)| {
                if name.starts_with("int") {
                    format!("{:.1}", 100. * param.max)
                } else {
                    format!("{:.1}", param.max)
                }
            })
            .collect::<Vec<String>>();

        let fixed = params
            .iter()
            .map(|param| {
                if param.vary {
                    "".to_string()
                } else {
                    "Fixed".to_string()
                }
            })
            .collect::<Vec<String>>();

        let at_bnd = params
            .iter()
            .map(|param| {
                if (param.val - param.min).abs() < 1e-3 || (param.val - param.max).abs() < 1e-3 {
                    "At Boundary".to_string()
                } else {
                    "".to_string()
                }
            })
            .collect::<Vec<String>>();

        let columns = [
            names, values, abs_errors, rel_errors, mins, maxs, fixed, at_bnd
        ];

        let col_widths = columns
            .iter()
            .map(|arr| arr.iter().map(|x| x.len()).max().unwrap_or(0))
            .collect::<Vec<usize>>();

        let headers = [
            "Parameter", "Value", "Stderr (abs)", "Stderr (rel)", "Min", "Max",
            "Fixed", ""
        ];

        let head_widths = headers
            .iter()
            .map(|head| head.len())
            .collect::<Vec<usize>>();

        let widths = col_widths
            .iter()
            .zip(&head_widths)
            .map(|(c, h)| *c.max(h))
            .collect::<Vec<usize>>();

        let formats = ["<", ">", ">", ">", ">", ">", "<", "<"];

        for (i, head) in headers.iter().enumerate() {
            text.push_str(&format!(" {:^1$} ", head, widths[i]));
        }

        text.push('\n');

        for wid in &widths {
            text.push_str(&format!(" {:-^1$} ", "", wid));
        }

        text.push('\n');

        for row in 0..n_rows {
            for (col, column) in columns.iter().enumerate() {
                if formats[col] == "<" {
                    text.push_str(&format!(" {:<1$} ", column[row], widths[col]));
                } else {
                    text.push_str(&format!(" {:>1$} ", column[row], widths[col]));
                }
            }
            text.push('\n');
        }

        text
    }


    pub fn fit_report(
        &self,
    ) -> Result<String, ReportError> {
        let mut text = String::new();

        let Some(result) = &self.fit_result else {
            return Err(ReportError::NoFitPerformed);
        };

        let model = &result.model;

        let n_l = model.lifetime_components.len();
        let n_r = model.resolution_components.len();

        let thick_hline = format!("{:#^100}\n", "");
        let thin_hline = format!("{:-^100}\n", "");

        text.push_str(&thick_hline);
        text.push_str(&format!("{:^100}\n", "Fit report"));
        text.push_str(&thick_hline);
        text.push('\n');
        text.push_str(&format!("Statistics: {}\n", self.counts));
        text.push_str(&format!("Number of Data Points: {}\n", result.n_dpoints));
        text.push('\n');
        text.push_str(&format!("Lifetime Components:   {}\n", n_l));
        text.push_str(&format!("Resolution Components: {}\n", n_r));
        if model.background_component.default_initialized {
            text.push_str("Background Components: 0\n");
        } else {
            text.push_str("Background Components: 1\n");
        }
        text.push('\n');
        text.push_str("Fit Method: Levenberg-Marquardt\n");
        text.push_str("Fit Status: Success\n");
        text.push_str(&format!("Reason for Termination: {}\n", result.fit_status));
        text.push_str(&format!("Number of Function Evaluations: {}\n", result.n_eval));
        if let Some(red_chi_2) = result.red_chi_2 {
            text.push_str(&format!("Reduced Chi squared: {}\n", red_chi_2));
        } else {
            text.push_str("Reduced Chi squared: NaN\n");
        }
        text.push('\n');
        text.push_str(&thin_hline);
        text.push_str(&format!("{:^100}\n", "Parameters"));
        text.push_str(&thin_hline);
        text.push('\n');
        text.push_str(&format!("{:=^100}\n", " General Parameters "));
        text.push('\n');

        text.push_str(&PALSpectrum::print_params(
            &[&model.scale_component, &model.shift_component, &model.background_component],
            vec!["N".into(), "t0".into(), "background".into()],
        ));

        text.push('\n');
        text.push_str(&format!("{:=^100}\n", " Lifetime Parameters "));
        text.push('\n');

        let lt_norm = model.lifetime_components.iter().map(|comp| comp.intensity.val).sum::<f64>();

        if (lt_norm - 1.).abs() > 1e-3 {
            text.push_str("!!! intensities don't add up to 100% !!!\n\n");
        }

        let mut params = Vec::new();
        let mut names = Vec::new();

        for i in 0..n_l {
            params.push(&model.lifetime_components[i].lifetime);
            names.push(format!("lifetime_{}", i+1));
        }
        for i in 0..n_l {
            params.push(&model.lifetime_components[i].intensity);
            names.push(format!("intensity_{}", i+1));
        }

        text.push_str(&PALSpectrum::print_params(&params, names));

        text.push('\n');
        text.push_str(&format!("{:=^100}\n", " Resolution Parameters "));
        text.push('\n');

        let res_norm = model.resolution_components.iter().map(|comp| comp.intensity.val).sum::<f64>();

        if (res_norm - 1.).abs() > 1e-3 {
            text.push_str("!!! Intensities don't add up to 100% !!!\n");
        }

        let mut params = Vec::new();
        let mut names = Vec::new();

        for i in 0..n_r {
            params.push(&model.resolution_components[i].fwhm);
            names.push(format!("fwhm_{}", i+1));
        }
        for i in 0..n_r {
            params.push(&model.resolution_components[i].intensity);
            names.push(format!("intensity_{}", i+1));
        }
        for i in 0..n_r {
            params.push(&model.resolution_components[i].t0);
            names.push(format!("t0_{}", i+1));
        }

        text.push_str(&PALSpectrum::print_params(&params, names));
        text.push('\n');
        text.push_str(&thin_hline);

        Ok(text)
    }
}
