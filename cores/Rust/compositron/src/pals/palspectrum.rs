use std::collections::HashMap;

use crate::core::utils::{Spectrum, spectrum_match};
use crate::core::fitting::FitParam;
use crate::pals::model::*;

#[derive(Debug)]
pub enum AnalysisError {

}

pub struct PALSpectrum {
    pub spectrum: Spectrum,
    pub detpair: String,
    pub tcal: (f64, f64),
    pub counts: u64,
    pub dcounts: f64,
    pub lt_model: Option<LifetimeModel>,
    pub res_model: Option<LifetimeModel>,
    pub peak_bnds: Option<(usize, usize)>,
    pub params: HashMap<String, FitParam>,
}

impl PALSpectrum {
    pub fn new(
        spectrum: Spectrum,
        detpair: String,
        tcal: (f64, f64),
    ) -> Self {
        let counts = spectrum_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );
        let dcounts = (counts as f64).sqrt();

        PALSpectrum {
            spectrum,
            detpair,
            tcal,
            counts,
            dcounts,
            lt_model: None,
            res_model: None,
            peak_bnds: None,
            params: HashMap::new(),
        }
    }

    fn fit(
        &mut self,
        lt_model: LifetimeModel,
        res_model: LifetimeModel,
    ) -> Result<(), AnalysisError> {
        Ok(())
    }

    fn fit_report(
        &self,
        file: &mut impl std::fmt::Write,
    ) {

    }
}
