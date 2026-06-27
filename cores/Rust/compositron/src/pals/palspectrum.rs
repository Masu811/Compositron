use std::collections::HashMap;
use crate::core::utils::{Spectrum, spectrum_match};
use crate::pals::model::*;

#[derive(Debug)]
pub enum AnalysisError {

}

pub struct LifetimeParam {
    pub val: f64,
    pub err: f64,
    pub min: f64,
    pub max: f64,
    pub varied: bool,
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
    pub params: HashMap<String, LifetimeParam>,
}

impl PALSpectrum {
    pub fn new(
        spectrum: Spectrum,
        detpair: String,
        tcal: (f64, f64),
        lt_model: Option<LifetimeModel>,
        res_model: Option<LifetimeModel>,
        autocompute: bool,
    ) -> Result<Self, AnalysisError> {
        let counts = spectrum_match!(
            &spectrum, arr => arr.iter().map(|&x| x as u64).sum::<u64>()
        );
        let dcounts = (counts as f64).sqrt();

        let mut p = PALSpectrum {
            spectrum,
            detpair,
            tcal,
            counts,
            dcounts,
            lt_model,
            res_model,
            peak_bnds: None,
            params: HashMap::new(),
        };

        if autocompute {
            PALSpectrum::init(&mut p)?;
        }

        Ok(p)
    }

    fn init(&mut self) -> Result<(), AnalysisError> {
        Ok(())
    }
}
