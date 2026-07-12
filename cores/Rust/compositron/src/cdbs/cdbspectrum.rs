use std::collections::HashMap;

use crate::core::fitting;
use crate::core::utils::Spectrum2D;
use crate::dbs::dbspectrum::LineshapeParam;

pub struct CDBSpectrum {
    pub spectrum: Spectrum2D,
    pub detpair: String,
    pub ecal: ((f64, f64), (f64, f64)),
    pub counts: u64,
    pub dcounts: f64,
    pub eres: Option<f64>,
    pub roi_counts: Option<f64>,
    pub droi_counts: Option<f64>,
    pub roi_params: HashMap<&'static str, fitting::SimpleFitParam>,
    pub lineshape_params: HashMap<String, LineshapeParam>,
}

impl CDBSpectrum {
    pub fn new(
        spectrum: Spectrum2D, detpair: String, ecal: ((f64, f64), (f64, f64))
    ) -> Self {
        CDBSpectrum {
            spectrum: spectrum,
            detpair: detpair,
            ecal: ecal,
            counts: 0,
            dcounts: f64::NAN,
            eres: None,
            roi_counts: None,
            droi_counts: None,
            roi_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }

    pub fn project_axes(&self) {}

    pub fn rebin(&mut self) {}

    pub fn fit_2d_peak(&mut self) {}

    pub fn cut_digonal_slice(&self) {}

    pub fn project_digonal(&self) {}

    pub fn calc_lineshape_param(&mut self) {}
}
