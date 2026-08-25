use std::collections::{HashMap, BTreeMap};

use nalgebra::{DVector, DMatrix};

use crate::importers::{DataFormat, ImportError};
use crate::core::utils::{
    Spectrum, Spectrum2D, EnergyDetector, EnergyDetectorPair, TimingDetectorPair
};
use crate::core::fitting;
use crate::dbs::dbspectrum;
use crate::cdbs::cdbspectrum;
use crate::pals::fitting::FitResult;

pub struct DBSpectra {
    pub detector: EnergyDetector,
    pub spectra: Vec<Spectrum>,
    pub counts: Vec<u64>,
    pub dcounts: Vec<f64>,
    pub peaks: Option<DVector<f64>>,
    pub peak_bnds: Option<Vec<(usize, usize)>>,
    pub peak_counts: Vec<f64>,
    pub dpeak_counts: Vec<f64>,
    pub peak_params: HashMap<&'static str, Vec<fitting::SimpleFitParam>>,
    pub lineshape_params: HashMap<String, Vec<dbspectrum::LineshapeParam>>,
}

pub struct CDBSpectra {
    pub detpair: EnergyDetectorPair,
    pub spectra: Vec<Spectrum2D>,
    pub counts: Vec<u64>,
    pub dcounts: Vec<f64>,
    pub peaks: Option<Vec<DMatrix<f64>>>,
    pub peak_bnds: Option<Vec<((usize, usize), (usize, usize))>>,
    pub peak_counts: Option<Vec<f64>>,
    pub dpeak_counts: Option<Vec<f64>>,
    pub peak_params: HashMap<&'static str, Vec<fitting::SimpleFitParam>>,
    pub lineshape_params: HashMap<String, Vec<cdbspectrum::LineshapeParam>>,
}

pub struct PALSpectra {
    pub detpair: TimingDetectorPair,
    pub spectra: Vec<Spectrum>,
    pub counts: Vec<u64>,
    pub dcounts: Vec<f64>,
    pub peak_bnds: Option<Vec<(usize, usize)>>,
    pub peak_centers: Option<Vec<f64>>,
    pub dpeak_centers: Option<Vec<f64>>,
    pub peak_fwhms: Option<Vec<f64>>,
    pub dpeak_fwhms: Option<Vec<f64>>,
    pub fit_results: Option<Vec<FitResult>>,
}

pub struct MeasurementArray {
    pub name: Option<String>,
    pub dbs: BTreeMap<String, DBSpectra>,
    pub cdbs: BTreeMap<String, CDBSpectra>,
    pub pals: BTreeMap<String, PALSpectra>,
    pub metadata: HashMap<String, Vec<String>>,
}

impl MeasurementArray {
    pub fn new() -> Self {
        MeasurementArray {
            name: None,
            dbs: BTreeMap::new(),
            cdbs: BTreeMap::new(),
            pals: BTreeMap::new(),
            metadata: HashMap::new(),
        }
    }
}
