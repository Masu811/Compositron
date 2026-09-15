use std::collections::HashMap;

use nalgebra::{DVector, DMatrix};

use crate::cdbs::cdbspectrum::Orientation;
use crate::core::fitting::SimpleFitParam;
use crate::core::utils::{
    Spectrum, Spectrum2D, EnergyDetector, EnergyDetectorPair, TimingDetectorPair
};
use crate::core::{Measurement, fitting};
use crate::dbs::{DBSpectrum, dbspectrum};
use crate::cdbs::{CDBSpectrum, cdbspectrum};
use crate::importers::{DataFormat, ImportError};
use crate::pals::PALSpectrum;
use crate::pals::fitting::FitResult;


pub struct DBSpectra {
    pub detector: EnergyDetector,
    pub spectra: Vec<Spectrum>,
    pub counts: Vec<u64>,
    pub peaks: Option<DVector<f64>>,
    pub peak_bnds: Option<Vec<(usize, usize)>>,
    pub peak_counts: Option<Vec<f64>>,
    pub peak_params: HashMap<&'static str, Vec<fitting::SimpleFitParam>>,
    pub lineshape_params: HashMap<String, Vec<dbspectrum::LineshapeParam>>,
}


impl DBSpectra {
    pub fn new(detector: EnergyDetector) -> Self {
        DBSpectra {
            detector,
            spectra: Vec::new(),
            counts: Vec::new(),
            peaks: None,
            peak_bnds: None,
            peak_counts: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }


    pub fn push(&mut self, d: DBSpectrum) {
        self.spectra.push(d.spectrum);
        self.counts.push(d.counts);
    }


    pub fn len(&self) -> usize {
        self.spectra.len()
    }
}


pub struct CDBSpectra {
    pub detpair: EnergyDetectorPair,
    pub spectra: Vec<Spectrum2D>,
    pub orientation: Orientation,
    pub counts: Vec<u64>,
    pub peaks: Option<Vec<DMatrix<f64>>>,
    pub peak_bnds: Option<Vec<((usize, usize), (usize, usize))>>,
    pub peak_counts: Option<Vec<f64>>,
    pub peak_params: HashMap<&'static str, Vec<fitting::SimpleFitParam>>,
    pub lineshape_params: HashMap<String, Vec<cdbspectrum::LineshapeParam>>,
}


impl CDBSpectra {
    pub fn new(detpair: EnergyDetectorPair, orientation: Orientation) -> Self {
        CDBSpectra {
            detpair,
            spectra: Vec::new(),
            orientation,
            counts: Vec::new(),
            peaks: None,
            peak_bnds: None,
            peak_counts: None,
            peak_params: HashMap::new(),
            lineshape_params: HashMap::new(),
        }
    }


    pub fn push(&mut self, c: CDBSpectrum) {
        self.spectra.push(c.spectrum);
        self.counts.push(c.counts);
    }


    pub fn len(&self) -> usize {
        self.spectra.len()
    }
}


pub struct PALSpectra {
    pub detpair: TimingDetectorPair,
    pub spectra: Vec<Spectrum>,
    pub counts: Vec<u64>,
    pub peak_params: HashMap<&'static str, Vec<SimpleFitParam>>,
    pub fit_results: Option<Vec<FitResult>>,
}


impl PALSpectra {
    pub fn new(detpair: TimingDetectorPair) -> Self {
        PALSpectra {
            detpair,
            spectra: Vec::new(),
            counts: Vec::new(),
            peak_params: HashMap::new(),
            fit_results: None,
        }
    }


    pub fn push(&mut self, p: PALSpectrum) {
        self.spectra.push(p.spectrum);
        self.counts.push(p.counts);
    }


    pub fn len(&self) -> usize {
        self.spectra.len()
    }
}


pub struct MeasurementArray {
    pub name: Option<String>,
    pub dbs: HashMap<String, DBSpectra>,
    pub cdbs: HashMap<String, CDBSpectra>,
    pub pals: HashMap<String, PALSpectra>,
    pub metadata: HashMap<String, Vec<String>>,
}


impl MeasurementArray {
    pub fn new() -> Self {
        MeasurementArray {
            name: None,
            dbs: HashMap::new(),
            cdbs: HashMap::new(),
            pals: HashMap::new(),
            metadata: HashMap::new(),
        }
    }


    pub fn from_dir(
        path: &str, format: &DataFormat
    ) -> Result<Self, ImportError> {
        let files = std::fs::read_dir(path)?;

        let n_files = files.count();

        let mut marr = MeasurementArray::new();

        let files = std::fs::read_dir(path)?;

        for (i, file) in files.enumerate() {
            let m = Measurement::from_file(
                file?
                    .path()
                    .to_str()
                    .ok_or(std::io::Error::from(std::io::ErrorKind::InvalidFilename))?,
                format
            )?;

            for (det_name, d) in m.dbs.into_iter() {
                if !marr.dbs.contains_key(&det_name) {
                    let det = d.detector.clone();
                    marr.dbs.insert(det_name.clone(), DBSpectra::new(det));
                }
                marr.dbs.get_mut(&det_name).unwrap().push(d);
            }

            for (det_name, c) in m.cdbs.into_iter() {
                if !marr.cdbs.contains_key(&det_name) {
                    let det = c.detpair.clone();
                    marr.cdbs.insert(
                        det_name.clone(), CDBSpectra::new(det, c.orientation)
                    );
                }
                marr.cdbs.get_mut(&det_name).unwrap().push(c);
            }

            for (det_name, p) in m.pals.into_iter() {
                if !marr.pals.contains_key(&det_name) {
                    let det = p.detpair.clone();
                    marr.pals.insert(det_name.clone(), PALSpectra::new(det));
                }
                marr.pals.get_mut(&det_name).unwrap().push(p);
            }

            let bar_size = 20;
            let frac_complete = (i as f64 + 1.) / n_files as f64;
            let progress = "=".repeat((bar_size as f64 * frac_complete) as usize);
            let empty = " ".repeat(bar_size - progress.len());

            print!(
                "\rImporting: [{}{}] {}/{} ({}%)",
                progress,
                empty,
                i + 1,
                n_files,
                (100. * frac_complete) as usize
            );
        }

        Ok(marr)
    }
}
