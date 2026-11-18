use std::collections::HashMap;

use crate::cdbs::cdbspectrum::{CDBSpectrumArray, Projection, ProjectionRatio};
use crate::dbs::dbspectrum::DBSpectrumArray;
use crate::eventlist::eventlist::{CoincidenceEventList, CorrelatedEventLists, EventList};
use crate::importers::slope_n42_importer::import_n42;
use crate::importers::meps_dat_importer::import_meps_dat;
use crate::importers::{DataFormat, ImportError};
use crate::pals::palspectrum::PALSpectrumArray;
use crate::waveforms::waveform::WaveformArray;
use crate::{cdbs::CDBSpectrum, dbs::DBSpectrum, pals::PALSpectrum};


pub struct DBSMeasurement {
    pub spectra: HashMap<String, DBSpectrum>,
}


pub struct CDBSMeasurement {
    pub spectra: HashMap<String, CDBSpectrum>,
    pub projections: HashMap<String, Projection>,
    pub ratios: HashMap<String, ProjectionRatio>,
}


pub struct PALSMeasurement {
    pub spectra: HashMap<String, PALSpectrum>,
}


pub struct EventListMeasurement {
    pub single_eventlists: HashMap<String, EventList>,
    pub corr_eventlists: HashMap<String, CorrelatedEventLists>,
    pub coinc_eventlists: HashMap<String, CoincidenceEventList>,
}


pub struct WaveformMeasurement {
    pub arrays: HashMap<String, WaveformArray>,
}


pub struct Measurement {
    pub name: Option<String>,
    pub dbs: DBSMeasurement,
    pub cdbs: CDBSMeasurement,
    pub pals: PALSMeasurement,
    pub eventlists: EventListMeasurement,
    pub waveforms: WaveformMeasurement,
    pub metadata: HashMap<String, String>,
}


pub struct MeasurementList {
    pub name: Option<String>,
    pub measurements: Vec<Measurement>,
    pub metadata: HashMap<String, Vec<String>>,
}


pub struct DBSMeasurementArray {
    pub spectra: HashMap<String, DBSpectrumArray>,
}


pub struct CDBSMeasurementArray {
    pub spectra: HashMap<String, CDBSpectrumArray>,
    pub projections: HashMap<String, Vec<Projection>>,
    pub ratios: HashMap<String, Vec<ProjectionRatio>>,
}


pub struct PALSMeasurementArray {
    pub spectra: HashMap<String, PALSpectrumArray>,
}


pub struct EventListMeasurementArray {
    pub single_eventlists: HashMap<String, Vec<EventList>>,
    pub corr_eventlists: HashMap<String, Vec<CorrelatedEventLists>>,
    pub coinc_eventlists: HashMap<String, Vec<CoincidenceEventList>>,
}


pub struct WaveformMeasurementArray {
    pub arrays: HashMap<String, Vec<WaveformArray>>,
}


pub struct MeasurementArray {
    pub name: Option<String>,
    pub dbs: DBSMeasurementArray,
    pub cdbs: CDBSMeasurementArray,
    pub pals: PALSMeasurementArray,
    pub eventlists: EventListMeasurementArray,
    pub waveforms: WaveformMeasurementArray,
    pub metadata: HashMap<String, String>,
}


impl Measurement {
    pub fn new() -> Self {
        Measurement {
            name: None,
            dbs: DBSMeasurement {
                spectra: HashMap::new(),
            },
            cdbs: CDBSMeasurement {
                spectra: HashMap::new(),
                projections: HashMap::new(),
                ratios: HashMap::new(),
            },
            pals: PALSMeasurement {
                spectra: HashMap::new(),
            },
            eventlists: EventListMeasurement {
                single_eventlists: HashMap::new(),
                corr_eventlists: HashMap::new(),
                coinc_eventlists: HashMap::new(),
            },
            waveforms: WaveformMeasurement {
                arrays: HashMap::new(),
            },
            metadata: HashMap::new(),
        }
    }

    pub fn from_file(
        filepath: &str, format: &DataFormat
    ) -> Result<Self, ImportError> {
        Ok(match format {
            DataFormat::SlopeN42 => import_n42(filepath)?,
            DataFormat::MePSDat => import_meps_dat(filepath)?,
            DataFormat::Custom { importer, extension: _ } => {
                importer(filepath).map_err(
                    |err| ImportError::CustomImporterError { source: err }
                )?
            }
        })
    }
}


impl MeasurementList {
    pub fn new() -> Self {
        MeasurementList {
            name: None,
            measurements: Vec::new(),
            metadata: HashMap::new(),
        }
    }

    pub fn from_dir(
        path: &str, format: &DataFormat
    ) -> Result<Self, ImportError> {
        let ext = match format {
            DataFormat::SlopeN42 => "n42",
            DataFormat::MePSDat => "dat",
            DataFormat::Custom { importer: _, extension } => extension,
        };

        let files = std::fs::read_dir(path)?
            .filter_map(Result::ok)
            .filter(|x| x.path().extension().and_then(|ext| ext.to_str()) == Some(ext));

        let n_files = files.count();

        let mut mlist = MeasurementList::new();
        mlist.measurements.reserve(n_files);

        let mut files = std::fs::read_dir(path)?
            .filter_map(Result::ok)
            .filter(|x| x.path().extension().and_then(|ext| ext.to_str()) == Some(ext))
            .collect::<Vec<_>>();

        files.sort_by(|a, b| natord::compare(
            &a.file_name().to_string_lossy(),
            &b.file_name().to_string_lossy(),
        ));

        for (i, file) in files.iter().enumerate() {
            let m = Measurement::from_file(
                file
                    .path()
                    .to_str()
                    .ok_or(std::io::Error::from(std::io::ErrorKind::InvalidFilename))?,
                format
            )?;

            mlist.measurements.push(m);

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

        println!();

        mlist.name = std::path::Path::new(path)
            .file_name()
            .and_then(|name| Some(name.to_string_lossy().to_string()));

        Ok(mlist)
    }


    pub fn transpose(self) -> MeasurementArray {
        let mut marr = MeasurementArray::new();

        marr.name = self.name;

        for m in self.measurements {
            marr.push(m);
        }

        marr
    }
}


impl MeasurementArray {
    pub fn new() -> Self {
        MeasurementArray {
            name: None,
            dbs: DBSMeasurementArray {
                spectra: HashMap::new(),
            },
            cdbs: CDBSMeasurementArray {
                spectra: HashMap::new(),
                projections: HashMap::new(),
                ratios: HashMap::new(),
            },
            pals: PALSMeasurementArray {
                spectra: HashMap::new(),
            },
            eventlists: EventListMeasurementArray {
                single_eventlists: HashMap::new(),
                corr_eventlists: HashMap::new(),
                coinc_eventlists: HashMap::new(),
            },
            waveforms: WaveformMeasurementArray {
                arrays: HashMap::new(),
            },
            metadata: HashMap::new(),
        }
    }


    pub fn from_dir(
        path: &str, format: &DataFormat
    ) -> Result<Self, ImportError> {
        Ok(MeasurementList::from_dir(path, format)?.transpose())
    }


    pub fn push(&mut self, m: Measurement) {
        for (det, dbspectrum) in m.dbs.spectra.into_iter() {
            if !self.dbs.spectra.contains_key(&det) {
                self.dbs.spectra.insert(
                    det.clone(),
                    DBSpectrumArray::new(dbspectrum.detector.clone())
                );
            }
            self.dbs.spectra.get_mut(&det).unwrap().push(dbspectrum).unwrap();
        }
        for (det, cdbspectrum) in m.cdbs.spectra.into_iter() {
            if !self.cdbs.spectra.contains_key(&det) {
                self.cdbs.spectra.insert(
                    det.clone(),
                    CDBSpectrumArray::new(cdbspectrum.detpair.clone(), cdbspectrum.orientation)
                );
            }
            self.cdbs.spectra.get_mut(&det).unwrap().push(cdbspectrum).unwrap();
        }
        for (det, proj) in m.cdbs.projections.into_iter() {
            if !self.cdbs.projections.contains_key(&det) {
                self.cdbs.projections.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.cdbs.projections.get_mut(&det).unwrap().push(proj);
        }
        for (det, ratio) in m.cdbs.ratios.into_iter() {
            if !self.cdbs.ratios.contains_key(&det) {
                self.cdbs.ratios.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.cdbs.ratios.get_mut(&det).unwrap().push(ratio);
        }
        for (det, palspectrum) in m.pals.spectra.into_iter() {
            if !self.pals.spectra.contains_key(&det) {
                self.pals.spectra.insert(
                    det.clone(),
                    PALSpectrumArray::new(palspectrum.detpair.clone())
                );
            }
            self.pals.spectra.get_mut(&det).unwrap().push(palspectrum).unwrap();
        }
        for (det, elist) in m.eventlists.single_eventlists.into_iter() {
            if !self.eventlists.single_eventlists.contains_key(&det) {
                self.eventlists.single_eventlists.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.eventlists.single_eventlists.get_mut(&det).unwrap().push(elist);
        }
        for (det, elist) in m.eventlists.corr_eventlists.into_iter() {
            if !self.eventlists.corr_eventlists.contains_key(&det) {
                self.eventlists.corr_eventlists.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.eventlists.corr_eventlists.get_mut(&det).unwrap().push(elist);
        }
        for (det, elist) in m.eventlists.coinc_eventlists.into_iter() {
            if !self.eventlists.coinc_eventlists.contains_key(&det) {
                self.eventlists.coinc_eventlists.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.eventlists.coinc_eventlists.get_mut(&det).unwrap().push(elist);
        }
        for (det, warr) in m.waveforms.arrays.into_iter() {
            if !self.waveforms.arrays.contains_key(&det) {
                self.waveforms.arrays.insert(
                    det.clone(),
                    Vec::new(),
                );
            }
            self.waveforms.arrays.get_mut(&det).unwrap().push(warr);
        }
    }
}
