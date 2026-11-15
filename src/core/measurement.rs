use std::collections::HashMap;
use std::fmt::{self, Display};

use crate::cdbs::cdbspectrum::CDBSpectrumArray;
use crate::dbs::dbspectrum::DBSpectrumArray;
use crate::importers::slope_n42_importer::import_n42;
use crate::importers::meps_dat_importer::import_meps_dat;
use crate::importers::{DataFormat, ImportError};
use crate::pals::palspectrum::PALSpectrumArray;
use crate::{cdbs::CDBSpectrum, dbs::DBSpectrum, pals::PALSpectrum};


#[derive(Debug, Clone, Copy)]
pub struct Shape {
    pub dbs: usize,
    pub cdbs: usize,
    pub pals: usize,
}


impl Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "Shape(dbs={}, cdbs={}, pals={})", self.dbs, self.cdbs, self.pals
        )
    }
}


pub struct Measurement {
    pub filename: Option<String>,
    pub name: Option<String>,
    pub dbs: HashMap<String, DBSpectrum>,
    pub cdbs: HashMap<String, CDBSpectrum>,
    pub pals: HashMap<String, PALSpectrum>,
    pub metadata: HashMap<String, String>,
}


impl Measurement {
    pub fn new() -> Self {
        Measurement {
            filename: None,
            name: None,
            dbs: HashMap::new(),
            cdbs: HashMap::new(),
            pals: HashMap::new(),
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

    pub fn shape(&self) -> Shape {
        Shape {
            dbs: self.dbs.len(),
            cdbs: self.cdbs.len(),
            pals: self.pals.len(),
        }
    }
}


pub struct MeasurementList {
    pub measurements: Vec<Measurement>,
    pub name: Option<String>,
}


impl MeasurementList {
    pub fn new() -> Self {
        MeasurementList {
            name: None,
            measurements: Vec::new(),
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

        Ok(mlist)
    }
}


pub struct MeasurementArray {
    pub name: Option<String>,
    pub dbs: HashMap<String, DBSpectrumArray>,
    pub cdbs: HashMap<String, CDBSpectrumArray>,
    pub pals: HashMap<String, PALSpectrumArray>,
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
        let ext = match format {
            DataFormat::SlopeN42 => "n42",
            DataFormat::MePSDat => "dat",
            DataFormat::Custom { importer: _, extension } => extension,
        };

        let files = std::fs::read_dir(path)?
            .filter_map(Result::ok)
            .filter(|x| x.path().extension().and_then(|ext| ext.to_str()) == Some(ext));

        let n_files = files.count();

        let mut marr = MeasurementArray::new();

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

            for (det_name, d) in m.dbs.into_iter() {
                if !marr.dbs.contains_key(&det_name) {
                    let det = d.detector.clone();
                    marr.dbs.insert(det_name.clone(), DBSpectrumArray::new(det));
                }
                marr.dbs.get_mut(&det_name).unwrap().push(d);
            }

            for (det_name, c) in m.cdbs.into_iter() {
                if !marr.cdbs.contains_key(&det_name) {
                    let det = c.detpair.clone();
                    marr.cdbs.insert(
                        det_name.clone(), CDBSpectrumArray::new(det, c.orientation)
                    );
                }
                marr.cdbs.get_mut(&det_name).unwrap().push(c);
            }

            for (det_name, p) in m.pals.into_iter() {
                if !marr.pals.contains_key(&det_name) {
                    let det = p.detpair.clone();
                    marr.pals.insert(det_name.clone(), PALSpectrumArray::new(det));
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
