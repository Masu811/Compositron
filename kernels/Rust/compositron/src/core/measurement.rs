use std::collections::{BTreeMap, HashMap};
use std::fmt::{self, Display};

use thiserror::Error;

use crate::importers::slope_n42_importer::{self, import_n42};
use crate::importers::meps_dat_importer::{self, import_meps_dat};
use crate::{cdbs::CDBSpectrum, dbs::DBSpectrum, pals::PALSpectrum};

#[derive(Debug, Error)]
pub enum ImportError {
    #[error("Error during data import")]
    SlopeN42ImportError {
        #[from]
        source: slope_n42_importer::ImportError,
    },

    #[error("Error during data import")]
    MePSDatImportError {
        #[from]
        source: meps_dat_importer::ImportError,
    },

    #[error("Error during data import")]
    CustomImporterError {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

#[derive(Debug, Clone, Copy)]
pub enum DataFormat {
    SlopeN42,
    MePSDat,
    Custom {
        importer: fn (&str) -> Result<Measurement, Box<dyn std::error::Error + Send + Sync>>,
    }
}

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
    pub dbs: BTreeMap<String, DBSpectrum>,
    pub cdbs: BTreeMap<String, CDBSpectrum>,
    pub pals: BTreeMap<String, PALSpectrum>,
    pub metadata: HashMap<String, String>,
}

impl Measurement {
    pub fn new() -> Self {
        Measurement {
            filename: None,
            name: None,
            dbs: BTreeMap::new(),
            cdbs: BTreeMap::new(),
            pals: BTreeMap::new(),
            metadata: HashMap::new(),
        }
    }

    pub fn import(
        &mut self, filepath: &str, format: DataFormat
    ) -> Result<(), ImportError> {
        let m = match format {
            DataFormat::SlopeN42 => import_n42(filepath)?,
            DataFormat::MePSDat => import_meps_dat(filepath)?,
            DataFormat::Custom { importer } => importer(filepath).map_err(
                |err| ImportError::CustomImporterError { source: err }
            )?
        };

        self.filename = Some(filepath.into());
        self.name = self.filename.clone();
        self.dbs = m.dbs;
        self.cdbs = m.cdbs;
        self.pals = m.pals;
        self.metadata = m.metadata;

        Ok(())
    }

    pub fn shape(&self) -> Shape {
        Shape {
            dbs: self.dbs.len(),
            cdbs: self.cdbs.len(),
            pals: self.pals.len(),
        }
    }
}
