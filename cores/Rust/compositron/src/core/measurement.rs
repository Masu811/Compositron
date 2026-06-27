use std::collections::{BTreeMap, HashMap};
use std::fmt::{self, Display};

use thiserror::Error;

use crate::importers::slope_n42_importer::{self, import_n42};
use crate::{cdbs::CDBSpectrum, dbs::DBSpectrum};

#[derive(Debug, Error)]
pub enum ImportError {
    #[error("Error during data import")]
    SlopeN42ImportError {
        #[from]
        source: slope_n42_importer::ImportError,
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
    Custom {
        importer: fn (&str) -> Result<Measurement, Box<dyn std::error::Error + Send + Sync>>,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Shape {
    pub dbs: usize,
    pub cdbs: usize,
}

impl Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f, "Shape(dbs={}, cdbs={})", self.dbs, self.cdbs,
        )
    }
}

pub struct Measurement {
    pub filename: Option<String>,
    pub name: Option<String>,
    pub dbs: BTreeMap<String, DBSpectrum>,
    pub cdbs: BTreeMap<String, CDBSpectrum>,
    pub metadata: HashMap<String, String>,
}

impl Measurement {
    pub fn new() -> Self {
        Measurement {
            filename: None,
            name: None,
            dbs: BTreeMap::new(),
            cdbs: BTreeMap::new(),
            metadata: HashMap::new(),
        }
    }

    pub fn import(
        &mut self, filepath: &str, format: DataFormat
    ) -> Result<(), ImportError> {
        let m = match format {
            DataFormat::SlopeN42 => import_n42(filepath)?,
            DataFormat::Custom { importer } => importer(filepath).map_err(
                |err| ImportError::CustomImporterError { source: err }
            )?
        };

        self.filename = Some(filepath.into());
        self.name = self.filename.clone();
        self.dbs = m.dbs;
        self.cdbs = m.cdbs;
        self.metadata = m.metadata;

        Ok(())
    }

    pub fn shape(&self) -> Shape {
        Shape {
            dbs: self.dbs.len(),
            cdbs: self.cdbs.len(),
        }
    }
}
