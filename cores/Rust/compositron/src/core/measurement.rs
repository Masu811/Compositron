use std::collections::{BTreeMap, HashMap};
use std::fmt::{self, Display};

use crate::importers::slope_n42_importer::{import_n42, ImportError};
use crate::{cdbs::CDBSpectrum, dbs::DBSpectrum};

#[derive(Debug, Clone, Copy)]
pub struct Shape {
    pub dbs: usize,
    pub cdbs: usize,
}

impl Display for Shape {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Shape(dbs={}, cdbs={})", self.dbs, self.cdbs)
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

    pub fn import(&mut self, filepath: &str) -> Result<(), ImportError> {
        let m = import_n42(filepath)?;

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
