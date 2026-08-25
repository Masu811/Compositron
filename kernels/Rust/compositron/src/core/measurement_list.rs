use crate::core::measurement::Measurement;
use crate::importers::{DataFormat, ImportError};

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

    pub fn import(
        &mut self, path: &str, format: DataFormat
    ) -> Result<(), ImportError> {
        Ok(())
    }
}
