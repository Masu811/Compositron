pub mod meps_dat_importer;
pub mod png_importer;
pub mod slope_n42_importer;

use thiserror::Error;
use crate::core::Measurement;

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
