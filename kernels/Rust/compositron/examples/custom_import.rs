use compositron::core::Measurement;
use compositron::importers::{DataFormat, ImportError};
use thiserror::Error;

#[derive(Debug, Error)]
enum MyCustomImportError {
    #[error("I/O Error")]
    IOError
}

fn import_csv(
    _: &str
) -> Result<Measurement, Box<dyn std::error::Error + Send + Sync>> {
    let m = Measurement::new();

    if true {
        return Err(Box::new(MyCustomImportError::IOError));
    }

    Ok(m)
}

fn main() -> anyhow::Result<(), ImportError> {
    let filename = "";

    let _ = Measurement::from_file(
        filename, DataFormat::Custom { importer: import_csv }
    )?;

    Ok(())
}
