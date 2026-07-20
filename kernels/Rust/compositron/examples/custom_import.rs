use compositron::core::{Measurement, measurement::{DataFormat, ImportError}};
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

    let mut m = Measurement::new();

    m.import(
        filename, DataFormat::Custom { importer: import_csv }
    )?;

    Ok(())
}
