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
