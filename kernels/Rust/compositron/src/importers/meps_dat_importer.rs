use std::fs::File;
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::core::Measurement;
use crate::pals::PALSpectrum;
use crate::dbs::DBSpectrum;
use crate::core::utils::Spectrum;

#[derive(Debug, Error)]
pub enum DATFormatError {
    #[error("Invalid spectrum format")]
    SpectrumParseError {
        #[from]
        inner: std::num::ParseIntError,
    },

    #[error("Invalid header field: {field}")]
    HeaderParseError {
        #[source]
        inner: std::num::ParseFloatError,
        field: String,
    },

    #[error("{info}")]
    HeaderFormatError {
        info: String,
    }
}

#[derive(Debug, Error)]
pub enum ImportError {
    #[error("IO Error opening {path}")]
    FileSystemError {
        #[source]
        inner: std::io::Error,
        path: PathBuf,
    },

    #[error("Invalid MePS file format")]
    InvalidMePSFormatError {
        #[from]
        inner: DATFormatError,
    }
}

struct MePSFileBundle {
    spectrum_file: String,
    param_file: String,
    start_phs_file: String,
    stop_phs_file: String,
}

struct HeaderInfo {
    index: u32,
    sample: String,
    energy: f64,
    temperature: f64,
    bin_width: f64,
}

fn split_filepath_componentes<'a>(filepath: &'a str) -> Option<Vec<&'a str>> {
    let provided_path = Path::new(filepath);
    let provided_filename = provided_path.file_name()?.to_str()?;
    let filename_components: Vec<&str> = provided_filename.split("_").collect();
    Some(filename_components)
}

fn get_all_filenames(filepath: &str) -> Result<MePSFileBundle, ImportError> {
    let filename_components = split_filepath_componentes(filepath).ok_or(
        ImportError::FileSystemError {
            inner: std::io::Error::from(std::io::ErrorKind::InvalidFilename),
            path: filepath.into(),
        }
    )?;

    let critical_idx = filename_components.len() - 3;

    match filename_components[critical_idx] {
        "PALS" => {
            let begin = filename_components[..critical_idx].join("_");
            let end = filename_components[critical_idx+1..].join("_");
            Ok(MePSFileBundle {
                spectrum_file: filename_components.join("_"),
                param_file: [&begin, "Measurement_Parameter", &end].join("_"),
                start_phs_file: [&begin, "PHS_Ch0", &end].join("_"),
                stop_phs_file: [&begin, "PHS_Ch1", &end].join("_")
            })
        },
        "Parameter" => {
            let begin = filename_components[..critical_idx-1].join("_");
            let end = filename_components[critical_idx+1..].join("_");
            Ok(MePSFileBundle {
                spectrum_file: [&begin, "PALS", &end].join("_"),
                param_file: filename_components.join("_"),
                start_phs_file: [&begin, "PHS_Ch0", &end].join("_"),
                stop_phs_file: [&begin, "PHS_Ch1", &end].join("_")
            })
        },
        "Ch0" => {
            let begin = filename_components[..critical_idx-1].join("_");
            let end = filename_components[critical_idx+1..].join("_");
            Ok(MePSFileBundle {
                spectrum_file: [&begin, "PALS", &end].join("_"),
                param_file: [&begin, "Measurement_Parameter", &end].join("_"),
                start_phs_file: filename_components.join("_"),
                stop_phs_file: [&begin, "PHS_Ch1", &end].join("_")
            })
        },
        "Ch1" => {
            let begin = filename_components[..critical_idx-1].join("_");
            let end = filename_components[critical_idx+1..].join("_");
            Ok(MePSFileBundle {
                spectrum_file: [&begin, "PALS", &end].join("_"),
                param_file: [&begin, "Measurement_Parameter", &end].join("_"),
                start_phs_file: [&begin, "PHS_Ch0", &end].join("_"),
                stop_phs_file: filename_components.join("_")
            })
        },
        _ => {
            Err(ImportError::FileSystemError {
                inner: std::io::Error::from(std::io::ErrorKind::InvalidFilename),
                path: filepath.into(),
            })
        },
    }
}

fn parse_header_field(
    field: &str, expected_unit: &str
) -> Result<f64, DATFormatError> {
    let l = field.len();
    let u = expected_unit.len();

    let value = field[..l-u].parse::<f64>().map_err(|err|
        DATFormatError::HeaderParseError {
            inner: err, field: field[..l-u].into()
        },
    )?;
    let unit = &field[l-u..];

    if unit == expected_unit {
        Ok(value)
    } else {
        Err(DATFormatError::HeaderFormatError {
            info: format!(
                "Unexpected unit {unit} in spectrum file header \
                (expected {expected_unit})"
            )
        })
    }
}

fn get_header_info(
    header: &str, m: &mut Measurement
) -> Result<HeaderInfo, DATFormatError> {
    let header_fields: Vec<&str> = header
        .split("_")
        .map(|field| field.trim())
        .collect();

    let header_size = header_fields.len();

    if header_size < 6 {
        return Err(DATFormatError::HeaderFormatError {
            info: "Spectrum header has unexpected format".into()
        });
    }

    let index = header_fields[0].parse::<u32>()?;
    let bin_width = parse_header_field(header_fields[header_size-1], "ns")?;
    let temperature = parse_header_field(header_fields[header_size-3], "K")?;
    let energy = parse_header_field(header_fields[header_size-4], "keV")?;
    let sample = header_fields[1..header_size-5].join("_");

    m.metadata.insert("Sample".into(), sample.clone());
    m.metadata.insert("PositronImplantationEnergy".into(), format!("{energy}"));
    m.metadata.insert("SampleTemperature".into(), format!("{temperature}"));

    Ok(HeaderInfo {
        index,
        sample,
        energy,
        temperature,
        bin_width,
    })
}

fn compare_filename_and_header(filepath: &PathBuf, header: &HeaderInfo) {
    let Some(filename) = filepath.file_name() else { return; };
    let Some(filename) = filename.to_str() else { return; };

    let Some(filename_components) = split_filepath_componentes(filename) else {
        return;
    };

    let l = filename_components.len();

    if let Ok(index) = filename_components[0].parse::<u32>()
        && index != header.index
    {
        // log::warn
    }

    if let Ok(energy) = parse_header_field(filename_components[l-6], "keV")
        && energy != header.energy
    {

    }

    if let Ok(temperature) = parse_header_field(filename_components[l-5], "K")
        && temperature != header.temperature
    {
        // log::warn
    }

    let sample = filename_components[1..l-7].join("_");

    if sample != header.sample {
        // log::warn
    }
}

fn import_spectrum_data(
    channel_data: &str
) -> Result<Spectrum, DATFormatError> {
    let channel_data = channel_data.split_whitespace();

    let mut spectrum: Vec<u32> = Vec::new();

    for x in channel_data {
        let x = x.parse::<u32>()?;
        spectrum.push(x);
    }

    Ok(Spectrum::from(spectrum))
}

fn import_palspectrum(filepath: &PathBuf, m: &mut Measurement) -> Result<PALSpectrum, ImportError> {
    let spectrum_file = File::open(filepath).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: filepath.into(),
        }
    )?;

    let mut reader = std::io::BufReader::new(spectrum_file);

    let mut header = String::new();
    reader.read_line(&mut header).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: filepath.into(),
        }
    )?;

    let mut channel_data = String::new();
    reader.read_to_string(&mut channel_data).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: filepath.into(),
        }
    )?;

    let header_info = get_header_info(&header, m)?;

    compare_filename_and_header(filepath, &header_info);

    let spectrum = import_spectrum_data(&channel_data)?;

    Ok(PALSpectrum::new(
        spectrum, "A".into(), (0., 1000. * header_info.bin_width)
    ))
}

fn import_params(filepath: &PathBuf, m: &mut Measurement) {
    let Ok(param_file) = File::open(filepath) else {
        // log::warn
        return;
    };

    let mut reader = std::io::BufReader::new(param_file);

    let mut content = String::new();
    let Ok(_) = reader.read_to_string(&mut content) else {
        return;
    };

    let lines = content.split("\n").map(|line| line.trim());

    for line in lines {
        if line.len() < 30 {
            continue;
        }

        let key = line[..29].trim();
        let value = line[29..].trim();

        m.metadata.insert(key.into(), value.into());
    }
}

fn import_ph_spectrum(filepath: &PathBuf) -> Option<DBSpectrum> {
    let Ok(param_file) = File::open(filepath) else {
        // log::warn
        return None;
    };

    let mut reader = std::io::BufReader::new(param_file);

    let mut spectrum_data = String::new();
    let Ok(_) = reader.read_to_string(&mut spectrum_data) else {
        return None;
    };

    let Ok(spectrum) = import_spectrum_data(&spectrum_data) else {
        return None;
    };

    Some(DBSpectrum::new(spectrum, "A".into(), (f64::NAN, f64::NAN), None))
}

pub fn import_meps_dat(filepath: &str) -> Result<Measurement, ImportError> {
    let datafiles = get_all_filenames(filepath)?;

    let dir = Path::new(filepath).parent().unwrap_or(Path::new(""));

    let mut m = Measurement::new();

    let p = import_palspectrum(&dir.join(datafiles.spectrum_file), &mut m)?;

    m.pals.insert("A".into(), p);

    import_params(&dir.join(datafiles.param_file), &mut m);

    if let Some(mut d1) = import_ph_spectrum(&dir.join(datafiles.start_phs_file)) {
        d1.detname = "Start".into();
        m.dbs.insert("Start".into(), d1);
    }

    if let Some(mut d2) = import_ph_spectrum(&dir.join(datafiles.stop_phs_file)) {
        d2.detname = "Stop".into();
        m.dbs.insert("Stop".into(), d2);
    }

    Ok(m)
}
