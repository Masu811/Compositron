use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use png;
use thiserror::Error;

use crate::core::utils::{Spectrum2D, FromRowMajor};

#[derive(Debug, Error)]
pub enum ImportError {
    #[error("Error while opening PNG file {path}")]
    FileSystemError {
        #[source]
        inner: std::io::Error,
        path: PathBuf,
    },

    #[error("Could not decode PNG file {path}")]
    DecodingError {
        #[source]
        inner: png::DecodingError,
        path: PathBuf,
    },

    #[error("PNG file {path} has unsupported bit depth {depth}")]
    PNGFormatError {
        depth: u8,
        path: PathBuf,
    },

    #[error("Insufficient memory to open PNG file {path}")]
    MemoryError {
        path: PathBuf,
    }
}

fn bytes_to_spectrum(bytes: &[u8], n_channels: usize) -> Vec<u32> {
    let mut spectrum = Vec::with_capacity(bytes.len() / n_channels);

    match n_channels {
        1 => {
            spectrum.extend(bytes.iter().map(|&b| b as u32));
        }
        2 => {
            let mut i = 0;
            while i + 1 <= bytes.len() {
                let b0 = bytes[i] as u32;
                let b1 = bytes[i + 1] as u32;
                spectrum.push(b0 | (b1 << 8));
                i += 2;
            }
        }
        3 => {
            let mut i = 0;
            while i + 3 <= bytes.len() {
                let b0 = bytes[i] as u32;
                let b1 = bytes[i + 1] as u32;
                let b2 = bytes[i + 2] as u32;
                spectrum.push(b0 | (b1 << 8) | (b2 << 16));
                i += 3;
            }
        }
        4 => {
            let mut i = 0;
            while i + 4 <= bytes.len() {
                let b0 = bytes[i] as u32;
                let b1 = bytes[i + 1] as u32;
                let b2 = bytes[i + 2] as u32;
                let b3 = bytes[i + 3] as u32;
                spectrum.push(b0 | (b1 << 8) | (b2 << 16) | (b3 << 24));
                i += 4;
            }
        }
        _ => unreachable!(),
    }

    spectrum
}

pub fn import_png(path: &Path) -> Result<Spectrum2D, ImportError> {
    let file = File::open(path).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: path.into(),
        }
    )?;

    let decoder = png::Decoder::new(BufReader::new(file));

    let mut reader = decoder.read_info().map_err(
        |err| ImportError::DecodingError {
            inner: err,
            path: path.into(),
        }
    )?;

    let Some(buffer_size) = reader.output_buffer_size() else {
        return Err(ImportError::MemoryError {
            path: path.into(),
        });
    };

    let mut buf = vec![0; buffer_size];

    let info = reader.next_frame(&mut buf).map_err(
        |err| ImportError::DecodingError {
            inner: err,
            path: path.into(),
        }
    )?;

    let n_channels = match info.color_type {
        png::ColorType::Grayscale | png::ColorType::Indexed => 1,
        png::ColorType::GrayscaleAlpha => 2,
        png::ColorType::Rgb => 3,
        png::ColorType::Rgba => 4,
    };

    if info.bit_depth != png::BitDepth::Eight {
        return Err(ImportError::PNGFormatError {
            depth: info.bit_depth as u8,
            path: path.into(),
        });
    }

    let bytes: &[u8] = &buf[..info.buffer_size()];

    let spectrum = bytes_to_spectrum(bytes, n_channels);

    Ok(Spectrum2D::from_row_major(
        &spectrum, reader.info().height as usize, reader.info().width as usize
    ))
}
