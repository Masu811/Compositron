use std::io::{BufRead, Write};

use compositron::cdbs::{CDBSpectrum, CDBSpectrumR};
use compositron::core::Measurement;
use compositron::core::utils::{EnergyDetector, EnergyDetectorPair, FromRowMajor, LinearCalibration, Spectrum2D, Unit};
// use compositron::cdbs::cdbspectrum::{
//     Area, LineshapeParamDefinition, ProjectionBins
// };
use compositron::cdbs::cdbspectrumr::{
    Area, LineshapeParamDefinition, ProjectionBins
};
use compositron::importers::{DataFormat, ImportError};
use thiserror::Error;

#[derive(Debug, Error)]
enum MyCustomImportError {
    #[error("I/O Error")]
    IOError {
        #[from]
        inner: std::io::Error
    }
}

fn import_csv(
    filepath: &str
) -> Result<Measurement, Box<dyn std::error::Error + Send + Sync>> {
    let mut m = Measurement::new();

    let path = std::path::Path::new(filepath);

    let filename_components = path.file_stem().unwrap().to_str().unwrap().split("_").collect::<Vec<_>>();

    let l = filename_components.len();

    let scale: f64 = filename_components[l - 1].parse().unwrap();
    let offset: f64 = filename_components[l - 2].parse().unwrap();

    let spectrum_file = std::fs::File::open(filepath)?;
    let reader = std::io::BufReader::new(spectrum_file);

    let mut spectrum = Vec::new();

    for line in reader.lines() {
        let Ok(line) = line else { break; };

        for elem in line.split_whitespace() {
            let x: u64 = elem.parse().unwrap();
            spectrum.push(x);
        }
    }

    let nrows = spectrum.len().isqrt();

    let c = CDBSpectrumR::new(
        Spectrum2D::from_row_major(&spectrum, nrows, nrows),
        EnergyDetectorPair {
            name: "A x B".into(),
            first_det: EnergyDetector {
                name: "A".into(),
                ecal: LinearCalibration { offset, scale },
                corrected_ecal: None,
                eres: None
            },
            second_det: EnergyDetector {
                name: "B".into(),
                ecal: LinearCalibration { offset, scale },
                corrected_ecal: None,
                eres: None
            },
            eres: None
        }
    );

    m.cdbr.insert("A x B".into(), c);

    Ok(m)
}

fn main() -> anyhow::Result<()> {
    // Import

    let mut m = Measurement::from_file(
        // "../../../../testdata/coinc_spectrum_100.000000_0.100000.txt",
        "../../../../testdata/coinc_spectrum_rot_-159.111803_0.070716.txt",
        DataFormat::Custom { importer: import_csv }
    )?;

    let c = m.cdbr.get_mut("A x B").unwrap();

    c.correct_ecal(compositron::core::utils::EcalCorrectionOrder::First)?;

    // let ls_param = LineshapeParamDefinition {
    //     name: "S",
    //     num: &[
    //         Area::Diagonal {
    //             width_cel: Unit::KeV(2.),
    //             width_cml: Unit::KeV(1.),
    //         },
    //     ],
    //     denom: &[
    //         Area::Diagonal {
    //             width_cel: Unit::KeV(10.),
    //             width_cml: Unit::KeV(1.),
    //         },
    //     ],
    //     in_corrected_peak: false,
    // };
    let ls_param = LineshapeParamDefinition {
        name: "S",
        num: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(2.),
                width_cml: Unit::KeV(1.),
            },
        ],
        denom: &[
            Area::AxisAligned {
                width_cel: Unit::KeV(10.),
                width_cml: Unit::KeV(1.),
            },
        ],
        in_corrected_peak: false,
    };

    let s = c.calc_lineshape_param(&ls_param)?;

    println!("Calculated S parameter: {} +/- {}", s.val, s.err);

    // Projection onto diagonal

    // let p = c.project_digonal(
    //     ProjectionBins::Linear(Unit::KeV(0.1)),
    //     Unit::KeV(2.),
    //     Unit::KeV(10.),
    //     false
    // )?;

    // // Write projection to file for inspection

    // let mut f = std::fs::File::create("projection.csv").unwrap();

    // for elem in p.spectrum.iter() {
    //     write!(f, "{elem}\n").unwrap();
    // }

    Ok(())
}
