use std::io::Write;
use std::time::Instant;

use compositron::core::Measurement;
use compositron::core::utils::{Unit, EcalCorrectionOrder};
use compositron::cdbs::cdbspectrum::{
    Area, LineshapeParamDefinition, ProjectionBins
};
use compositron::importers::DataFormat;

fn main() -> anyhow::Result<()> {
    // Import

    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/depth-profile_Copper_0000.n42",
        DataFormat::SlopeN42
    )?;

    let c = m.cdbs.get_mut("OAA x OAB").unwrap();

    c.correct_ecal(EcalCorrectionOrder::First)?;

    // S parameter calculation

    let ls_param = LineshapeParamDefinition {
        name: "S",
        num: &[
            Area::Diagonal {
                width_cel: Unit::KeV(2.),
                width_cml: Unit::KeV(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        denom: &[
            Area::Diagonal {
                width_cel: Unit::KeV(10.),
                width_cml: Unit::KeV(1.),
                offset_cel: Unit::KeV(0.),
                offset_cml: Unit::KeV(0.),
            },
        ],
        in_corrected_peak: false,
    };

    let start = Instant::now();
    let s = c.calc_lineshape_param(&ls_param)?;
    let duration = start.elapsed();

    println!("Calculated S parameter: {} +/- {}", s.val, s.err);
    println!("Calculation duration: {:?}", duration);

    // Projection onto diagonal

    let start = Instant::now();
    let p = c.project_digonal(
        ProjectionBins::Linear(Unit::KeV(0.1)),
        Unit::KeV(2.),
        Unit::KeV(10.),
        false
    )?;
    let duration = start.elapsed();

    println!("Projection duration: {:?}", duration);

    // Write projection to file for inspection

    let mut f = std::fs::File::create("projection.csv").unwrap();

    for elem in p.spectrum.iter() {
        write!(f, "{elem}\n").unwrap();
    }

    println!("{:?}", p.ecal);

    Ok(())
}
