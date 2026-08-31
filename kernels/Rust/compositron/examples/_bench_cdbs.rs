use std::time::Instant;

use compositron::core::Measurement;
use compositron::core::utils::{Unit, EcalCorrectionOrder};
use compositron::cdbs::cdbspectrum::{
    Area, Axis, LineshapeParamDefinition, ProjectionBins
};
use compositron::importers::DataFormat;


fn print_mean_and_std(name: &str, data: &Vec<u128>) {
    let n = data.len();
    let mean = data.iter().sum::<u128>() as f64 / n as f64;
    let std = (
        data
            .iter()
            .map(|&x| (x as f64 - mean).powi(2))
            .sum::<f64>() / n as f64
    ).sqrt();

    println!("{name}: {mean} ns +/- {std} ns");
}


fn main() -> anyhow::Result<()> {
    let n_iter = 10000;

    let mut imports = Vec::with_capacity(n_iter);
    let mut corr = Vec::with_capacity(n_iter);
    let mut s_calc = Vec::with_capacity(n_iter);
    let mut proj_axes = Vec::with_capacity(n_iter);
    let mut proj_diag = Vec::with_capacity(n_iter);

    for _ in 0..n_iter {
        let start = Instant::now();
        let mut m = Measurement::from_file(
            "../../../../testdata/depth-profile_Copper_0000.n42",
            DataFormat::SlopeN42
        )?;
        imports.push(start.elapsed().as_nanos());

        let c = m.cdbs.get_mut("OAA x OAB").unwrap();

        let start = Instant::now();
        c.correct_ecal(EcalCorrectionOrder::First)?;
        corr.push(start.elapsed().as_nanos());

        let ls_param = LineshapeParamDefinition {
            name: "S",
            num: &[
                Area::Diagonal {
                    width_cel: Unit::KeV(2.),
                    width_cml: Unit::KeV(1.),
                },
            ],
            denom: &[
                Area::Diagonal {
                    width_cel: Unit::KeV(10.),
                    width_cml: Unit::KeV(1.),
                },
            ],
            in_corrected_peak: false,
        };

        let start = Instant::now();
        let s = c.calc_lineshape_param(&ls_param)?;
        s_calc.push(start.elapsed().as_nanos());

        let start = Instant::now();
        let p = c.project_digonal(
            ProjectionBins::Linear(Unit::KeV(0.1)),
            Unit::KeV(2.),
            Unit::KeV(100.),
            false
        )?;
        proj_diag.push(start.elapsed().as_nanos());

        let start = Instant::now();
        let p = c.project_axes(Axis::FirstDetector);
        proj_axes.push(start.elapsed().as_nanos());
    }

    print_mean_and_std("Import", &imports);
    print_mean_and_std("Ecal Corr", &corr);
    print_mean_and_std("S Param", &s_calc);
    print_mean_and_std("Proj Diag", &proj_diag);
    print_mean_and_std("Proj Axes", &proj_axes);

    Ok(())
}
