use std::time::Instant;

use compositron::core::Measurement;
use compositron::importers::DataFormat;
use compositron::core::utils::{EcalCorrectionOrder, Unit};
use compositron::dbs::dbspectrum::{Area, LineshapeParamDefinition, PeakModel};


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
    let n_iter = 1000;

    let mut imports = Vec::with_capacity(n_iter);
    let mut corr = Vec::with_capacity(n_iter);
    let mut bg_sub = Vec::with_capacity(n_iter);
    let mut s_calc = Vec::with_capacity(n_iter);

    for _ in 0..n_iter {
        let start = Instant::now();
        let mut m = Measurement::from_file(
            "../../../../testdata/depth-profile_Copper_0000.n42",
            DataFormat::SlopeN42
        )?;
        imports.push(start.elapsed().as_nanos());

        let d = m.dbs.get_mut("OAB").unwrap();

        let start = Instant::now();
        d.correct_ecal(EcalCorrectionOrder::First)?;
        corr.push(start.elapsed().as_nanos());

        let start = Instant::now();
        d.extract_peak(30., PeakModel::ErfLinear2Gauss)?;
        bg_sub.push(start.elapsed().as_nanos());

        let param = LineshapeParamDefinition {
            name: "S",
            num: &[Area::Peak { width: Unit::KeV(2.) }],
            denom: &[Area::Peak { width: Unit::KeV(20.) }],
            in_corrected_peak: true,
        };

        let start = Instant::now();
        d.calc_lineshape_param(&param, compositron::dbs::dbspectrum::BinIntegrationScheme::Const)?;
        s_calc.push(start.elapsed().as_nanos());
    }

    print_mean_and_std("Import", &imports);
    print_mean_and_std("Ecal Corr", &corr);
    print_mean_and_std("BG Sub", &bg_sub);
    print_mean_and_std("S Param", &s_calc);

    Ok(())
}
