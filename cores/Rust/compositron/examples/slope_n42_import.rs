use compositron::core::{Measurement, measurement::DataFormat};
use compositron::dbs::dbspectrum::EcalCorrectionOrder;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/depth-profile_Copper_0000.n42",
        DataFormat::SlopeN42
    )?;

    for (_, spec) in &mut m.dbs {
        spec.analyze(
            EcalCorrectionOrder::First,
            60.,
            true,
            compositron::dbs::dbspectrum::PeakModel::ErfLinear2Gauss,
            compositron::dbs::dbspectrum::STD_LINESHAPE_PARAMS,
            compositron::dbs::dbspectrum::BinIntegrationScheme::Const,
            false,
        )?;

        println!("{:?}", spec.peak_params);
    }

    println!("{}", m.shape());

    Ok(())
}
