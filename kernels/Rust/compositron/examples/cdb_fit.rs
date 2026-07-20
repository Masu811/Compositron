use compositron::core::Measurement;
use compositron::cdbs::cdbspectrum::EcalCorrectionOrder;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/depth-profile_Copper_0000.n42",
        compositron::core::measurement::DataFormat::SlopeN42
    )?;

    let s = m.cdbs.get_mut("OAA x OAB").unwrap();

    println!("{:?}", s.ecal);

    s.correct_ecal(EcalCorrectionOrder::Zeroth)?;

    println!("{:?}", s.corrected_ecal);

    Ok(())
}
