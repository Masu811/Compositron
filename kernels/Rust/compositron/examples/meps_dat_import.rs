use compositron::core::Measurement;
use compositron::importers::DataFormat;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat",
        DataFormat::MePSDat,
    )?;

    println!("{}", m.shape());

    Ok(())
}
