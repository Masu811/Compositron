use compositron::core::Measurement;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    )?;

    println!("{}", m.shape());

    Ok(())
}
