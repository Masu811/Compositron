use compositron::core::Measurement;
use compositron::importers::DataFormat;


fn main() -> anyhow::Result<()> {
    let mut m = Measurement::from_file(
        "../../../../testdata/depth-profile_Copper_0000.n42",
        DataFormat::SlopeN42
    )?;

    let d = m.dbs.get_mut("OAB").unwrap();

    d.default_analyze()?;

    let s = d.lineshape_params.get("S").unwrap();
    let w = d.lineshape_params.get("W").unwrap();
    let v = d.lineshape_params.get("V/P").unwrap();
    let p = d.lineshape_params.get("P/T").unwrap();

    println!("  S   = {:.6} +/- {:.6}", s.val, s.err);
    println!("  W   = {:.6} +/- {:.6}", w.val, w.err);
    println!("  V/P = {:.6} +/- {:.6}", v.val, v.err);
    println!("  P/T = {:.6} +/- {:.6}", p.val, p.err);

    Ok(())
}
