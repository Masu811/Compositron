use compositron::core::{Measurement, measurement::DataFormat};

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/depth-profile_Copper_0000.n42",
        DataFormat::SlopeN42
    )?;

    println!("{}", m.shape());

    for (det, spec) in &mut m.dbs {
        spec.default_analyze()?;

        let s = spec.lineshape_params.get("S").unwrap();
        let w = spec.lineshape_params.get("W").unwrap();
        let v = spec.lineshape_params.get("V/P").unwrap();
        let p = spec.lineshape_params.get("P/T").unwrap();

        println!("{det}:");
        println!("  S   = {:.3} +/- {:.3}", s.val, s.err);
        println!("  W   = {:.3} +/- {:.3}", w.val, w.err);
        println!("  V/P = {:.3} +/- {:.3}", v.val, v.err);
        println!("  P/T = {:.3} +/- {:.3}", p.val, p.err);
    }

    Ok(())
}
