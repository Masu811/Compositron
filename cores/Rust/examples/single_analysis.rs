use cstacs::importer::n42_importer::import_n42;

fn main() -> anyhow::Result<()> {
    let mut m = import_n42("../../../testdata/depth-profile_Copper_0000.n42")?;

    println!("{}", m.shape());

    let s = m.singles.get_mut("OAA").unwrap();

    s.analyze(
        1.,
        1.,
        1.,
        true,
        30.,
        true,
        (400., 500., 506., 516.),
        cstacs::single::FollowPeakOrder::First,
    )?;

    println!("S: {} +/- {}", s.s(), s.ds());
    println!("W: {} +/- {}", s.w(), s.dw());
    println!("V2P: {} +/- {}", s.v2p(), s.dv2p());

    Ok(())
}
