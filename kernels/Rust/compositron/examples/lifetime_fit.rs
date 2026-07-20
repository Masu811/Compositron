use compositron::core::Measurement;
use compositron::lifetime_model;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat",
        compositron::core::measurement::DataFormat::MePSDat
    )?;

    let lt_model = lifetime_model!(
        "lifetime_1": (180, 100, 5000),
        "intensity_1": (0.79, 0, 1),
        "lifetime_2": (350, 100, 5000),
        "intensity_2": (0.2, 0, 1),
        "lifetime_3": (2700, 100, 5000),
        "intensity_3": (0.01, 0, 1),
        "res_fwhm_1": (210, 50, 1000),
        "res_intensity_1": (0.8, 0, 1),
        "res_fwhm_2": (300, 50, 1000),
        "res_intensity_2": (0.2, 0, 1),
        "res_t0_2": (0, -200, 200)
    )?;

    let s = m.pals.get_mut("A").unwrap();

    s.fit(lt_model, 16500, 18000)?;

    print!("{}", s.fit_report()?);

    Ok(())
}
