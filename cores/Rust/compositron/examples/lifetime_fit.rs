use compositron::core::Measurement;
use compositron::lifetime_model;
use compositron::pals::fitting::fit_lifetime_spectrum;

fn main() -> anyhow::Result<()> {
    let mut m = Measurement::new();

    m.import(
        "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat",
        compositron::core::measurement::DataFormat::MePSDat
    )?;

    let lt_model = lifetime_model!(
        "N": (1e6, 0),
        "background": (0, false),
        "t0": (500),
        "lifetime_1": (180, 100, 5000),
        "intensity_1": (0.79, 0, 1),
        "lifetime_2": (350, 100, 5000),
        "intensity_2": (0.2, 0, 1),
        "lifetime_3": (2700, 100, 5000),
        "intensity_3": (0.01, 0, 1),
        "res_fwhm_1": (90, 30, 500),
        "res_intensity_1": (0.8, 0, 1),
        "res_t0_1": (0, false),
        "res_fwhm_2": (126, 30, 500),
        "res_intensity_2": (0.2, 0, 1),
        "res_t0_2": (0, -500, 500)
    )?;

    let s = m.pals.get("A").unwrap();

    let y = compositron::spectrum_match!(
        &s.spectrum, arr => arr[16500..17000].iter().map(|&x| x as f64).collect::<Vec<f64>>()
    );

    let x = (0..y.len()).map(|i| i as f64 * s.tcal.1 + s.tcal.0).collect::<Vec<f64>>();

    fit_lifetime_spectrum(&x, &y, lt_model)?;

    Ok(())
}
