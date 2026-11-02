use compositron::{core::measurement_list::MeasurementList, importers::DataFormat};

fn main() -> anyhow::Result<()> {
    let dir = "/home/max/2025-04-10_CDB_spectra_W3Re_irr_18keV/";

    let mut mlist = MeasurementList::from_dir(dir, &DataFormat::SlopeN42)?;

    for m in mlist.measurements.iter_mut() {
        for d in m.dbs.values_mut() {
            d.default_analyze()?;
        }

        for c in m.cdbs.values_mut() {
            c.default_analyze()?;
        }
    }

    Ok(())
}
