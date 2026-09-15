use compositron::{core::measurement_array::MeasurementArray, importers::DataFormat};

fn main() -> anyhow::Result<()> {
    let dir = "/home/max/2026-09-05_depth-profile_Li_Reference/";

    let marr = MeasurementArray::from_dir(dir, &DataFormat::SlopeN42)?;

    Ok(())
}
