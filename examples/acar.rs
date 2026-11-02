use std::io::Write;
use compositron::{core::utils::Detector, elistcol_match, eventlist::eventlist::{Bins, EventListColumnType}, importers::eventlist::bin_importer::import_bin_eventlist_row_major, spectrum_match, spectrum2d_match};

fn main() -> anyhow::Result<()> {
    let columns: Vec<(String, EventListColumnType)> = vec![
        ("x1".into(), EventListColumnType::I32),
        ("y1".into(), EventListColumnType::I32),
        ("e1".into(), EventListColumnType::I32),
        ("x2".into(), EventListColumnType::I32),
        ("y2".into(), EventListColumnType::I32),
        ("e2".into(), EventListColumnType::I32),
    ];

    let keys = columns.iter().map(|x| x.0.clone()).collect::<Vec<_>>();

    let elist = import_bin_eventlist_row_major(
        "/home/max/measurements/old_2d_acar/Cu-poly03-annealed-14K-minus60A-2020-10-09.coi",
        columns,
        "A".into(),
        compositron::core::utils::ByteOrder::LittleEndian,
        0
    )?;

    let mut bins_x = Bins::from_bin_width(1000., None, None)?;
    let mut bins_y = bins_x;

    let hist = elist.to_2d_histogram("x2", "y2", &mut bins_x, &mut bins_y)?;

    println!("{:?}", bins_x);

    let mut f = std::fs::File::create("acar_2d_positions.csv")?;

    spectrum2d_match!(
        hist,
        arr => {
            for row in arr.row_iter() {
                for x in row.iter() {
                    write!(f, "{x} ")?;
                }
                writeln!(f, "")?;
            }
        }
    );

    let mut bins = Bins::from_bin_width(1000., None, None)?;

    let hist = elist.to_histogram("x1", &mut bins)?;

    let mut f = std::fs::File::create("acar_espectrum.csv")?;

    spectrum_match!(
        hist,
        arr => {
            for x in arr.iter() {
                writeln!(f, "{x}")?;
            }
        }
    );

    Ok(())
}
