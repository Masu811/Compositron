use std::io::Write;

use compositron::{elistcol_match, spectrum_match, spectrum2d_match};
use compositron::eventlist::eventlist::Bins;
use compositron::eventlist::coinc_eval::CoincEvalMode;
use compositron::importers::eventlist::nepomuc_cdbs_bin_importer::CDBSBinFile;


fn main() -> anyhow::Result<()> {
    let file = CDBSBinFile::from_file("/home/max/Raw_29889.bin")?;

    let mut cbins = Bins::from_bin_width(5., Some(10500.), Some(11500.))?;
    let mut rbins = cbins;

    let coinc_time_window = 300.;
    let coinc_eval_mode = CoincEvalMode::FirstPivot;

    let via_celist = true;

    let coinc_spectrum = if via_celist {
        let celist = file.eventlists.to_coinc_eventlist(
            "Det 0",
            "Det 2",
            "Det 0 x Det 2".into(),
            "timestamp",
            "timestamp",
            coinc_time_window,
            coinc_eval_mode,
        )?;

        celist.to_2d_histogram(
            Some("Det 0"),
            Some("Det 2"),
            "energy",
            "energy",
            &mut cbins,
            &mut rbins,
        )?
    } else {
        file.eventlists.to_2d_coinc_histogram(
            "Det 0",
            "Det 2",
            "timestamp",
            "timestamp",
            "energy",
            "energy",
            &mut cbins,
            &mut rbins,
            coinc_time_window,
            coinc_eval_mode,
        )?
    };

    spectrum2d_match!(
        &coinc_spectrum,
        arr => println!("{}", arr.iter().map(|&x| x as u64).sum::<u64>())
    );

    let mut f = std::fs::File::create("coinc_spectrum.csv")?;

    spectrum2d_match!(
        coinc_spectrum,
        arr => {
            for row in arr.row_iter() {
                for x in row.iter() {
                    write!(f, "{x} ")?;
                }
                writeln!(f, "")?;
            }
        }
    );

//     for (det, elist) in file.eventlists.eventlists.iter() {
//         let mut bins = Bins::from_bin_width(10., None, None)?;
//
//         let hist = elist.to_histogram("energy", &mut bins)?;
//
//         let mut file = std::fs::File::create(format!("./espectrum_{}.csv", det))?;
//
//         spectrum_match!(
//             hist,
//             arr => {
//                 for x in arr.iter() {
//                     writeln!(file, "{x}")?;
//                 }
//             }
//         )
//     }

    Ok(())
}
