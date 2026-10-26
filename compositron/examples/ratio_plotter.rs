use compositron::{cdbs::cdbspectrum::{Axis, ProjectionBins}, core::{measurement_list::MeasurementList, utils::Unit}, importers::DataFormat};
use std::io::Write;


fn main() -> anyhow::Result<()> {
    let mut mlist = MeasurementList::from_dir(
        "/home/max/2026-09-15_CDB_spectra_Li_Reference_PostIssues_newCal_25kv",
        &DataFormat::SlopeN42
    )?;

    let mut p_li_rod = Vec::new();
    let mut p_li_foil = Vec::new();
    let mut p_cu = Vec::new();

    for m in mlist.measurements.iter_mut() {
        let c = m.cdbs.get_mut("OAA x OAB").unwrap();

        c.correct_ecal(compositron::core::utils::EcalCorrectionOrder::First)?;

        let p = c.project(
            Axis::CEL,
            ProjectionBins::Linear(Unit::KeV(0.1)),
            Unit::Eres(1.),
            Unit::KeV(10.),
            false
        )?;

        let p = p.fold()?;

        match m.metadata.get("Sample").unwrap().as_str() {
            "Li_foil" => {
                if p_li_foil.len() == 0 {
                    for x in p.spectrum.iter() {
                        p_li_foil.push(*x);
                    }
                } else {
                    for i in 0..p.spectrum.len() {
                        p_li_foil[i] += p.spectrum[i];
                    }
                }
            },
            "Li_rod" => {
                if p_li_rod.len() == 0 {
                    for x in p.spectrum.iter() {
                        p_li_rod.push(*x);
                    }
                } else {
                    for i in 0..p.spectrum.len() {
                        p_li_rod[i] += p.spectrum[i];
                    }
                }
            },
            "Cu_Foil" => {
                if p_cu.len() == 0 {
                    for x in p.spectrum.iter() {
                        p_cu.push(*x);
                    }
                } else {
                    for i in 0..p.spectrum.len() {
                        p_cu[i] += p.spectrum[i];
                    }
                }
            },
            _ => {
                return Err(anyhow::Error::msg("Unexpected Sample name"));
            }
        }
    }

    let mut f = std::fs::File::create("./projection_li_foil.csv")?;

    for x in p_li_foil.iter() {
        writeln!(f, "{x}").unwrap();
    }

    let mut f = std::fs::File::create("./projection_li_rod.csv")?;

    for x in p_li_rod.iter() {
        writeln!(f, "{x}").unwrap();
    }

    let mut f = std::fs::File::create("./projection_cu.csv")?;

    for x in p_cu.iter() {
        writeln!(f, "{x}").unwrap();
    }

    Ok(())
}
