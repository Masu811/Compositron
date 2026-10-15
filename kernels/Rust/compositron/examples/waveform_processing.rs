use std::io::Write;

use compositron::eventlist::eventlist::{EventListColumn, EventListColumnType};
use compositron::importers::waveforms::lecroy_trc_importer;
use compositron::waveforms::filters::PulsePresenceFilter;
use compositron::waveforms::waveform::{
    Amplitude, CrossThresholdFrom, PulseHeightDeterminationMethod,
    PulsePolarity, TimeRange, TimeStep, Trigger, TriggerTimeDeterminationMethod
};


fn main() -> anyhow::Result<()> {
    let mut warr = lecroy_trc_importer::import_lecroy_trc(
        "../../../../testdata/waveforms/C1--Kapton--00000.trc",
        PulsePolarity::Negative,
    )?;

    warr.get_baseline(TimeRange::Samples(0, 10))?;

    warr.check_for_pulse(Amplitude::LSB(5.))?;

    let mut warr = warr.filter(PulsePresenceFilter);

    warr.get_pulse_height_and_time(
        PulseHeightDeterminationMethod::CubicFit { n_points: 6 }
    )?;

    warr.get_pulse_area(TimeStep::Time(10000.), TimeStep::Time(10000.))?;

    warr.get_trigger_time(
        Trigger::Simple {
            trigger_threshold: Amplitude::CFLevel(0.2),
            crossing_direction: CrossThresholdFrom::Above,
        },
        TriggerTimeDeterminationMethod::LinearFit { n_points: 4 },
    )?;

    let elist = warr.to_eventlist(&[
        EventListColumnType::Timestamp,
        EventListColumnType::PulseHeight,
        EventListColumnType::PulseArea,
    ])?;

    let mut f = std::fs::File::create("eventlist.csv")?;

    let nrows = match &elist.columns[0] {
        EventListColumn::Timestamp(items) => items.len(),
        EventListColumn::PulseHeight(items) => items.len(),
        EventListColumn::PulseArea(items) => items.len(),
        EventListColumn::PulseRiseTime(items) => items.len(),
        EventListColumn::PulseDecayTime(items) => items.len(),
    };

    for row in 0..nrows {
        for col in elist.columns.iter() {
            match col {
                EventListColumn::Timestamp(items) => {
                    write!(f, "{},", items[row])?;
                },
                EventListColumn::PulseHeight(items) => {
                    write!(f, "{},", items[row])?;
                },
                EventListColumn::PulseArea(items) => {
                    write!(f, "{},", items[row])?;
                },
                EventListColumn::PulseRiseTime(items) => {
                    write!(f, "{},", items[row])?;
                },
                EventListColumn::PulseDecayTime(items) => {
                    write!(f, "{},", items[row])?;
                },
            };
        }
        writeln!(f, "")?;
    }

    Ok(())
}
