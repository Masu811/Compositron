use std::collections::HashMap;
use std::io::Write;

use compositron::elistcol_match;
use compositron::importers::waveforms::lecroy_trc_importer;
use compositron::waveforms::filters::PulsePresenceFilter;
use compositron::waveforms::waveform::{
    Amplitude, CrossThresholdFrom, PulseHeightDeterminationMethod, PulsePolarity, TimeRange, TimeStep, Trigger, TriggerTimeDeterminationMethod, WaveformParam
};


fn main() -> anyhow::Result<()> {
    let mut warr = lecroy_trc_importer::import_lecroy_trc(
        "../../../../testdata/waveforms/C1--Kapton--00000.trc",
        "C1".into(),
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

    let elist = warr.to_eventlist(HashMap::from([
        (WaveformParam::TriggerTime, "timestamp".into()),
        (WaveformParam::PulseHeight, "energy".into()),
    ]))?;

    let mut f = std::fs::File::create("eventlist.csv")?;

    let output_columns = vec!["timestamp", "energy"];

    let output_columns = output_columns
        .iter()
        .map(|&key| elist.columns.get(key).unwrap())
        .collect::<Vec<_>>();

    for row in 0..elist.nrows() {
        for col in &output_columns {
            elistcol_match!(
                col,
                arr => write!(f, "{},", arr[row])?
            );
        }
        writeln!(f, "")?;
    }

    Ok(())
}
