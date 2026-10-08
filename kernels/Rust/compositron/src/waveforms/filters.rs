use crate::{trace_match, waveforms::waveform::{Waveform, WaveformRef}};


pub trait FilterWaveform {
    fn filter(&self, waveform: &WaveformRef) -> bool;
}


pub struct HeightFilter {
    pub all_higher_than: Option<f64>,
    pub all_lower_than: Option<f64>,
    pub any_higher_than: Option<f64>,
    pub any_lower_than: Option<f64>,
}


impl HeightFilter {
    fn check_all_higher_than(&self, waveform: &WaveformRef) -> bool {
        trace_match!(
            &waveform.trace, trace => self.all_higher_than
                .map_or(true, |min| trace.iter().all(|x| *x as f64 >= min))
        )
    }

    fn check_all_lower_than(&self, waveform: &WaveformRef) -> bool {
        trace_match!(
            &waveform.trace, trace => self.all_lower_than
                .map_or(true, |min| trace.iter().all(|x| *x as f64 <= min))
        )
    }

    fn check_any_higher_than(&self, waveform: &WaveformRef) -> bool {
        trace_match!(
            &waveform.trace, trace => self.any_higher_than
                .map_or(true, |min| trace.iter().any(|x| *x as f64 >= min))
        )
    }

    fn check_any_lower_than(&self, waveform: &WaveformRef) -> bool {
        trace_match!(
            &waveform.trace, trace => self.any_lower_than
                .map_or(true, |min| trace.iter().any(|x| *x as f64 <= min))
        )
    }
}


impl FilterWaveform for HeightFilter {
    fn filter(&self, waveform: &WaveformRef) -> bool {
        self.check_all_higher_than(waveform)
            && self.check_all_lower_than(waveform)
            && self.check_any_higher_than(waveform)
            && self.check_any_lower_than(waveform)
    }
}


pub struct AreaFilter {
    pub min_area: Option<f64>,
    pub max_area: Option<f64>,
}


pub struct HeightAreaRatioFilter {
    pub min_height_to_area_ratio: Option<f64>,
    pub max_height_to_area_ratio: Option<f64>,
}


pub struct IntegralShapeCheckpoint {
    pub x: f64,
    pub y_min: f64,
    pub y_max: f64,
}


pub struct IntegralShapeFilter {
    pub checkpoints: Vec<IntegralShapeCheckpoint>,
}


pub struct RiseAndDecayTimeFilter {
    pub mean_rise_time: f64,
    pub mean_decay_time: f64,
    pub radius_rise_time: f64,
    pub radius_decay_time: f64,
}
