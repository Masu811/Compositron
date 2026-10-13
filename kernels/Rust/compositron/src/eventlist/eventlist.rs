use std::sync::Arc;

use thiserror::Error;

use crate::cdbs::{CDBSpectrum, cdbspectrum::Orientation};
use crate::dbs::DBSpectrum;
use crate::pals::PALSpectrum;


#[derive(Debug, Error)]
pub enum AnalysisError {

}


pub enum EventListColumnType {
    Timestamp,
    PulseHeight,
    PulseArea,
    PulseRiseTime,
    PulseDecayTime,
}


pub enum EventListColumn {
    Timestamp(Arc<[u64]>),
    PulseHeight(Arc<[f64]>),
    PulseArea(Arc<[f64]>),
    PulseRiseTime(Arc<[f64]>),
    PulseDecayTime(Arc<[f64]>),
}


pub struct Bins {
    pub bin_width: f64,
    pub min_bin_edge: Option<f64>,
    pub max_bin_edge: Option<f64>,
}


pub struct EventList {
    pub columns: Vec<EventListColumn>,
}


impl EventList {
    pub fn new() -> Self {
        EventList { columns: Vec::new() }
    }

    pub fn add_column(&mut self, column: EventListColumn) {
        self.columns.push(column);
    }

    pub fn to_histogram(
        &self,
        column_index: usize,
        bins: Bins,
    ) {
        todo!();
    }

    pub fn to_dbspectrum(&self) -> DBSpectrum {
        todo!();
    }
}


pub struct CoincidenceEventLists<'a> {
    pub events_first_det: &'a EventList,
    pub events_second_det: &'a EventList,
    pub timestamp_offset: f64,
}


impl<'a> CoincidenceEventLists<'a> {
    pub fn new() -> Self {
        todo!();
    }

    pub fn to_cdbspectrum(
        &self,
        coinc_time_window: f64,
        bins_first_det: Bins,
        bins_second_det: Bins,
        orientation: Orientation,
        reject_pileups: bool,
    ) -> CDBSpectrum {
        todo!();
    }

    pub fn to_palspectrum(
        &self,
        bins: Bins,
        start_energy_window_first_det: Option<(f64, f64)>,
        stop_energy_window_first_det: Option<(f64, f64)>,
        start_energy_window_second_det: Option<(f64, f64)>,
        stop_energy_window_second_det: Option<(f64, f64)>,
        reject_pileups: bool,
    ) -> PALSpectrum {
        todo!();
    }
}


pub struct TripleCoincidenceEventLists<'a> {
    pub events_first_det: &'a EventList,
    pub events_second_det: &'a EventList,
    pub events_third_det: &'a EventList,
    pub timestamp_offset_second_det: f64,
    pub timestamp_offset_third_det: f64,
}


impl<'a> TripleCoincidenceEventLists<'a> {
    pub fn new() -> Self {
        todo!();
    }

    pub fn to_amocspectrum(&self) {
        todo!();
    }
}
