use std::sync::Arc;

use thiserror::Error;

use crate::core::utils::LinearCalibration;
use crate::eventlist::eventlist::{EventList, EventListColumn, EventListColumnType};
use crate::waveforms::analysis;
use crate::waveforms::filters::FilterWaveform;
use crate::waveforms::transformations::TransformWaveform;


#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("Provided trace has 0 length")]
    EmptyTrace,

    #[error("Attempted access on index {idx}, but array has size {size}")]
    IndexError {
        idx: isize,
        size: usize,
    },

    #[error(
        "Cannot convert units of {from_unit} to LSB due to missing {missing_factor}"
    )]
    MissingConversionFactor {
        from_unit: String,
        missing_factor: String,
    },

    #[error("Pulse must be present at this point, which has not been checked")]
    NoPulseCheck,

    #[error("No pulse detected")]
    NoPulseDetected,

    #[error(
        "Baseline must have been calculated at this point, \
        which has not been done"
    )]
    NoBaseline,

    #[error("Pulse appears to be lying outside of trace")]
    PulseOutsideOfTrace,

    #[error("Fit result produced nan or inf values")]
    FitError,

    #[error(
        "Error exporting Eventlist column {missing_column}: \
        Parameter has not been calculated yet"
    )]
    MissingEventListColumn {
        missing_column: String,
    },
}


#[derive(Debug, Clone)]
pub enum Trace {
    I8(Arc<[i8]>),
    I16(Arc<[i16]>),
    I32(Arc<[i32]>),
    I64(Arc<[i64]>),
    F32(Arc<[f32]>),
    F64(Arc<[f64]>),
}


#[macro_export]
macro_rules! trace_match {
    ($s:expr, $v:ident => $body:expr) => {
        match $s {
            $crate::waveforms::waveform::Trace::I8($v) => $body,
            $crate::waveforms::waveform::Trace::I16($v) => $body,
            $crate::waveforms::waveform::Trace::I32($v) => $body,
            $crate::waveforms::waveform::Trace::I64($v) => $body,
            $crate::waveforms::waveform::Trace::F32($v) => $body,
            $crate::waveforms::waveform::Trace::F64($v) => $body,
        }
    };
}


#[macro_export]
macro_rules! split_trace_match {
    ($s:expr, $v:ident => $body_v:expr, $w:ident => $body_w:expr) => {
        match $s {
            $crate::waveforms::waveform::Trace::I8($v) => $body_v,
            $crate::waveforms::waveform::Trace::I16($v) => $body_v,
            $crate::waveforms::waveform::Trace::I32($v) => $body_v,
            $crate::waveforms::waveform::Trace::I64($v) => $body_v,
            $crate::waveforms::waveform::Trace::F32($w) => $body_w,
            $crate::waveforms::waveform::Trace::F64($w) => $body_w,
        }
    };
}


impl Trace {
    pub fn len(&self) -> usize {
        trace_match!(self, trace => trace.len())
    }
}


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PulsePolarity {
    Positive,
    Negative,
}


#[derive(Debug, Clone, Copy)]
pub enum Amplitude {
    Volts(f64),
    LSB(f64),
    CFLevel(f64),
    BaselineStdDev(f64),
}


#[derive(Debug, Clone)]
pub enum AmplitudeArray {
    Volts(Arc<[f64]>),
    LSB(Arc<[f64]>),
    CFLevel(Arc<[f64]>),
    BaselineStdDev(Arc<[f64]>),
}


impl Amplitude {
    pub fn to_lsb(
        &self,
        vcal: LinearCalibration,
        pulse_polarity: PulsePolarity,
        pulse_height: Option<f64>,
        baseline: Option<Average>,
    ) -> Result<f64, AnalysisError> {
        Ok(match self {
            Amplitude::Volts(x) => vcal.to_index_f64(*x),
            Amplitude::LSB(x) => *x,
            Amplitude::CFLevel(x) => {
                let Some(baseline) = baseline else {
                    return Err(AnalysisError::MissingConversionFactor {
                        from_unit: "CF Level".into(),
                        missing_factor: "baseline".into()
                    });
                };
                let Some(height) = pulse_height else {
                    return Err(AnalysisError::MissingConversionFactor {
                        from_unit: "CF Level".into(),
                        missing_factor: "pulse height".into()
                    });
                };

                baseline.mean + x * (height - baseline.mean)
            },
            Amplitude::BaselineStdDev(x) => {
                let Some(baseline) = baseline else {
                    return Err(AnalysisError::MissingConversionFactor {
                        from_unit: "baseline standard deviations".into(),
                        missing_factor: "baseline".into()
                    });
                };

                let sign = match pulse_polarity {
                    PulsePolarity::Positive => 1.,
                    PulsePolarity::Negative => -1.,
                };

                baseline.mean + x * sign * baseline.stddev
            },
        })
    }
}


impl AmplitudeArray {
    pub fn get(&self, index: usize) -> Amplitude {
        match self {
            AmplitudeArray::Volts(items) => Amplitude::Volts(
                items[index]
            ),
            AmplitudeArray::LSB(items) => Amplitude::LSB(
                items[index]
            ),
            AmplitudeArray::CFLevel(items) => Amplitude::CFLevel(
                items[index]
            ),
            AmplitudeArray::BaselineStdDev(items) => Amplitude::BaselineStdDev(
                items[index]
            ),
        }
    }

    pub fn slice(&self, idcs: &[usize]) -> AmplitudeArray {
        match self {
            AmplitudeArray::Volts(values) => AmplitudeArray::Volts(
                idcs.iter().map(|&i| values[i]).collect()
            ),
            AmplitudeArray::LSB(values) => AmplitudeArray::LSB(
                idcs.iter().map(|&i| values[i]).collect()
            ),
            AmplitudeArray::CFLevel(values) => AmplitudeArray::CFLevel(
                idcs.iter().map(|&i| values[i]).collect()
            ),
            AmplitudeArray::BaselineStdDev(values) => AmplitudeArray::BaselineStdDev(
                idcs.iter().map(|&i| values[i]).collect()
            ),
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub struct Average {
    pub mean: f64,
    pub stddev: f64,
}


#[derive(Debug, Clone)]
pub struct AverageArray {
    pub means: Arc<[f64]>,
    pub stddevs: Arc<[f64]>,
}


#[derive(Debug, Clone, Copy)]
pub struct Element {
    pub idx: usize,
    pub value: f64,
}


#[derive(Debug, Clone)]
pub struct ElementArray {
    pub idcs: Arc<[usize]>,
    pub values: Arc<[f64]>,
}


#[derive(Debug, Clone, Copy)]
pub enum TimeStep {
    Time(f64),
    Channels(usize),
}


impl TimeStep {
    pub fn to_channels(&self, time_channel_disp: f64) -> usize {
        match self {
            TimeStep::Time(t) => (t / time_channel_disp) as usize,
            TimeStep::Channels(i) => *i,
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub enum TimeRange {
    Time(f64, f64),
    Channels(usize, usize),
}


impl TimeRange {
    pub fn to_channels(&self, time_channel_disp: f64) -> (usize, usize) {
        match self {
            TimeRange::Time(t1, t2) => (
                (t1 / time_channel_disp) as usize, (t2 / time_channel_disp) as usize
            ),
            TimeRange::Channels(i1, i2) => (*i1, *i2),
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub enum CrossThresholdFrom {
    Below,
    Above,
    Any,
}


#[derive(Debug, Clone, Copy)]
pub enum Trigger {
    Simple {
        trigger_threshold: Amplitude,
        crossing_direction: CrossThresholdFrom,
    },
    Runlength {
        trigger_threshold: Amplitude,
        crossing_direction: CrossThresholdFrom,
        arm_runlength: usize,
        trigger_runlength: usize,
    },
}


#[derive(Debug, Clone)]
pub enum TriggerArray {
    Simple {
        trigger_thresholds: AmplitudeArray,
        crossing_direction: CrossThresholdFrom,
    },
    Runlength {
        trigger_thresholds: AmplitudeArray,
        crossing_direction: CrossThresholdFrom,
        arm_runlength: usize,
        trigger_runlength: usize,
    },
}


impl Trigger {
    fn simple_trigger(
        trace: &Trace,
        trigger_threshold: f64,
        crossing_direction: CrossThresholdFrom,
    ) -> Option<(usize, f64)> {
        trace_match!(
            trace, arr => {
                let mut armed = false;

                match crossing_direction {
                    CrossThresholdFrom::Below => {
                        for (i, &x) in arr.iter().enumerate() {
                            let x = x as f64;
                            if armed {
                                if x > trigger_threshold {
                                    return Some((i, trigger_threshold));
                                }
                            } else if x < trigger_threshold {
                                armed = true;
                            }
                        }
                    },
                    CrossThresholdFrom::Above => {
                        for (i, &x) in arr.iter().enumerate() {
                            let x = x as f64;
                            if armed {
                                if x < trigger_threshold {
                                    return Some((i, trigger_threshold));
                                }
                            } else if x > trigger_threshold {
                                armed = true;
                            }
                        }
                    },
                    CrossThresholdFrom::Any => return Trigger::simple_trigger(
                        trace,
                        trigger_threshold,
                        if (arr[0] as f64) < trigger_threshold {
                            CrossThresholdFrom::Below
                        } else {
                            CrossThresholdFrom::Above
                        }
                    ),
                }
            }
        );

        None
    }


    fn runlength_trigger(
        trace: &Trace,
        trigger_threshold: f64,
        crossing_direction: CrossThresholdFrom,
        arm_runlength: usize,
        trigger_runlength: usize,
    ) -> Option<(usize, f64)> {
        trace_match!(
            trace, arr => {
                let mut armed = false;
                let mut arm_count = 0;
                let mut trigger_count = 0;

                match crossing_direction {
                    CrossThresholdFrom::Below => {
                        for (i, &x) in arr.iter().enumerate() {
                            let x = x as f64;
                            if armed {
                                if x > trigger_threshold {
                                    if trigger_count >= trigger_runlength {
                                        return Some((i, trigger_threshold));
                                    } else {
                                        trigger_count += 1;
                                    }
                                } else {
                                    trigger_count = 0;
                                }
                            } else if x < trigger_threshold {
                                if arm_count >= arm_runlength {
                                    armed = true;
                                } else {
                                    arm_count += 1;
                                }
                            } else {
                                arm_count = 0;
                            }
                        }
                    },
                    CrossThresholdFrom::Above => {
                        for (i, &x) in arr.iter().enumerate() {
                            let x = x as f64;
                            if armed {
                                if x < trigger_threshold {
                                    if trigger_count >= trigger_runlength {
                                        return Some((i, trigger_threshold));
                                    } else {
                                        trigger_count += 1;
                                    }
                                } else {
                                    trigger_count = 0;
                                }
                            } else if x > trigger_threshold {
                                if arm_count >= arm_runlength {
                                    armed = true;
                                } else {
                                    arm_count += 1;
                                }
                            } else {
                                arm_count = 0;
                            }
                        }
                    },
                    CrossThresholdFrom::Any => return Trigger::runlength_trigger(
                        trace,
                        trigger_threshold,
                        if (arr[0] as f64) < trigger_threshold {
                            CrossThresholdFrom::Below
                        } else {
                            CrossThresholdFrom::Above
                        },
                        arm_runlength,
                        trigger_runlength,
                    ),
                }
            }
        );

        None
    }


    pub fn get_first_trigger_element(
        self,
        waveform: &Waveform,
    ) -> Result<Option<(usize, f64)>, AnalysisError> {
        match self {
            Trigger::Simple { trigger_threshold, crossing_direction } => {
                Ok(Trigger::simple_trigger(
                    &waveform.trace,
                    trigger_threshold.to_lsb(
                        waveform.vcal,
                        waveform.pulse_polarity,
                        waveform.params.pulse_height,
                        waveform.params.baseline,
                    )?,
                    crossing_direction
                ))
            },
            Trigger::Runlength {
                trigger_threshold, crossing_direction, arm_runlength,
                trigger_runlength,
            } => {
                Ok(Trigger::runlength_trigger(
                    &waveform.trace,
                    trigger_threshold.to_lsb(
                        waveform.vcal,
                        waveform.pulse_polarity,
                        waveform.params.pulse_height,
                        waveform.params.baseline,
                    )?,
                    crossing_direction,
                    arm_runlength,
                    trigger_runlength,
                ))
            },
        }
    }
}


impl TriggerArray {
    pub fn get(&self, index: usize) -> Trigger {
        match self {
            TriggerArray::Simple {
                trigger_thresholds, crossing_direction
            } => Trigger::Simple {
                trigger_threshold: trigger_thresholds.get(index),
                crossing_direction: *crossing_direction,
            },
            TriggerArray::Runlength {
                trigger_thresholds, crossing_direction,
                arm_runlength, trigger_runlength
            } => Trigger::Runlength {
                trigger_threshold: trigger_thresholds.get(index),
                crossing_direction: *crossing_direction,
                arm_runlength: *arm_runlength,
                trigger_runlength: *trigger_runlength,
            },
        }
    }


    pub fn slice(&self, idcs: &[usize]) -> TriggerArray {
        match self {
            TriggerArray::Simple {
                trigger_thresholds, crossing_direction
            } => TriggerArray::Simple {
                trigger_thresholds: trigger_thresholds.slice(idcs),
                crossing_direction: *crossing_direction,
            },
            TriggerArray::Runlength {
                trigger_thresholds, crossing_direction,
                arm_runlength, trigger_runlength,
            } => TriggerArray::Runlength {
                trigger_thresholds: trigger_thresholds.slice(idcs),
                crossing_direction: *crossing_direction,
                arm_runlength: *arm_runlength,
                trigger_runlength: *trigger_runlength,
            },
        }
    }
}


#[derive(Debug, Clone, Copy)]
pub struct TriggerTime {
    pub value: f64,
    pub threshold: f64,
}


#[derive(Debug, Clone)]
pub struct TriggerTimeArray {
    pub values: Arc<[f64]>,
    pub thresholds: Arc<[f64]>,
}


#[derive(Debug, Clone, Copy)]
pub struct TriggerTimeRange {
    pub value: f64,
    pub start_trigger: Trigger,
    pub stop_trigger: Trigger,
}


#[derive(Debug, Clone)]
pub struct TriggerTimeRangeArray {
    pub values: Arc<[f64]>,
    pub start_triggers: TriggerArray,
    pub stop_triggers: TriggerArray,
}


#[derive(Debug, Clone, Copy)]
pub enum PulseHeightDeterminationMethod {
    ExtremalElement,
    QuadraticInterpolation,
    CubicFit { n_points: usize },
}


#[derive(Debug, Clone, Copy)]
pub enum TriggerTimeDeterminationMethod {
    FirstThresholdCrossingElement,
    LinearInterpolation,
    QuarticInterpolation,
    LinearFit { n_points: usize },
    QuadraticFit { n_points: usize },
}


#[derive(Debug, Clone, Copy)]
pub struct WaveformParams {
    pub contains_pulse: Option<bool>,
    pub baseline: Option<Average>,
    pub extremal_element: Option<Element>,
    pub pulse_height: Option<f64>,
    pub pulse_time: Option<f64>,
    pub pulse_area: Option<f64>,
    pub trigger_time: Option<TriggerTime>,
    pub pulse_rise_time: Option<TriggerTimeRange>,
    pub pulse_decay_time: Option<TriggerTimeRange>,
}


#[derive(Debug, Clone)]
pub struct WaveformParamArray {
    pub contains_pulse: Option<Arc<[bool]>>,
    pub baselines: Option<AverageArray>,
    pub extremal_elements: Option<ElementArray>,
    pub pulse_heights: Option<Arc<[f64]>>,
    pub pulse_times: Option<Arc<[f64]>>,
    pub pulse_areas: Option<Arc<[f64]>>,
    pub trigger_times: Option<TriggerTimeArray>,
    pub pulse_rise_times: Option<TriggerTimeRangeArray>,
    pub pulse_decay_times: Option<TriggerTimeRangeArray>,
}


#[derive(Debug, Clone)]
pub struct Waveform {
    pub trace: Trace,
    pub vcal: LinearCalibration,
    pub time_channel_disp: f64,
    pub global_time_offset: u64,
    pub pulse_polarity: PulsePolarity,
    pub params: WaveformParams,
}


#[derive(Debug, Clone)]
pub struct WaveformList {
    pub data: Arc<[Waveform]>,
    active_indices: Vec<usize>,
}


#[derive(Debug, Clone)]
pub struct WaveformArray {
    pub traces: Arc<[Trace]>,
    pub time_channel_disp: f64,
    pub global_time_offsets: Arc<[u64]>,
    pub vcal: LinearCalibration,
    pub pulse_polarity: PulsePolarity,
    pub params: WaveformParamArray,
    active_indices: Vec<usize>,
}


impl Waveform {
    pub fn new(
        trace: Trace,
        vcal: LinearCalibration,
        time_channel_disp: f64,
        global_time_offset: u64,
        pulse_polarity: PulsePolarity,
    ) -> Self {
        Waveform {
            trace,
            vcal,
            time_channel_disp,
            global_time_offset,
            pulse_polarity,
            params: WaveformParams {
                contains_pulse: None,
                baseline: None,
                extremal_element: None,
                pulse_height: None,
                pulse_time: None,
                pulse_area: None,
                pulse_rise_time: None,
                pulse_decay_time: None,
                trigger_time: None,
            },
        }
    }


    pub fn average_range(
        &mut self, time_range: TimeRange
    ) -> Result<Average, AnalysisError> {
        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let n = self.trace.len();

        if i1 >= n {
            return Err(AnalysisError::IndexError { idx: i1 as isize, size: n });
        }

        if i2 >= n {
            return Err(AnalysisError::IndexError { idx: i2 as isize, size: n });
        }

        let integration_range = (i2 - i1) as f64;

        Ok(analysis::average_range((i1, i2), integration_range, &self.trace))
    }


    pub fn get_baseline(
        &mut self, time_range: TimeRange
    ) -> Result<Average, AnalysisError> {
        let baseline = self.average_range(time_range)?;

        self.params.baseline = Some(baseline);

        Ok(baseline)
    }


    pub fn get_extremal_element(
        &mut self
    ) -> Result<Element, AnalysisError> {
        if let Some(extremal_element) = self.params.extremal_element {
            return Ok(extremal_element);
        }

        let Some(element) = analysis::get_extremal_element(
            &self.trace, self.pulse_polarity
        ) else { return Err(AnalysisError::EmptyTrace); };

        let elem = Element { idx: element.idx, value: element.value };

        self.params.extremal_element = Some(elem);

        Ok(elem)
    }


    pub fn check_for_pulse(
        &mut self, min_required_pulse_height: Amplitude,
    ) -> Result<bool, AnalysisError> {
        let Some(baseline) = self.params.baseline else {
            return Err(AnalysisError::NoBaseline);
        };

        let extr_elem = self.get_extremal_element()?;

        let extr = extr_elem.value;

        let threshold = min_required_pulse_height.to_lsb(
            self.vcal,
            self.pulse_polarity,
            self.params.pulse_height,
            self.params.baseline,
        )?;

        let contains_pulse = analysis::check_for_pulse(
            baseline.mean, extr, threshold
        );

        self.params.contains_pulse = Some(contains_pulse);

        Ok(contains_pulse)
    }


    pub fn get_pulse_height_and_time(
        &mut self,
        pulse_height_determination_method: PulseHeightDeterminationMethod,
    ) -> Result<(f64, f64), AnalysisError> {
        let Some(contains_pulse) = self.params.contains_pulse else {
            return Err(AnalysisError::NoPulseCheck);
        };

        if !contains_pulse {
            return Err(AnalysisError::NoPulseDetected);
        }

        let Some(baseline) = self.params.baseline else {
            return Err(AnalysisError::NoBaseline);
        };

        let extr_elem = self.get_extremal_element()?;

        let (extr_idx, extr) = analysis::interpolate_pulse_height_and_time(
            &self.trace,
            extr_elem.idx,
            extr_elem.value,
            pulse_height_determination_method,
        )?;

        let pulse_height = (baseline.mean - extr).abs();

        self.params.pulse_height = Some(pulse_height);
        self.params.pulse_time = Some(extr_idx);

        Ok((extr_idx, pulse_height))
    }


    pub fn get_area(
        &mut self,
        time_range: TimeRange,
    ) -> Result<f64, AnalysisError> {
        let Some(baseline) = self.params.baseline else {
            return Err(AnalysisError::NoBaseline);
        };

        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let integration_range = (i2 - i1) as f64;

        let area = analysis::get_area(
            &self.trace, i1, i2, integration_range, baseline.mean,
        )?;

        Ok(area)
    }


    pub fn get_pulse_area(
        &mut self,
        integration_range_before_pulse_max: TimeStep,
        integration_range_after_pulse_max: TimeStep,
    ) -> Result<f64, AnalysisError> {
        let Some(contains_pulse) = self.params.contains_pulse else {
            return Err(AnalysisError::NoPulseCheck);
        };

        if !contains_pulse {
            return Err(AnalysisError::NoPulseDetected);
        }

        let Some(baseline) = self.params.baseline else {
            return Err(AnalysisError::NoBaseline);
        };

        let extr_elem = self.get_extremal_element()?;

        let extr_idx = extr_elem.idx;

        let i1 = integration_range_before_pulse_max.to_channels(self.time_channel_disp);
        let i2 = integration_range_after_pulse_max.to_channels(self.time_channel_disp);

        let n = self.trace.len();

        if i1 > extr_idx {
            return Err(AnalysisError::IndexError {
                idx: extr_idx as isize - i1 as isize, size: n
            });
        }

        if i2 + extr_idx >= n {
            return Err(AnalysisError::IndexError {
                idx: (i2 + extr_idx) as isize, size: n
            });
        }

        let i1 = extr_idx - i1;
        let i2 = extr_idx + i2;

        let integration_range = (i2 + i1) as f64;

        let area = analysis::get_area(
            &self.trace, i1, i2, integration_range, baseline.mean,
        )?;

        self.params.pulse_area = Some(area);

        Ok(area)
    }


    pub fn get_trigger_time(
        &mut self,
        trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<f64, AnalysisError> {
        let Some(contains_pulse) = self.params.contains_pulse else {
            return Err(AnalysisError::NoPulseCheck);
        };

        if !contains_pulse {
            return Err(AnalysisError::NoPulseDetected);
        }

        let Some((trig_idx, trig)) =
            trigger.get_first_trigger_element(&self)? else {
                return Err(AnalysisError::NoPulseDetected);
            };

        let trigger_time = analysis::interpolate_trigger_time(
            &self.trace,
            trig_idx,
            trig,
            trigger_time_determination_method,
        )?;

        self.params.trigger_time = Some(TriggerTime {
            value: trigger_time, threshold: trig
        });

        Ok(trigger_time)
    }


    pub fn get_trigger_time_range(
        &mut self,
        start_trigger: Trigger,
        stop_trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<f64, AnalysisError> {
        let Some(contains_pulse) = self.params.contains_pulse else {
            return Err(AnalysisError::NoPulseCheck);
        };

        if !contains_pulse {
            return Err(AnalysisError::NoPulseDetected);
        }

        let Some((start_trig_idx, start_trig)) =
            start_trigger.get_first_trigger_element(&self)? else {
                return Err(AnalysisError::NoPulseDetected);
            };

        let start_time = analysis::interpolate_trigger_time(
            &self.trace,
            start_trig_idx,
            start_trig,
            trigger_time_determination_method,
        )?;

        let Some((stop_trig_idx, stop_trig)) =
            stop_trigger.get_first_trigger_element(&self)? else {
                return Err(AnalysisError::NoPulseDetected);
            };

        let stop_time = analysis::interpolate_trigger_time(
            &self.trace,
            stop_trig_idx,
            stop_trig,
            trigger_time_determination_method,
        )?;

        let time_range = stop_time - start_time;

        Ok(time_range)
    }
}


impl WaveformList {
    pub fn new(data: Vec<Waveform>) -> Self {
        let l = data.len();

        WaveformList {
            data: Arc::from(data),
            active_indices: (0..l).collect(),
        }
    }


    pub fn transform(
        &self, transform: impl TransformWaveform,
    ) -> Result<WaveformList, AnalysisError> {
        todo!();
    }


    pub fn filter(
        &self, filter: impl FilterWaveform,
    ) -> Result<WaveformList, AnalysisError> {
        Ok(WaveformList {
            data: self.data.clone(),
            active_indices: self.active_indices
                .iter()
                .copied()
                .filter(|&i| {
                    let waveform = &self.data[i];

                    filter.filter(&Waveform {
                        trace: waveform.trace.clone(),
                        vcal: waveform.vcal,
                        time_channel_disp: waveform.time_channel_disp,
                        global_time_offset: waveform.global_time_offset,
                        pulse_polarity: waveform.pulse_polarity,
                        params: waveform.params,
                    })
                })
                .collect()
        })
    }


    pub fn to_eventlist(
        &self,
        eventlist: &mut EventList,
        pre_filter_transforms: Option<&[impl TransformWaveform]>,
        filters: Option<&[impl FilterWaveform]>,
    ) {
        todo!();
    }
}


impl WaveformArray {
    pub fn new(
        traces: Vec<Trace>,
        time_channel_disp: f64,
        global_time_offsets: Vec<u64>,
        vcal: LinearCalibration,
        pulse_polarity: PulsePolarity,
    ) -> Self {
        let l = traces.len();
        WaveformArray {
            traces: Arc::from(traces),
            vcal,
            time_channel_disp,
            global_time_offsets: Arc::from(global_time_offsets),
            pulse_polarity,
            params: WaveformParamArray {
                contains_pulse: None,
                baselines: None,
                extremal_elements: None,
                pulse_heights: None,
                pulse_areas: None,
                pulse_times: None,
                pulse_rise_times: None,
                pulse_decay_times: None,
                trigger_times: None,
            },
            active_indices: (0..l).collect(),
        }
    }


    pub fn get(&self, index: usize) -> Waveform {
        let i = self.active_indices[index];
        Waveform {
            trace: self.traces[i].clone(),
            vcal: self.vcal,
            time_channel_disp: self.time_channel_disp,
            global_time_offset: self.global_time_offsets[i],
            pulse_polarity: self.pulse_polarity,
            params: WaveformParams {
                contains_pulse:
                    self.params.contains_pulse.as_deref().map(|arr| arr[index]),
                baseline:
                    self.params.baselines.as_ref().map(|avgs| Average {
                        mean: avgs.means[index], stddev: avgs.stddevs[index],
                    }),
                extremal_element:
                    self.params.extremal_elements.as_ref().map(|elems| Element {
                        idx: elems.idcs[index], value: elems.values[index],
                    }),
                pulse_height:
                    self.params.pulse_heights.as_deref().map(|arr| arr[index]),
                pulse_time:
                    self.params.pulse_times.as_deref().map(|arr| arr[index]),
                pulse_area:
                    self.params.pulse_areas.as_deref().map(|arr| arr[index]),
                trigger_time:
                    self.params.trigger_times.as_ref().map(|trig_times| TriggerTime {
                        value: trig_times.values[index],
                        threshold: trig_times.thresholds[index],
                    }),
                pulse_rise_time:
                    self.params.pulse_rise_times.as_ref().map(|rise_times| TriggerTimeRange {
                        value: rise_times.values[index],
                        start_trigger: rise_times.start_triggers.get(index),
                        stop_trigger: rise_times.stop_triggers.get(index),
                    }),
                pulse_decay_time:
                self.params.pulse_decay_times.as_ref().map(|decay_times| TriggerTimeRange {
                    value: decay_times.values[index],
                    start_trigger: decay_times.start_triggers.get(index),
                    stop_trigger: decay_times.stop_triggers.get(index),
                }),
            },
        }
    }


    pub fn average_range(
        &mut self, time_range: TimeRange
    ) -> Result<AverageArray, AnalysisError> {
        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let n = self.traces[0].len();

        if i1 >= n {
            return Err(AnalysisError::IndexError { idx: i1 as isize, size: n });
        }

        if i2 >= n {
            return Err(AnalysisError::IndexError { idx: i2 as isize, size: n });
        }

        let integration_range = (i2 - i1) as f64;

        let mut means = Vec::with_capacity(self.active_indices.len());
        let mut stddevs = Vec::with_capacity(self.active_indices.len());

        for &i in &self.active_indices {
            let avg = analysis::average_range((i1, i2), integration_range, &self.traces[i]);

            means.push(avg.mean);
            stddevs.push(avg.stddev);
        }

        Ok(AverageArray { means: Arc::from(means), stddevs: Arc::from(stddevs) })
    }


    pub fn get_baseline(
        &mut self, time_range: TimeRange
    ) -> Result<(), AnalysisError> {
        let baselines = self.average_range(time_range)?;

        self.params.baselines = Some(baselines);

        Ok(())
    }


    pub fn get_extremal_element(
        &mut self
    ) -> Result<(), AnalysisError> {
        if self.params.extremal_elements.is_some() {
            return Ok(());
        }

        let mut indices = Vec::with_capacity(self.active_indices.len());
        let mut values = Vec::with_capacity(self.active_indices.len());

        for &i in &self.active_indices {
            let Some(element) = analysis::get_extremal_element(
                &self.traces[i], self.pulse_polarity
            ) else { return Err(AnalysisError::EmptyTrace); };

            indices.push(element.idx);
            values.push(element.value);
        }

        self.params.extremal_elements = Some(ElementArray {
            idcs: Arc::from(indices), values: Arc::from(values),
        });

        Ok(())
    }


    pub fn check_for_pulse(
        &mut self, min_required_pulse_height: Amplitude,
    ) -> Result<(), AnalysisError> {
        if self.params.baselines.is_none() {
            return Err(AnalysisError::NoBaseline);
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let baselines = self.params.baselines.as_ref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_ref().unwrap();

        let mut contains_pulse = Vec::with_capacity(self.active_indices.len());

        for i in 0..self.active_indices.len() {
            let threshold = min_required_pulse_height.to_lsb(
                self.vcal,
                self.pulse_polarity,
                self.params.pulse_heights.as_ref().map(|arr| arr[i]),
                self.params.baselines.as_ref().map(|arr| Average {
                    mean: arr.means[i], stddev: arr.stddevs[i],
                }),
            )?;

            contains_pulse.push(analysis::check_for_pulse(
                baselines.means[i], extremal_elements.values[i], threshold
            ));
        }

        self.params.contains_pulse = Some(Arc::from(contains_pulse));

        Ok(())
    }


    pub fn get_pulse_height_and_time(
        &mut self,
        pulse_height_determination_method: PulseHeightDeterminationMethod,
    ) -> Result<(), AnalysisError> {
        if self.params.contains_pulse.is_none() {
            return Err(AnalysisError::NoPulseCheck);
        }

        if self.params.baselines.is_none() {
            return Err(AnalysisError::NoBaseline);
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let contains_pulse = self.params.contains_pulse.as_deref().unwrap();
        let baselines = self.params.baselines.as_ref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_ref().unwrap();

        let mut pulse_heights = Vec::with_capacity(self.active_indices.len());
        let mut pulse_times = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_pulse[i] {
                pulse_heights.push(f64::NAN);
                pulse_times.push(f64::NAN);
                continue;
            }

            let baseline = baselines.means[i];

            let extr_idx = extremal_elements.idcs[i];
            let extr = extremal_elements.values[i];

            let (extr_idx, extr) = analysis::interpolate_pulse_height_and_time(
                &self.traces[idx],
                extr_idx,
                extr,
                pulse_height_determination_method,
            )?;

            let pulse_height = (baseline- extr).abs();

            pulse_heights.push(pulse_height);
            pulse_times.push(self.time_channel_disp * extr_idx);
        }

        self.params.pulse_heights = Some(Arc::from(pulse_heights));
        self.params.pulse_times = Some(Arc::from(pulse_times));

        Ok(())
    }


    pub fn get_area(
        &mut self,
        time_range: TimeRange,
    ) -> Result<(), AnalysisError> {
        if self.params.baselines.is_none() {
            return Err(AnalysisError::NoBaseline);
        }

        let baselines = self.params.baselines.as_ref().unwrap();

        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let integration_range = (i2 - i1) as f64;

        let mut areas = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            areas.push(analysis::get_area(
                &self.traces[idx], i1, i2, integration_range, baselines.means[i],
            )?);
        }

        Ok(())
    }


    pub fn get_pulse_area(
        &mut self,
        integration_range_before_pulse_max: TimeStep,
        integration_range_after_pulse_max: TimeStep,
    ) -> Result<(), AnalysisError> {
        if self.params.contains_pulse.is_none() {
            return Err(AnalysisError::NoPulseCheck);
        }

        if self.params.baselines.is_none() {
            return Err(AnalysisError::NoBaseline);
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let contains_pulse = self.params.contains_pulse.as_deref().unwrap();
        let baselines = self.params.baselines.as_ref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_ref().unwrap();

        let i1 = integration_range_before_pulse_max.to_channels(self.time_channel_disp);
        let i2 = integration_range_after_pulse_max.to_channels(self.time_channel_disp);

        let mut areas = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_pulse[i] {
                areas.push(f64::NAN);
                continue;
            }

            let extr_idx = extremal_elements.idcs[i];

            let n = self.traces[idx].len();

            if i1 > extr_idx {
                return Err(AnalysisError::IndexError {
                    idx: extr_idx as isize - i1 as isize, size: n
                });
            }

            if i2 + extr_idx >= n {
                return Err(AnalysisError::IndexError {
                    idx: (i2 + extr_idx) as isize, size: n
                });
            }

            let i1 = extr_idx - i1;
            let i2 = extr_idx + i2;

            let integration_range = (i2 + i1) as f64;

            areas.push(analysis::get_area(
                &self.traces[idx], i1, i2, integration_range, baselines.means[i],
            )?);
        }

        self.params.pulse_areas = Some(Arc::from(areas));

        Ok(())
    }


    pub fn get_trigger_time(
        &mut self,
        trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<(), AnalysisError> {
        if self.params.contains_pulse.is_none() {
            return Err(AnalysisError::NoPulseCheck);
        }

        let contains_pulse = self.params.contains_pulse.as_deref().unwrap();

        let mut trig_values = Vec::with_capacity(self.active_indices.len());
        let mut trig_thresholds = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_pulse[i] {
                trig_values.push(f64::NAN);
                trig_thresholds.push(f64::NAN);
                continue;
            }

            let Some((trig_idx, trig)) =
                trigger.get_first_trigger_element(&self.get(i))? else {
                    trig_values.push(f64::NAN);
                    trig_thresholds.push(f64::NAN);
                    continue;
            };

            let Ok(trigger_time) = analysis::interpolate_trigger_time(
                &self.traces[idx],
                trig_idx,
                trig,
                trigger_time_determination_method,
            ) else {
                trig_values.push(f64::NAN);
                trig_thresholds.push(f64::NAN);
                continue;
            };

            trig_values.push(trigger_time);
            trig_thresholds.push(trig);
        }

        self.params.trigger_times = Some(TriggerTimeArray {
            values: Arc::from(trig_values),
            thresholds: Arc::from(trig_thresholds),
        });

        Ok(())
    }


    pub fn get_trigger_time_range(
        &mut self,
        start_trigger: Trigger,
        stop_trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<Vec<f64>, AnalysisError> {
        if self.params.contains_pulse.is_none() {
            return Err(AnalysisError::NoPulseCheck);
        }

        let contains_pulse = self.params.contains_pulse.as_deref().unwrap();

        let mut time_ranges = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_pulse[i] {
                time_ranges.push(f64::NAN);
                continue;
            }

            let Some((start_trig_idx, start_trig)) =
                start_trigger.get_first_trigger_element(&self.get(i))? else {
                time_ranges.push(f64::NAN);
                continue;
            };

            let Ok(start_trigger_time) = analysis::interpolate_trigger_time(
                &self.traces[idx],
                start_trig_idx,
                start_trig,
                trigger_time_determination_method,
            ) else {
                time_ranges.push(f64::NAN);
                continue;
            };

            let Some((stop_trig_idx, stop_trig)) =
                stop_trigger.get_first_trigger_element(&self.get(i))? else {
                time_ranges.push(f64::NAN);
                continue;
            };

            let Ok(stop_trigger_time) = analysis::interpolate_trigger_time(
                &self.traces[idx],
                stop_trig_idx,
                stop_trig,
                trigger_time_determination_method,
            ) else {
                time_ranges.push(f64::NAN);
                continue;
            };

            time_ranges.push(stop_trigger_time - start_trigger_time);
        }

        Ok(time_ranges)
    }


    pub fn transform(&self, transforms: impl TransformWaveform) {
        todo!();
    }


    pub fn filter(&self, filter: impl FilterWaveform) -> WaveformArray {
        let active_indices = self.active_indices
            .iter()
            .filter(|&&i| filter.filter(&self.get(i)))
            .map(|&i| i)
            .collect::<Vec<_>>();

        WaveformArray {
            traces: self.traces.clone(),
            vcal: self.vcal,
            time_channel_disp: self.time_channel_disp,
            global_time_offsets: self.global_time_offsets.clone(),
            pulse_polarity: self.pulse_polarity,
            params: WaveformParamArray {
                contains_pulse: self.params.contains_pulse.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                baselines: self.params.baselines.as_ref().and_then(
                    |arr| Some(AverageArray {
                        means: active_indices.iter().map(|&i| arr.means[i]).collect(),
                        stddevs: active_indices.iter().map(|&i| arr.stddevs[i]).collect(),
                    })
                ),
                extremal_elements: self.params.extremal_elements.as_ref().and_then(
                    |arr| Some(ElementArray {
                        idcs: active_indices.iter().map(|&i| arr.idcs[i]).collect(),
                        values: active_indices.iter().map(|&i| arr.values[i]).collect(),
                    })
                ),
                pulse_heights: self.params.pulse_heights.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                pulse_areas: self.params.pulse_areas.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                pulse_times: self.params.pulse_times.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                trigger_times: self.params.trigger_times.as_ref().and_then(
                    |arr| Some(TriggerTimeArray {
                        values: active_indices.iter().map(|&i| arr.values[i]).collect(),
                        thresholds: active_indices.iter().map(|&i| arr.thresholds[i]).collect(),
                    })
                ),
                pulse_rise_times: self.params.pulse_rise_times.as_ref().and_then(
                    |arr| Some(TriggerTimeRangeArray {
                        values: active_indices.iter().map(|&i| arr.values[i]).collect(),
                        start_triggers: arr.start_triggers.slice(&active_indices),
                        stop_triggers: arr.stop_triggers.slice(&active_indices),
                    })
                ),
                pulse_decay_times: self.params.pulse_decay_times.as_ref().and_then(
                    |arr| Some(TriggerTimeRangeArray {
                        values: active_indices.iter().map(|&i| arr.values[i]).collect(),
                        start_triggers: arr.start_triggers.slice(&active_indices),
                        stop_triggers: arr.stop_triggers.slice(&active_indices),
                    })
                ),
            },
            active_indices,
        }
    }


    pub fn to_eventlist(
        &self, columns: &[EventListColumnType],
    ) -> Result<EventList, AnalysisError> {
        let mut eventlist = EventList {
            columns: Vec::new(),
        };

        for col in columns {
            match col {
                EventListColumnType::Timestamp => {
                    let Some(trigger_times) = &self.params.trigger_times else {
                        return Err(AnalysisError::MissingEventListColumn {
                            missing_column: "trigger_times".into()
                        });
                    };

                    eventlist.columns.push(
                        EventListColumn::Timestamp(
                            Arc::from(self.global_time_offsets
                                .iter()
                                .zip(trigger_times.values.iter())
                                .map(|(&t_g, &t_l)| t_g + t_l as u64)
                                .collect::<Vec<_>>()
                            )
                        )
                    )
                },
                EventListColumnType::PulseHeight => {
                    let Some(pulse_heights) = &self.params.pulse_heights else {
                        return Err(AnalysisError::MissingEventListColumn {
                            missing_column: "pulse_heights".into()
                        });
                    };

                    eventlist.columns.push(
                        EventListColumn::PulseHeight(pulse_heights.clone())
                    );
                },
                EventListColumnType::PulseArea => {
                    let Some(pulse_areas) = &self.params.pulse_areas else {
                        return Err(AnalysisError::MissingEventListColumn {
                            missing_column: "pulse_areas".into()
                        });
                    };

                    eventlist.columns.push(
                        EventListColumn::PulseHeight(pulse_areas.clone())
                    );
                },
                EventListColumnType::PulseRiseTime => {
                    let Some(pulse_rise_times) = &self.params.pulse_rise_times else {
                        return Err(AnalysisError::MissingEventListColumn {
                            missing_column: "pulse_rise_times".into()
                        });
                    };

                    eventlist.columns.push(
                        EventListColumn::PulseHeight(
                            pulse_rise_times.values.clone()
                        )
                    );
                },
                EventListColumnType::PulseDecayTime => {
                    let Some(pulse_decay_times) = &self.params.pulse_decay_times else {
                        return Err(AnalysisError::MissingEventListColumn {
                            missing_column: "pulse_decay_times".into()
                        });
                    };

                    eventlist.columns.push(
                        EventListColumn::PulseHeight(
                            pulse_decay_times.values.clone()
                        )
                    );
                },
            }
        }

        Ok(eventlist)
    }
}
