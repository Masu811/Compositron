use std::sync::Arc;

use thiserror::Error;

use crate::core::utils::LinearCalibration;
use crate::eventlist::eventlist::EventList;
use crate::waveforms::analysis;
use crate::waveforms::filters::FilterWaveform;
use crate::waveforms::transformations::TransformWaveform;


#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("Provided trace has 0 length")]
    EmptyTrace,

    #[error("Provided range exceeds size of trace")]
    SizeError,

    #[error(
        "Cannot convert units of {from_unit} to LSB due to missing {missing_factor}"
    )]
    MissingConversionFactor {
        from_unit: String,
        missing_factor: String,
    },

    #[error("Peak must be present at this point, which has not been checked")]
    NoPeakCheck,

    #[error("No peak detected")]
    NoPeakDetected,

    #[error("Peak appears to be lying outside of the trace")]
    PeakOutsideOfTrace,

    #[error("Fit error")]
    FitError,
}


#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PeakPolarity {
    Positive,
    Negative,
}


#[derive(Debug, Clone)]
pub enum Trace {
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
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


#[derive(Debug, Clone, Copy)]
pub struct Average {
    pub mean: f64,
    pub stddev: f64,
}


#[derive(Debug, Clone, Copy)]
pub struct TriggerTimeRange {
    pub value: f64,
    pub start_trigger: Trigger,
    pub stop_trigger: Trigger,
}


#[derive(Debug, Clone, Copy)]
pub struct TriggerTime {
    pub value: f64,
    pub threshold_lsb: f64,
}


#[derive(Debug, Clone, Copy)]
pub enum PeakHeightDeterminationMethod {
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
pub enum Amplitude {
    Volts(f64),
    LSB(f64),
    CFLevel(f64),
    BaselineStdDev(f64),
}


impl Amplitude {
    pub fn to_lsb(
        &self,
        vcal: LinearCalibration,
        peak_polarity: PeakPolarity,
        peak_height: Option<f64>,
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
                let Some(height) = peak_height else {
                    return Err(AnalysisError::MissingConversionFactor {
                        from_unit: "CF Level".into(),
                        missing_factor: "peak height".into()
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

                let sign = match peak_polarity {
                    PeakPolarity::Positive => 1.,
                    PeakPolarity::Negative => -1.,
                };

                baseline.mean + x * sign * baseline.stddev
            },
        })
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
        waveform: WaveformRef,
    ) -> Result<Option<(usize, f64)>, AnalysisError> {
        match self {
            Trigger::Simple { trigger_threshold, crossing_direction } => {
                Ok(Trigger::simple_trigger(
                    &waveform.trace,
                    trigger_threshold.to_lsb(
                        waveform.vcal,
                        waveform.peak_polarity,
                        waveform.params.peak_height,
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
                        waveform.peak_polarity,
                        waveform.params.peak_height,
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


#[derive(Debug, Clone, Copy)]
pub struct WaveformParams {
    pub contains_peak: Option<bool>,
    pub baseline: Option<Average>,
    pub extremal_element: Option<(usize, f64)>,
    pub peak_height: Option<f64>,
    pub peak_time: Option<f64>,
    pub peak_area: Option<f64>,
    pub trigger_time: Option<TriggerTime>,
    pub peak_rise_time: Option<TriggerTimeRange>,
    pub peak_decay_time: Option<TriggerTimeRange>,
}


#[derive(Debug, Clone)]
pub struct Waveform {
    pub trace: Trace,
    pub vcal: LinearCalibration,
    pub time_channel_disp: f64,
    pub global_time_offset: f64,
    pub peak_polarity: PeakPolarity,
    pub params: WaveformParams,
}


#[derive(Debug, Clone)]
pub struct WaveformRef<'a> {
    pub trace: &'a Trace,
    pub vcal: LinearCalibration,
    pub time_channel_disp: f64,
    pub global_time_offset: f64,
    pub peak_polarity: PeakPolarity,
    pub params: WaveformParams,
}


impl Waveform {
    pub fn new(
        trace: Trace,
        vcal: LinearCalibration,
        time_channel_disp: f64,
        global_time_offset: f64,
        peak_polarity: PeakPolarity,
    ) -> Self {
        Waveform {
            trace,
            vcal,
            time_channel_disp,
            global_time_offset,
            peak_polarity,
            params: WaveformParams {
                contains_peak: None,
                baseline: None,
                extremal_element: None,
                peak_height: None,
                peak_time: None,
                peak_area: None,
                peak_rise_time: None,
                peak_decay_time: None,
                trigger_time: None,
            },
        }
    }


    pub fn as_waveform_ref(&self) -> WaveformRef {
        WaveformRef {
            trace: &self.trace,
            vcal: self.vcal,
            time_channel_disp: self.time_channel_disp,
            global_time_offset: self.global_time_offset,
            peak_polarity: self.peak_polarity,
            params: self.params,
        }
    }


    pub fn average_range(
        &mut self, time_range: TimeRange
    ) -> Result<Average, AnalysisError> {
        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let n = self.trace.len();

        if i1 < 0 || i1 >= n || i2 < 0 || i2 >= n || i2 > i1 {
            return Err(AnalysisError::SizeError);
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
    ) -> Result<(usize, f64), AnalysisError> {
        if let Some(extremal_element) = self.params.extremal_element {
            return Ok(extremal_element);
        }

        let Some((extr_idx, extr)) = analysis::get_extremal_element(
            &self.trace, self.peak_polarity
        ) else { return Err(AnalysisError::EmptyTrace); };

        self.params.extremal_element = Some((extr_idx, extr));

        Ok((extr_idx, extr))
    }


    pub fn check_for_peak(
        &mut self, min_required_peak_height: Amplitude,
    ) -> Result<bool, AnalysisError> {
        let baseline = self.params.baseline.unwrap_or(
            self.get_baseline(TimeRange::Channels(0, 10))?
        );

        let (_, extr) = self.get_extremal_element()?;

        let threshold = min_required_peak_height.to_lsb(
            self.vcal,
            self.peak_polarity,
            self.params.peak_height,
            self.params.baseline,
        )?;

        let contains_peak = analysis::check_for_peak(baseline.mean, extr, threshold);

        self.params.contains_peak = Some(contains_peak);

        Ok(contains_peak)
    }


    pub fn get_peak_height_and_time(
        &mut self,
        peak_height_determination_method: PeakHeightDeterminationMethod,
    ) -> Result<(f64, f64), AnalysisError> {
        let contains_peak = self.params.contains_peak.unwrap_or(
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?
        );

        if !contains_peak {
            return Err(AnalysisError::NoPeakDetected);
        }

        let baseline = self.params.baseline.unwrap_or(
            self.get_baseline(TimeRange::Channels(0, 10))?
        );

        let (extr_idx, extr) = self.get_extremal_element()?;

        let (extr_idx, extr) = analysis::interpolate_peak_height_and_time(
            &self.trace, extr_idx, extr, peak_height_determination_method
        )?;

        let peak_height = (baseline.mean - extr).abs();

        self.params.peak_height = Some(peak_height);
        self.params.peak_time = Some(extr_idx);

        Ok((extr_idx, peak_height))
    }


    pub fn get_area(
        &mut self,
        time_range: TimeRange,
    ) -> Result<f64, AnalysisError> {
        let baseline = self.params.baseline.unwrap_or(
            self.get_baseline(TimeRange::Channels(0, 10))?
        );

        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        if i2 > i1 {
            return Err(AnalysisError::SizeError);
        }

        let integration_range = (i2 - i1) as f64;

        let area = analysis::get_area(
            &self.trace, i1, i2, integration_range, baseline.mean,
        )?;

        Ok(area)
    }


    pub fn get_peak_area(
        &mut self,
        integration_range_before_peak_max: TimeStep,
        integration_range_after_peak_max: TimeStep,
    ) -> Result<f64, AnalysisError> {
        let contains_peak = self.params.contains_peak.unwrap_or(
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?
        );

        if !contains_peak {
            return Err(AnalysisError::NoPeakDetected);
        }

        let baseline = self.params.baseline.unwrap_or(
            self.get_baseline(TimeRange::Channels(0, 10))?
        );

        let (extr_idx, _) = self.get_extremal_element()?;

        let i1 = integration_range_before_peak_max.to_channels(self.time_channel_disp);
        let i2 = integration_range_after_peak_max.to_channels(self.time_channel_disp);

        if i1 > extr_idx || i2 + extr_idx >= self.trace.len() {
            return Err(AnalysisError::SizeError);
        }

        let i1 = extr_idx - i1;
        let i2 = extr_idx + i2;

        let integration_range = (i2 + i1) as f64;

        let area = analysis::get_area(
            &self.trace, i1, i2, integration_range, baseline.mean,
        )?;

        self.params.peak_area = Some(area);

        Ok(area)
    }


    pub fn interpolate_trigger_time(
        &mut self,
        trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<f64, AnalysisError> {
        let contains_peak = self.params.contains_peak.unwrap_or(
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?
        );

        if !contains_peak {
            return Err(AnalysisError::NoPeakDetected);
        }

        let Some((trig_idx, trig)) =
            trigger.get_first_trigger_element(self.as_waveform_ref())? else {
            return Err(AnalysisError::NoPeakDetected);
        };

        let trigger_time = analysis::interpolate_trigger_time(
            &self.trace,
            trig_idx,
            trig,
            trigger_time_determination_method,
        )?;

        self.params.trigger_time = Some(TriggerTime {
            value: trigger_time, threshold_lsb: trig
        });

        Ok(trigger_time)
    }


    pub fn interpolate_trigger_time_range(
        &mut self,
        start_trigger: Trigger,
        stop_trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<f64, AnalysisError> {
        let contains_peak = self.params.contains_peak.unwrap_or(
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?
        );

        if !contains_peak {
            return Err(AnalysisError::NoPeakDetected);
        }

        let Some((start_trig_idx, start_trig)) =
            start_trigger.get_first_trigger_element(self.as_waveform_ref())? else {
            return Err(AnalysisError::NoPeakDetected);
        };

        let start_time = analysis::interpolate_trigger_time(
            &self.trace,
            start_trig_idx,
            start_trig,
            trigger_time_determination_method,
        )?;

        let Some((stop_trig_idx, stop_trig)) =
            stop_trigger.get_first_trigger_element(self.as_waveform_ref())? else {
            return Err(AnalysisError::NoPeakDetected);
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


#[derive(Debug, Clone)]
pub struct WaveformList {
    data: Arc<[Waveform]>,
    active_indices: Vec<usize>,
}


impl WaveformList {
    pub fn new(data: Vec<Waveform>) -> Self {
        let l = data.len();

        WaveformList {
            data: Arc::from(data),
            active_indices: (0..l).collect(),
        }
    }


    fn get(&self, index: usize) -> WaveformRef {
        self.data[index].as_waveform_ref()
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

                    filter.filter(&WaveformRef {
                        trace: &waveform.trace,
                        vcal: waveform.vcal,
                        time_channel_disp: waveform.time_channel_disp,
                        global_time_offset: waveform.global_time_offset,
                        peak_polarity: waveform.peak_polarity,
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


#[derive(Debug, Clone)]
pub struct WaveformParamArray {
    pub contains_peak: Option<Vec<bool>>,
    pub baselines: Option<Vec<Average>>,
    pub extremal_elements: Option<Vec<(usize, f64)>>,
    pub peak_heights: Option<Vec<f64>>,
    pub peak_areas: Option<Vec<f64>>,
    pub peak_times: Option<Vec<f64>>,
    pub peak_rise_times: Option<Vec<TriggerTimeRange>>,
    pub peak_decay_times: Option<Vec<TriggerTimeRange>>,
    pub trigger_times: Option<Vec<TriggerTime>>,
}


#[derive(Debug, Clone)]
pub struct WaveformArray {
    traces: Arc<[Trace]>,
    time_channel_disp: f64,
    global_time_offsets: Arc<[f64]>,
    vcal: LinearCalibration,
    peak_polarity: PeakPolarity,
    params: WaveformParamArray,
    active_indices: Vec<usize>,
}


impl WaveformArray {
    pub fn new(
        traces: Vec<Trace>,
        time_channel_disp: f64,
        global_time_offsets: Vec<f64>,
        vcal: LinearCalibration,
        peak_polarity: PeakPolarity,
    ) -> Self {
        let l = traces.len();
        WaveformArray {
            traces: Arc::from(traces),
            vcal,
            time_channel_disp,
            global_time_offsets: Arc::from(global_time_offsets),
            peak_polarity,
            params: WaveformParamArray {
                contains_peak: None,
                baselines: None,
                extremal_elements: None,
                peak_heights: None,
                peak_areas: None,
                peak_times: None,
                peak_rise_times: None,
                peak_decay_times: None,
                trigger_times: None,
            },
            active_indices: (0..l).collect(),
        }
    }


    pub fn get(&self, index: usize) -> WaveformRef {
        let i = self.active_indices[index];
        WaveformRef {
            trace: &self.traces[i],
            vcal: self.vcal,
            time_channel_disp: self.time_channel_disp,
            global_time_offset: self.global_time_offsets[i],
            peak_polarity: self.peak_polarity,
            params: WaveformParams {
                contains_peak:
                    self.params.contains_peak.as_deref().map(|arr| arr[i]),
                baseline:
                    self.params.baselines.as_deref().map(|arr| arr[i]),
                extremal_element:
                    self.params.extremal_elements.as_deref().map(|arr| arr[i]),
                peak_height:
                    self.params.peak_heights.as_deref().map(|arr| arr[i]),
                peak_time:
                    self.params.peak_times.as_deref().map(|arr| arr[i]),
                peak_area:
                    self.params.peak_areas.as_deref().map(|arr| arr[i]),
                trigger_time:
                    self.params.trigger_times.as_deref().map(|arr| arr[i]),
                peak_rise_time:
                    self.params.peak_rise_times.as_deref().map(|arr| arr[i]),
                peak_decay_time:
                    self.params.peak_decay_times.as_deref().map(|arr| arr[i]),
            },
        }
    }


    pub fn average_range(
        &mut self, time_range: TimeRange
    ) -> Result<Vec<Average>, AnalysisError> {
        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        let n = self.traces[0].len();

        if i1 < 0 || i1 >= n || i2 < 0 || i2 >= n || i2 > i1 {
            return Err(AnalysisError::SizeError);
        }

        let integration_range = (i2 - i1) as f64;

        let mut averages = Vec::with_capacity(self.active_indices.len());

        for &i in &self.active_indices {
            averages.push(
                analysis::average_range((i1, i2), integration_range, &self.traces[i])
            );
        }

        Ok(averages)
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
        if !self.params.extremal_elements.is_none() {
            return Ok(());
        }

        let mut extremal_elements = Vec::with_capacity(self.active_indices.len());

        for &i in &self.active_indices {
            let Some(element) = analysis::get_extremal_element(
                &self.traces[i], self.peak_polarity
            ) else { return Err(AnalysisError::EmptyTrace); };

            extremal_elements.push(element);
        }

        self.params.extremal_elements = Some(extremal_elements);

        Ok(())
    }


    pub fn check_for_peak(
        &mut self, min_required_peak_height: Amplitude,
    ) -> Result<(), AnalysisError> {
        if self.params.baselines.is_none() {
            self.get_baseline(TimeRange::Channels(0, 10))?;
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let baselines = self.params.baselines.as_deref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_deref().unwrap();

        let mut contains_peak = Vec::with_capacity(self.active_indices.len());

        for i in 0..self.active_indices.len() {
            let threshold = min_required_peak_height.to_lsb(
                self.vcal,
                self.peak_polarity,
                self.params.peak_heights.as_ref().map(|arr| arr[i]),
                self.params.baselines.as_ref().map(|arr| arr[i]),
            )?;

            let (_, extr) = extremal_elements[i];

            contains_peak.push(analysis::check_for_peak(
                baselines[i].mean, extr, threshold
            ));
        }

        self.params.contains_peak = Some(contains_peak);

        Ok(())
    }


    pub fn get_peak_height_and_time(
        &mut self,
        peak_height_determination_method: PeakHeightDeterminationMethod,
    ) -> Result<(&[f64], &[f64]), AnalysisError> {
        if self.params.contains_peak.is_none() {
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?;
        }

        if self.params.baselines.is_none() {
            self.get_baseline(TimeRange::Channels(0, 10))?;
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let contains_peak = self.params.contains_peak.as_deref().unwrap();
        let baselines = self.params.baselines.as_deref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_deref().unwrap();

        let mut peak_heights = Vec::with_capacity(self.active_indices.len());
        let mut peak_times = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_peak[i] {
                peak_heights.push(f64::NAN);
                peak_times.push(f64::NAN);
                continue;
            }

            let baseline = baselines[i];

            let (extr_idx, extr) = extremal_elements[i];

            let (extr_idx, extr) = analysis::interpolate_peak_height_and_time(
                &self.traces[idx],
                extr_idx,
                extr,
                peak_height_determination_method,
            )?;

            let peak_height = (baseline.mean - extr).abs();

            peak_heights.push(peak_height);
            peak_times.push(self.time_channel_disp * extr_idx);
        }

        self.params.peak_heights = Some(peak_heights);
        self.params.peak_times = Some(peak_times);

        Ok((
            self.params.peak_heights.as_ref().unwrap(),
            self.params.peak_times.as_ref().unwrap(),
        ))
    }


    pub fn get_area(
        &mut self,
        time_range: TimeRange,
    ) -> Result<(), AnalysisError> {
        if self.params.baselines.is_none() {
            self.get_baseline(TimeRange::Channels(0, 10))?;
        }

        let baselines = self.params.baselines.as_deref().unwrap();

        let (i1, i2) = time_range.to_channels(self.time_channel_disp);

        if i2 > i1 {
            return Err(AnalysisError::SizeError);
        }

        let integration_range = (i2 - i1) as f64;

        let mut areas = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            areas.push(analysis::get_area(
                &self.traces[idx], i1, i2, integration_range, baselines[i].mean,
            )?);
        }

        Ok(())
    }


    pub fn get_peak_area(
        &mut self,
        integration_range_before_peak_max: TimeStep,
        integration_range_after_peak_max: TimeStep,
    ) -> Result<(), AnalysisError> {
        if self.params.contains_peak.is_none() {
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?;
        }

        if self.params.baselines.is_none() {
            self.get_baseline(TimeRange::Channels(0, 10))?;
        }

        if self.params.extremal_elements.is_none() {
            self.get_extremal_element()?;
        }

        let contains_peak = self.params.contains_peak.as_deref().unwrap();
        let baselines = self.params.baselines.as_deref().unwrap();
        let extremal_elements = self.params.extremal_elements.as_deref().unwrap();

        let i1 = integration_range_before_peak_max.to_channels(self.time_channel_disp);
        let i2 = integration_range_after_peak_max.to_channels(self.time_channel_disp);

        let mut areas = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_peak[i] {
                areas.push(f64::NAN);
                continue;
            }

            let (extr_idx, _) = extremal_elements[i];

            if i1 > extr_idx || i2 + extr_idx >= self.traces[idx].len() {
                return Err(AnalysisError::SizeError);
            }

            let i1 = extr_idx - i1;
            let i2 = extr_idx + i2;

            let integration_range = (i2 + i1) as f64;

            areas.push(analysis::get_area(
                &self.traces[idx], i1, i2, integration_range, baselines[i].mean,
            )?);
        }

        self.params.peak_areas = Some(areas);

        Ok(())
    }


    pub fn interpolate_trigger_time(
        &mut self,
        trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<(), AnalysisError> {
        if self.params.contains_peak.is_none() {
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?;
        }

        let contains_peak = self.params.contains_peak.as_deref().unwrap();

        let mut trigger_times = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_peak[i] {
                trigger_times.push(TriggerTime {
                    value: f64::NAN, threshold_lsb: f64::NAN,
                });
                continue;
            }

            let Some((trig_idx, trig)) =
                trigger.get_first_trigger_element(self.get(i))? else {
                trigger_times.push(TriggerTime {
                    value: f64::NAN, threshold_lsb: f64::NAN,
                });
                continue;
            };

            let trigger_time = analysis::interpolate_trigger_time(
                &self.traces[idx],
                trig_idx,
                trig,
                trigger_time_determination_method,
            )?;

            trigger_times.push(TriggerTime {
                value: trigger_time, threshold_lsb: trig
            });
        }

        self.params.trigger_times = Some(trigger_times);

        Ok(())
    }


    pub fn interpolate_trigger_time_range(
        &mut self,
        start_trigger: Trigger,
        stop_trigger: Trigger,
        trigger_time_determination_method: TriggerTimeDeterminationMethod,
    ) -> Result<Vec<f64>, AnalysisError> {
        if self.params.contains_peak.is_none() {
            self.check_for_peak(Amplitude::BaselineStdDev(3.))?;
        }

        let contains_peak = self.params.contains_peak.as_deref().unwrap();

        let mut time_ranges = Vec::with_capacity(self.active_indices.len());

        for (i, &idx) in self.active_indices.iter().enumerate() {
            if !contains_peak[i] {
                time_ranges.push(f64::NAN);
                continue;
            }

            let Some((start_trig_idx, start_trig)) =
                start_trigger.get_first_trigger_element(self.get(i))? else {
                time_ranges.push(f64::NAN);
                continue;
            };

            let start_trigger_time = analysis::interpolate_trigger_time(
                &self.traces[idx],
                start_trig_idx,
                start_trig,
                trigger_time_determination_method,
            )?;

            let Some((stop_trig_idx, stop_trig)) =
                stop_trigger.get_first_trigger_element(self.get(i))? else {
                time_ranges.push(f64::NAN);
                continue;
            };

            let stop_trigger_time = analysis::interpolate_trigger_time(
                &self.traces[idx],
                stop_trig_idx,
                stop_trig,
                trigger_time_determination_method,
            )?;

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
            .enumerate()
            .filter(|&(i, &wave_idx)| filter.filter(&WaveformRef {
                // TODO: It might be more efficient to have two kinds of
                //       filters: WaveformFilters and WaveformArrayFilters
                trace: &self.traces[wave_idx],
                vcal: self.vcal,
                time_channel_disp: self.time_channel_disp,
                global_time_offset: self.global_time_offsets[i],
                peak_polarity: self.peak_polarity,
                params: WaveformParams {
                    contains_peak: self.params.contains_peak.as_ref().map(|arr| arr[i]),
                    baseline: self.params.baselines.as_ref().map(|arr| arr[i]),
                    extremal_element: self.params.extremal_elements.as_ref().map(|arr| arr[i]),
                    peak_height: self.params.peak_heights.as_ref().map(|arr| arr[i]),
                    peak_time: self.params.peak_times.as_ref().map(|arr| arr[i]),
                    peak_area: self.params.peak_areas.as_ref().map(|arr| arr[i]),
                    trigger_time: self.params.trigger_times.as_ref().map(|arr| arr[i]),
                    peak_rise_time: self.params.peak_rise_times.as_ref().map(|arr| arr[i]),
                    peak_decay_time: self.params.peak_decay_times.as_ref().map(|arr| arr[i]),
                },
            }))
            .map(|(i, _)| i)
            .collect::<Vec<_>>();

        WaveformArray {
            traces: self.traces.clone(),
            vcal: self.vcal,
            time_channel_disp: self.time_channel_disp,
            global_time_offsets: self.global_time_offsets.clone(),
            peak_polarity: self.peak_polarity,
            params: WaveformParamArray {
                contains_peak: self.params.contains_peak.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                baselines: self.params.baselines.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                extremal_elements: self.params.extremal_elements.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                peak_heights: self.params.peak_heights.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                peak_areas: self.params.peak_areas.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                peak_times: self.params.peak_times.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                trigger_times: self.params.trigger_times.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                peak_rise_times: self.params.peak_rise_times.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
                peak_decay_times: self.params.peak_decay_times.as_ref().and_then(
                    |arr| Some(active_indices.iter().map(|&i| arr[i]).collect())
                ),
            },
            active_indices,
        }
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
