use nalgebra::{DMatrix, DVector};
use statrs::statistics::Statistics;

use crate::{split_trace_match, trace_match};
use crate::waveforms::waveform::*;


pub fn average_range(
    time_range: (usize, usize),
    integration_range: f64,
    trace: &Trace,
) -> Average {
    let (i1, i2) = time_range;

    trace_match!(trace, arr => {
        let view = arr[i1..i2].iter().map(|x| *x as f64);

        let mean = view.clone().sum::<f64>() / integration_range;
        let stddev = view.std_dev();

        Average { mean, stddev }
    })
}


pub fn get_extremal_element(
    trace: &Trace,
    pulse_polarity: PulsePolarity,
) -> Option<Element> {
    split_trace_match!(
        trace,
        int_arr => {
            match pulse_polarity {
                PulsePolarity::Negative => int_arr
                    .iter()
                    .enumerate()
                    .min_by(|&(_, a), &(_, b)| a.cmp(b)),
                PulsePolarity::Positive => int_arr
                    .iter()
                    .enumerate()
                    .max_by(|&(_, a), &(_, b)| a.cmp(b)),
            }.map(|(extr_idx, extr)| Element {
                idx: extr_idx, value: *extr as f64,
            })
        },
        float_arr => {
            match pulse_polarity {
                PulsePolarity::Negative => float_arr
                    .iter()
                    .enumerate()
                    .min_by(|&(_, a), &(_, b)| a.total_cmp(b)),
                PulsePolarity::Positive => float_arr
                    .iter()
                    .enumerate()
                    .max_by(|&(_, a), &(_, b)| a.total_cmp(b)),
            }.map(|(extr_idx, extr)| Element {
                idx: extr_idx, value: *extr as f64,
            })
        }
    )
}


pub fn check_for_pulse(
    baseline: f64, extremum: f64, threshold: f64,
) -> bool {
    (baseline - extremum).abs() > threshold
}


pub fn quadratic_fit(
    y: DVector<f64>
) -> Result<(f64, f64, f64), AnalysisError> {
    let n_points = y.len();

    let l = -(n_points as i32 / 2) as f64;

    let mut mat = DMatrix::from_element(n_points, 3, 1.);
    for i in 0..n_points {
        let xi = l + i as f64;
        mat[(i, 0)] = xi * xi;
        mat[(i, 1)] = xi;
    }

    let svd = mat.svd(true, true);
    let Ok(beta) = svd.solve(&y, 1e-12) else {
        return Err(AnalysisError::FitError);
    };

    Ok((beta[0], beta[1], beta[2]))
}


pub fn cubic_fit(
    y: DVector<f64>
) -> Result<(f64, f64, f64, f64), AnalysisError> {
    let n_points = y.len();

    let l = -(n_points as i32 / 2) as f64;

    let mut mat = DMatrix::from_element(n_points, 4, 1.);
    for i in 0..n_points {
        let xi = l + i as f64;
        let xi_sq = xi * xi;
        let xi_cb = xi_sq * xi;
        mat[(i, 0)] = xi_cb;
        mat[(i, 1)] = xi_sq;
        mat[(i, 2)] = xi;
    }

    let svd = mat.svd(true, true);
    let Ok(beta) = svd.solve(&y, 1e-12) else {
        return Err(AnalysisError::FitError);
    };

    Ok((beta[0], beta[1], beta[2], beta[3]))
}


pub fn interpolate_pulse_height_and_time(
    trace: &Trace,
    extr_idx: usize,
    extr: f64,
    pulse_height_determination_method: PulseHeightDeterminationMethod,
) -> Result<(f64, f64), AnalysisError> {
    match pulse_height_determination_method {
        PulseHeightDeterminationMethod::ExtremalElement => {
            Ok((extr_idx as f64, extr))
        },
        PulseHeightDeterminationMethod::QuadraticInterpolation => {
            let x1 = extr_idx.max(1);

            if x1 + 1 >= trace.len() {
                return Err(AnalysisError::IndexError {
                    idx: (x1 + 1) as isize, size: trace.len()
                });
            }

            let (y0, y1, y2) = trace_match!(trace, arr => (
                arr[x1 - 1] as f64, arr[x1] as f64, arr[x1 + 1] as f64
            ));

            let dx_extr = 0.5 * (y0 - y2) / ((y0 + y2) - 2. * y1);
            let y_extr = y1 - 0.25 * dx_extr * (y0 - y2);

            Ok((x1 as f64 + dx_extr, y_extr))
        },
        PulseHeightDeterminationMethod::CubicFit { n_points } => {
            let start_idx = extr_idx - n_points / 2;
            let stop_idx = start_idx + n_points;

            if extr_idx < n_points / 2 {
                return Err(AnalysisError::IndexError {
                    idx: extr_idx as isize - (n_points / 2) as isize,
                    size: trace.len()
                });
            }

            if stop_idx >= trace.len() {
                return Err(AnalysisError::IndexError {
                    idx: stop_idx as isize, size: trace.len()
                });
            }

            let y = trace_match!(trace, arr => DVector::from_iterator(
                n_points, (start_idx..stop_idx).map(|i| arr[i] as f64)
            ));

            let (a, b, c, d) = cubic_fit(y)?;

            let dx_extr  = (-b - (b * b - 3. * a * c).sqrt()) / (3. * a);

            if !dx_extr.is_finite() {
                return Err(AnalysisError::FitError);
            }

            let dx_extr_sq = dx_extr * dx_extr;
            let dx_extr_cb = dx_extr_sq * dx_extr;

            let pulse_height = a * dx_extr_cb
                + b * dx_extr_sq
                + c * dx_extr
                + d;

            Ok((extr_idx as f64 + dx_extr, pulse_height))
        },
    }
}


pub fn get_area(
    trace: &Trace,
    start_idx: usize,
    stop_idx: usize,
    integration_range: f64,
    baseline: f64,
) -> Result<f64, AnalysisError> {
    let integral = trace_match!(
        trace, arr => arr[start_idx..stop_idx]
            .iter()
            .map(|&x| x as f64)
            .sum::<f64>()
    );

    Ok((integral - baseline * integration_range).abs())
}


pub fn interpolate_trigger_time(
    trace: &Trace,
    trig_idx: usize,
    trigger_threshold: f64,
    trigger_time_determination_method: TriggerTimeDeterminationMethod,
) -> Result<f64, AnalysisError> {
    match trigger_time_determination_method {
        TriggerTimeDeterminationMethod::FirstThresholdCrossingElement => {
            Ok(trig_idx as f64)
        },
        TriggerTimeDeterminationMethod::LinearInterpolation => {
            trace_match!(trace, arr => {
                let n = arr.len();

                if trig_idx < 1 {
                    return Err(AnalysisError::IndexError { idx: -1, size: n });
                }

                if trig_idx >= n {
                    return Err(AnalysisError::IndexError {
                        idx: trig_idx as isize, size: n
                    });
                }

                let y_l = arr[trig_idx - 1] as f64;
                let y_r = arr[trig_idx] as f64;

                let m = y_r - y_l;

                let dx = (trigger_threshold - y_r) / m;

                Ok(trig_idx as f64 + dx)
            })
        },
        TriggerTimeDeterminationMethod::QuarticInterpolation => {
            trace_match!(trace, arr => {
                let n = trace.len();

                if trig_idx < 2 {
                    return Err(AnalysisError::IndexError {
                        idx: trig_idx as isize - 2, size: n
                    });
                }

                if trig_idx + 2 >= n {
                    return Err(AnalysisError::IndexError {
                        idx: trig_idx as isize + 2, size: n
                    });
                }

                let y0 = arr[trig_idx - 2] as f64;
                let y1 = arr[trig_idx - 1] as f64;
                let y2 = arr[trig_idx] as f64;
                let y3 = arr[trig_idx + 1] as f64;
                let y4 = arr[trig_idx + 2] as f64;

                let h = 1. / 24.;

                let a0 = y2 - trigger_threshold;
                let a1 = 2. * h * (y0 - 8. * y1 + 8. * y3 - y4);
                let a2 = -h * (y0 - 16. * y1 + 30. * y2 - 16. * y3 + y4);
                let a3 = 2. * h * (-y0 + 2. * y1 - 2. * y3 + y4);
                let a4 = h * (y0 - 4. * y2 + 6. * y2 - 4. * y3 + y4);

                let a0_sq = a0 * a0;
                let a1_sq = a1 * a1;
                let a1_cb = a1_sq * a1;

                let dx = (
                    a0 * (2. * a1 * a2 - a0_sq * a3 - a1_cb)
                ) / (
                    a1 * a1_cb
                    - 3. * a0 * a1_sq * a2
                    + 2. * a0_sq * a1 * a3
                    + a0_sq * (a2 * a2 - a0 * a4)
                );

                Ok(trig_idx as f64 + dx)
            })
        },
        TriggerTimeDeterminationMethod::LinearFit { n_points } => {
            let n = trace.len();
            let n_points_halves = n_points / 2;

            let start_idx = trig_idx - n_points_halves;
            let stop_idx = start_idx + n_points;

            if trig_idx < n_points_halves {
                return Err(AnalysisError::IndexError {
                    idx: trig_idx as isize - n_points_halves as isize,
                    size: trace.len()
                });
            }

            if stop_idx >= n {
                return Err(AnalysisError::IndexError {
                    idx: stop_idx as isize, size: n
                });
            }

            let n = n_points as f64;

            let mean_x = (n - 1.) / 2.;

            trace_match!(trace, arr => {
                let m = 12. / (n * (n * n - 1.)) * arr[start_idx..stop_idx]
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| (i as f64 - mean_x) * x as f64)
                    .sum::<f64>();

                let mean_y = arr[start_idx..stop_idx]
                    .iter()
                    .map(|&x| x as f64)
                    .sum::<f64>() / n;

                let t = mean_y - m * mean_x;

                let x0 = (trigger_threshold - t) / m;

                Ok(start_idx as f64 + x0)
            })
        },
        TriggerTimeDeterminationMethod::QuadraticFit { n_points } => {
            let n = trace.len();
            let n_points_halves = n_points / 2;

            let start_idx = trig_idx - n_points_halves;
            let stop_idx = start_idx + n_points;

            if trig_idx < n_points_halves {
                return Err(AnalysisError::IndexError {
                    idx: trig_idx as isize - n_points_halves as isize,
                    size: trace.len()
                });
            }

            if stop_idx >= n {
                return Err(AnalysisError::IndexError {
                    idx: stop_idx as isize, size: n
                });
            }

            let y = trace_match!(trace, arr => DVector::from_iterator(
                n_points, (start_idx..stop_idx).map(|i| arr[i] as f64)
            ));

            let (a, b, mut c) = quadratic_fit(y)?;

            c -= trigger_threshold;

            let dx = (-b + (b * b - 4. * a * c).sqrt()) / (2. * a);

            Ok(trig_idx as f64 + dx)
        },
    }
}
