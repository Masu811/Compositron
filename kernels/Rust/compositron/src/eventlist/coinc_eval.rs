use std::ops::Sub;


// "stable" if the output is independent of the order of detectors passed in
// "single match" if one event contributes to at most one coinc event
#[derive(Debug, Clone, Copy)]
pub enum CoincEvalMode {
    EachMatch,  // stable / multi match
    FirstMatch,  // stable / single match
    ClosestMatch,  // stable / single match
    UniqueMatch,  // stable / single match
    FirstPivot,  // unstable / multi match
}


// TODO:
// Everything for arbitrary many detectors
// Everything for multi-threading and streaming data from file


pub fn to_hist_each_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_c: &[T],
    tarr_r: &[T],
    barr_c: &[u16],
    barr_r: &[u16],
    tidcs_c: &[usize],
    tidcs_r: &[usize],
    hist: &mut [u64],
    nrows: usize,
    coinc_time_window: T,
) {
    let mut offset_r = 0;

    for (&tidx_c, &bin_c) in tidcs_c.iter().zip(barr_c) {
        if bin_c == u16::MAX {
            continue;
        }

        let t_c = tarr_c[tidx_c];

        for (&tidx_r, &bin_r) in tidcs_c[offset_r..].iter().zip(&barr_r[offset_r..]) {
            let t_r = tarr_r[tidx_r];

            if bin_r < u16::MAX && (t_r > t_c || t_c - t_r <= coinc_time_window) {
                break;
            }

            offset_r += 1;
        }

        let col = &mut hist[
            (bin_c as usize * nrows)..((bin_c as usize + 1) * nrows)
        ];

        for (&tidx_r, &bin_r) in tidcs_r[offset_r..].iter().zip(&barr_r[offset_r..]) {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c && t_r - t_c > coinc_time_window {
                break;
            }

            if bin_r < u16::MAX {
                col[bin_r as usize] += 1;
            }
        }
    }
}


pub fn to_celist_each_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_1: &[T],
    tarr_2: &[T],
    tidcs_1: &[usize],
    tidcs_2: &[usize],
    coinc_time_window: T,
) -> (Vec<usize>, Vec<usize>) {
    let mut idcs_1 = Vec::new();
    let mut idcs_2 = Vec::new();

    let mut offset_2 = 0;

    for &tidx_1 in tidcs_1 {
        let t_1 = tarr_1[tidx_1];

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 || t_1 - t_2 <= coinc_time_window {
                break;
            }

            offset_2 += 1;
        }

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 && t_2 - t_1 > coinc_time_window {
                break;
            }

            idcs_1.push(tidx_1);
            idcs_2.push(tidx_2);
        }
    }

    (idcs_1, idcs_2)
}


pub fn to_hist_first_pivot<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_c: &[T],
    tarr_r: &[T],
    barr_c: &[u16],
    barr_r: &[u16],
    tidcs_c: &[usize],
    tidcs_r: &[usize],
    hist: &mut [u64],
    nrows: usize,
    coinc_time_window: T,
) {
    let mut offset_r = 0;

    for (&tidx_c, &bin_c) in tidcs_c.iter().zip(barr_c) {
        if bin_c == u16::MAX {
            continue;
        }

        let t_c = tarr_c[tidx_c];

        for &tidx_r in &tidcs_r[offset_r..] {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c {
                break;
            }

            offset_r += 1;
        }

        if {
            offset_r < tidcs_r.len()
            && barr_r[offset_r] < u16::MAX
            && tarr_r[tidcs_r[offset_r]] - t_c <= coinc_time_window
        } {
            // pass
        } else if {
            offset_r > 0
            && barr_r[offset_r - 1] < u16::MAX
            && t_c - tarr_r[tidcs_r[offset_r - 1]] <= coinc_time_window
        } {
            offset_r -= 1;
        } else {
            continue;
        }

        let bin_r = barr_r[offset_r];

        if bin_r < u16::MAX {
            hist[bin_c as usize * nrows + bin_r as usize] += 1;
        }
    }
}


pub fn to_celist_first_pivot<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_1: &[T],
    tarr_2: &[T],
    tidcs_1: &[usize],
    tidcs_2: &[usize],
    coinc_time_window: T,
) -> (Vec<usize>, Vec<usize>) {
    let mut idcs_1 = Vec::new();
    let mut idcs_2 = Vec::new();

    let mut offset_2 = 0;

    for &tidx_1 in tidcs_1 {
        let t_1 = tarr_1[tidx_1];

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 {
                break;
            }

            offset_2 += 1;
        }

        if offset_2 < tidcs_2.len() && tarr_2[tidcs_2[offset_2]] - t_1 <= coinc_time_window {
            // pass
        } else if offset_2 > 0 && t_1 - tarr_2[tidcs_2[offset_2 - 1]] <= coinc_time_window {
            offset_2 -= 1;
        } else {
            continue;
        }

        idcs_1.push(tidx_1);
        idcs_2.push(tidcs_2[offset_2]);
    }

    (idcs_1, idcs_2)
}


pub fn to_hist_closest_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_c: &[T],
    tarr_r: &[T],
    barr_c: &[u16],
    barr_r: &[u16],
    tidcs_c: &[usize],
    tidcs_r: &[usize],
    hist: &mut [u64],
    nrows: usize,
    coinc_time_window: T,
) {
    let mut offset_r = 0;

    for (i_c, (&tidx_c, &bin_c)) in tidcs_c.iter().zip(barr_c).enumerate() {
        if bin_c == u16::MAX {
            continue;
        }

        let t_c = tarr_c[tidx_c];

        for &tidx_r in &tidcs_r[offset_r..] {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c {
                break;
            }

            offset_r += 1;
        }

        if offset_r == 0 {
            let t_r = tarr_r[tidcs_r[offset_r]];

            let my_dist_r = t_r - t_c;

            if my_dist_r <= coinc_time_window {
                if i_c + 1 >= tidcs_c.len() {
                    // pass
                } else if tarr_c[tidcs_c[i_c + 1]] < t_r {
                    continue;
                } else if tarr_c[tidcs_c[i_c + 1]] - t_r < my_dist_r {
                    continue;
                }
            } else {
                continue;
            }
        } else if offset_r == tidcs_r.len() {
            let t_r = tarr_r[tidcs_r[offset_r - 1]];

            let my_dist_l = t_c - t_r;

            if my_dist_l <= coinc_time_window {
                if i_c == 0 {
                    // pass
                } else if t_r < tarr_c[tidcs_c[i_c - 1]] {
                    continue;
                } else if t_r - tarr_c[tidcs_c[i_c - 1]] <= my_dist_l {
                    continue;
                }

                offset_r -= 1;
            } else {
                continue;
            }
        } else {
            let my_dist_r = tarr_r[tidcs_r[offset_r]] - t_c;
            let my_dist_l = t_c - tarr_r[tidcs_r[offset_r - 1]];

            if my_dist_l <= coinc_time_window {
                if my_dist_r < my_dist_l {
                    if i_c + 1 >= tarr_c.len() {
                        // pass
                    } else if tarr_c[i_c + 1] < tarr_r[tidcs_r[offset_r]] {
                        continue;
                    } else if tarr_c[i_c + 1] - tarr_r[tidcs_r[offset_r]] < my_dist_r {
                        continue;
                    }
                } else {
                    if i_c == 0 {
                        // pass
                    } else if tarr_r[tidcs_r[offset_r - 1]] < tarr_c[tidcs_c[i_c - 1]] {
                        continue;
                    } else if tarr_r[tidcs_r[offset_r - 1]] - tarr_c[tidcs_c[i_c - 1]] <= my_dist_l {
                        continue;
                    }

                    offset_r -= 1;
                }
            } else if my_dist_r <= coinc_time_window {
                if i_c + 1 >= tidcs_c.len() {
                    // pass
                } else if tarr_c[tidcs_c[i_c + 1]] < tarr_r[tidcs_r[offset_r]] {
                    continue;
                } else if tarr_c[i_c + 1] - tarr_r[tidcs_r[offset_r]] < my_dist_r {
                    continue;
                }
            } else {
                continue;
            }
        }

        let bin_r = barr_r[offset_r];

        if bin_r < u16::MAX {
            hist[bin_c as usize * nrows + bin_r as usize] += 1;
        }
        offset_r += 1;
    }
}


pub fn to_celist_closest_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_1: &[T],
    tarr_2: &[T],
    tidcs_1: &[usize],
    tidcs_2: &[usize],
    coinc_time_window: T,
) -> (Vec<usize>, Vec<usize>) {
    let mut idcs_1 = Vec::new();
    let mut idcs_2 = Vec::new();

    let mut offset_2 = 0;

    for (i_1, &tidx_1) in tidcs_1.iter().enumerate() {
        let t_1 = tarr_1[tidx_1];

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 {
                break;
            }

            offset_2 += 1;
        }

        if offset_2 == 0 {
            let t_2 = tarr_2[tidcs_2[offset_2]];

            let my_dist_r = t_2 - t_1;

            if my_dist_r <= coinc_time_window {
                if i_1 + 1 >= tidcs_1.len() {
                    // pass
                } else if tarr_1[tidcs_1[i_1 + 1]] < t_2 {
                    continue;
                } else if tarr_1[tidcs_1[i_1 + 1]] - t_2 < my_dist_r {
                    continue;
                }
            } else {
                continue;
            }
        } else if offset_2 == tidcs_2.len() {
            let t_2 = tarr_2[tidcs_2[offset_2 - 1]];

            let my_dist_l = t_1 - t_2;

            if my_dist_l <= coinc_time_window {
                if i_1 == 0 {
                    // pass
                } else if t_2 < tarr_1[tidcs_1[i_1 - 1]] {
                    continue;
                } else if t_2 - tarr_1[tidcs_1[i_1 - 1]] <= my_dist_l {
                    continue;
                }

                offset_2 -= 1;
            } else {
                continue;
            }
        } else {
            let my_dist_r = tarr_2[tidcs_2[offset_2]] - t_1;
            let my_dist_l = t_1 - tarr_2[tidcs_2[offset_2 - 1]];

            if my_dist_l <= coinc_time_window {
                if my_dist_r < my_dist_l {
                    if i_1 + 1 >= tarr_1.len() {
                        // pass
                    } else if tarr_1[tidcs_1[i_1 + 1]] < tarr_2[tidcs_2[offset_2]] {
                        continue;
                    } else if tarr_1[tidcs_1[i_1 + 1]] - tarr_2[tidcs_2[offset_2]] < my_dist_r {
                        continue;
                    }
                } else {
                    if i_1 == 0 {
                        // pass
                    } else if tarr_2[tidcs_2[offset_2 - 1]] < tarr_1[tidcs_1[i_1 - 1]] {
                        continue;
                    } else if tarr_2[tidcs_2[offset_2 - 1]] - tarr_1[tidcs_1[i_1 - 1]] <= my_dist_l {
                        continue;
                    }

                    offset_2 -= 1;
                }
            } else if my_dist_r <= coinc_time_window {
                if i_1 + 1 >= tarr_1.len() {
                    // pass
                } else if tarr_1[tidcs_1[i_1 + 1]] < tarr_2[tidcs_2[offset_2]] {
                    continue;
                } else if tarr_1[tidcs_1[i_1 + 1]] - tarr_2[tidcs_2[offset_2]] < my_dist_r {
                    continue;
                }
            } else {
                continue;
            }
        }

        idcs_1.push(tidx_1);
        idcs_2.push(tidcs_2[offset_2]);
        offset_2 += 1;
    }

    (idcs_1, idcs_2)
}


pub fn to_hist_first_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_c: &[T],
    tarr_r: &[T],
    barr_c: &[u16],
    barr_r: &[u16],
    tidcs_c: &[usize],
    tidcs_r: &[usize],
    hist: &mut [u64],
    nrows: usize,
    coinc_time_window: T,
) {
    let mut offset_r = 0;

    for (&tidx_c, &bin_c) in tidcs_c.iter().zip(barr_c) {
        let t_c = tarr_c[tidx_c];

        for &tidx_r in tidcs_r[offset_r..].iter() {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c {
                if t_r - t_c <= coinc_time_window {
                    let bin_r = barr_r[offset_r];
                    if bin_c < u16::MAX && bin_r < u16::MAX {
                        hist[bin_c as usize * nrows + bin_r as usize] += 1;
                    }
                    offset_r += 1;
                }
                break;
            } else if t_c - t_r <= coinc_time_window {
                let bin_r = barr_r[offset_r];
                if bin_c < u16::MAX && bin_r < u16::MAX {
                    hist[bin_c as usize * nrows + bin_r as usize] += 1;
                }
                offset_r += 1;
                break;
            }

            offset_r += 1;
        }
    }
}


pub fn to_celist_first_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_1: &[T],
    tarr_2: &[T],
    tidcs_1: &[usize],
    tidcs_2: &[usize],
    coinc_time_window: T,
) -> (Vec<usize>, Vec<usize>) {
    let mut idcs_1 = Vec::new();
    let mut idcs_2 = Vec::new();

    let mut offset_2 = 0;

    for &tidx_1 in tidcs_1 {
        let t_1 = tarr_1[tidx_1];

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 {
                if t_2 - t_1 <= coinc_time_window {
                    idcs_1.push(tidx_1);
                    idcs_2.push(tidx_2);
                    offset_2 += 1;
                }
                break;
            } else if t_1 - t_2 <= coinc_time_window {
                idcs_1.push(tidx_1);
                idcs_2.push(tidx_2);
                offset_2 += 1;
                break;
            }

            offset_2 += 1;
        }
    }

    (idcs_1, idcs_2)
}


pub fn to_hist_unique_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_c: &[T],
    tarr_r: &[T],
    barr_c: &[u16],
    barr_r: &[u16],
    tidcs_c: &[usize],
    tidcs_r: &[usize],
    hist: &mut [u64],
    nrows: usize,
    coinc_time_window: T,
) {
    let mut offset_r = 0;

    for (i_c, (&tidx_c, &bin_c)) in tidcs_c.iter().zip(barr_c).enumerate() {
        if bin_c == u16::MAX {
            continue;
        }

        let t_c = tarr_c[tidx_c];

        for &tidx_r in &tidcs_r[offset_r..] {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c || t_c - t_r <= coinc_time_window {
                break;
            }

            offset_r += 1;
        }

        let mut found_exactly_one_neighbor = false;

        for &tidx_r in tidcs_r[offset_r..].iter() {
            let t_r = tarr_r[tidx_r];

            if t_r > t_c && t_r - t_c > coinc_time_window {
                break;
            }

            if found_exactly_one_neighbor {
                found_exactly_one_neighbor = false;
                break;
            } else {
                found_exactly_one_neighbor = true;
            }
        }

        if !found_exactly_one_neighbor {
            continue;
        }

        let t_r = tarr_r[tidcs_r[offset_r]];

        if i_c > 0 {
            let prev_t_c = tarr_c[tidcs_c[i_c - 1]];
            if prev_t_c < t_r {
                if t_r - prev_t_c <= coinc_time_window {
                    continue;
                }
            } else if prev_t_c - t_r <= coinc_time_window {
                continue;
            }
        }

        if i_c + 1 < tidcs_c.len() {
            let next_t_c = tarr_c[tidcs_c[i_c + 1]];
            if next_t_c < t_r {
                if t_r - next_t_c <= coinc_time_window {
                    continue;
                }
            } else if next_t_c - t_r <= coinc_time_window {
                continue;
            }
        }

        let bin_r = barr_r[offset_r];

        if bin_r < u16::MAX {
            hist[bin_c as usize * nrows + bin_r as usize] += 1;
        }
    }
}


pub fn to_celist_unique_match<T: Copy + Sub<Output = T> + PartialOrd>(
    tarr_1: &[T],
    tarr_2: &[T],
    tidcs_1: &[usize],
    tidcs_2: &[usize],
    coinc_time_window: T,
) -> (Vec<usize>, Vec<usize>) {
    let mut idcs_1 = Vec::new();
    let mut idcs_2 = Vec::new();

    let mut offset_2 = 0;

    for (i_1, &tidx_1) in tidcs_1.iter().enumerate() {
        let t_1 = tarr_1[tidx_1];

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 || t_1 - t_2 <= coinc_time_window {
                break;
            }

            offset_2 += 1;
        }

        let mut found_exactly_one_neighbor = false;

        for &tidx_2 in &tidcs_2[offset_2..] {
            let t_2 = tarr_2[tidx_2];

            if t_2 > t_1 && t_2 - t_1 > coinc_time_window {
                break;
            }

            if found_exactly_one_neighbor {
                found_exactly_one_neighbor = false;
                break;
            } else {
                found_exactly_one_neighbor = true;
            }
        }

        if !found_exactly_one_neighbor {
            continue;
        }

        let t_2 = tarr_2[tidcs_2[offset_2]];

        if i_1 > 0 {
            let prev_t_1 = tarr_1[tidcs_1[i_1 - 1]];
            if prev_t_1 < t_2 {
                if t_2 - prev_t_1 <= coinc_time_window {
                    continue;
                }
            } else if prev_t_1 - t_2 <= coinc_time_window {
                continue;
            }
        }

        if i_1 + 1 < tidcs_1.len() {
            let next_t_1 = tarr_1[tidcs_1[i_1 + 1]];
            if next_t_1 < t_2 {
                if t_2 - next_t_1 <= coinc_time_window {
                    continue;
                }
            } else if next_t_1 - t_2 <= coinc_time_window {
                continue;
            }
        }

        idcs_1.push(tidx_1);
        idcs_2.push(tidcs_2[offset_2]);
    }

    (idcs_1, idcs_2)
}
