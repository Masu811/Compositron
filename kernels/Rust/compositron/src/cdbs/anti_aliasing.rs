use nalgebra::{DMatrix, Dyn, MatrixViewMut, Vector2};
use thiserror::Error;

use agg_rust::{
    color::Gray8, path_storage::PathStorage, pixfmt_gray::PixfmtGray8,
    rasterizer_scanline_aa::RasterizerScanlineAa, renderer_base::RendererBase,
    renderer_scanline::render_scanlines_aa_solid,
    rendering_buffer::RenderingBuffer, scanline_u::ScanlineU8,
};

use crate::constants::{SQRT_TWO_HALVES, PI_OVER_FOUR};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AntiAliasingScheme {
    None,
    AxisAligned,
    LargeParallelogram,
    Agg,
}

#[derive(Debug, Error)]
pub enum AntiAliasingError {
    #[error("Input data is too small in size")]
    MatrixTooSmall,

    #[error("Polygon vertices are out of matrix bounds")]
    OutOfBounds,

    #[error("Insufficient number of vertices")]
    InsufficientVertices,

    #[error("Invalid value for polygon vertices (nan or inf)")]
    InvalidVertices,
}

#[derive(Debug, Clone, Copy)]
pub struct Vertex(pub f64, pub f64);

pub trait ConvexToVerticesCounterClockwise {
    fn to_vertices(&self) -> Vec<Vertex>;
}

#[derive(Debug, Clone)]
pub struct Polygon {
    pub vertices: Vec<Vertex>,
}

impl ConvexToVerticesCounterClockwise for Polygon {
    fn to_vertices(&self) -> Vec<Vertex> {
        let x_0 = self.vertices
            .iter()
            .map(|v| v.0)
            .sum::<f64>() / self.vertices.len() as f64;
        let y_0 = self.vertices
            .iter()
            .map(|v| v.1)
            .sum::<f64>() / self.vertices.len() as f64;

        let angles = self.vertices
            .iter()
            .map(|v| (v.1 - y_0).atan2(v.0 - x_0))
            .collect::<Vec<f64>>();

        let mut pairs = self.vertices
            .iter()
            .zip(angles)
            .collect::<Vec<(&Vertex, f64)>>();

        pairs.sort_by(|(_, phi1), (_, phi2)| phi1.total_cmp(phi2));

        pairs.iter().map(|&(&v, _)| v).collect::<Vec<Vertex>>()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Rectangle {
    pub i_min: f64,
    pub i_max: f64,
    pub j_min: f64,
    pub j_max: f64,
}

impl ConvexToVerticesCounterClockwise for Rectangle {
    fn to_vertices(&self) -> Vec<Vertex> {
        vec![
            Vertex(self.j_min, self.i_min),
            Vertex(self.j_min, self.i_max),
            Vertex(self.j_max, self.i_max),
            Vertex(self.j_max, self.i_min),
        ]
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Parallelogram {
    pub center: (f64, f64),
    pub slope_1: f64,
    pub slope_2: f64,
    pub side_length_1: f64,
    pub side_length_2: f64,
}

impl ConvexToVerticesCounterClockwise for Parallelogram {
    fn to_vertices(&self) -> Vec<Vertex> {
        let half_width_1 = self.side_length_1 / 2.;
        let half_width_2 = self.side_length_2 / 2.;

        let w = Vector2::new(1., self.slope_1).normalize() * half_width_1;
        let h = Vector2::new(1., self.slope_2).normalize() * half_width_2;

        let x0 = Vector2::new(self.center.0, self.center.1);

        let v1 = x0 + (w + h);
        let v2 = x0 + (w - h);
        let v3 = x0 + (-w - h);
        let v4 = x0 + (-w + h);

        vec![
            Vertex(v1.x, v1.y), Vertex(v2.x, v2.y),
            Vertex(v3.x, v3.y), Vertex(v4.x, v4.y),
        ]
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Ellipse {
    pub center_i: f64,
    pub center_j: f64,
    pub radius_i: f64,
    pub radius_j: f64,
    pub phi: f64,
}

impl ConvexToVerticesCounterClockwise for Ellipse {
    fn to_vertices(&self) -> Vec<Vertex> {
        let (sin_phi, cos_phi) = self.phi.sin_cos();

        let n_vertices = 128;

        (0..n_vertices)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / n_vertices as f64;
                let (s, c) = t.sin_cos();

                Vertex(
                    self.center_i + self.radius_i * c * cos_phi
                        - self.radius_j * s * sin_phi,
                    self.center_j + self.radius_i * c * sin_phi
                        + self.radius_j * s * cos_phi,
                )
            })
            .collect()
    }
}

fn get_bounding_box(
    vertices: &[Vertex]
) -> Result<(usize, usize, usize, usize), AntiAliasingError> {
    let (Some(left_col), Some(right_col), Some(lower_row), Some(upper_row)) = (
        vertices.iter().map(|v| v.0).min_by(|x1, x2| x1.total_cmp(x2)),
        vertices.iter().map(|v| v.0).max_by(|x1, x2| x1.total_cmp(x2)),
        vertices.iter().map(|v| v.1).max_by(|y1, y2| y1.total_cmp(y2)),
        vertices.iter().map(|v| v.1).min_by(|y1, y2| y1.total_cmp(y2)),
    ) else {
        return Err(AntiAliasingError::InvalidVertices);
    };

    let left_col = left_col.round() as usize;
    let right_col = right_col.round() as usize;
    let lower_row = lower_row.round() as usize;
    let upper_row = upper_row.round() as usize;

    Ok((left_col, right_col, lower_row, upper_row))
}

/// Flip the value of all pixels whose center is below the given edge
fn invert_from_edge(
    view: &mut MatrixViewMut<'_, f64, Dyn, Dyn>, v1: Vertex, v2: Vertex
) {
    let m = (v1.1 - v2.1) / (v1.0 - v2.0);
    let t = v1.1 - m * v1.0;

    let left_col = v1.0.min(v2.0).ceil() as usize;
    let right_col = v1.0.max(v2.0).floor() as usize;

    let nrows = view.nrows();

    for (j, mut col) in view
        .columns_range_mut(left_col..right_col + 1)
        .column_iter_mut()
        .enumerate()
    {
        let row = (m * (j + left_col) as f64 + t).ceil() as usize;
        let view_row = row.min(nrows);
        col.rows_range_mut(view_row..)
            .iter_mut()
            .for_each(|z: &mut f64| *z = 1. - *z);
    }
}

/// Set weight to 1 if pixel center is inside the polygon, otherwise 0.
pub fn no_aa(
    nrows: usize,
    ncols: usize,
    polygon: &impl ConvexToVerticesCounterClockwise,
) -> Result<DMatrix<f64>, AntiAliasingError> {
    if nrows < 3 || ncols < 3 {
        return Err(AntiAliasingError::MatrixTooSmall);
    }

    let vertices = polygon.to_vertices();

    if vertices.len() < 3 {
        return Err(AntiAliasingError::InsufficientVertices);
    }

    let (left_col, right_col, lower_row, upper_row) = get_bounding_box(&vertices)?;

    if right_col > ncols - 1 || lower_row > nrows - 1 {
        return Err(AntiAliasingError::OutOfBounds);
    }

    let nrows_view = lower_row - upper_row + 1;
    let ncols_view = right_col - left_col + 1;

    let mut weights = DMatrix::zeros(nrows, ncols);

    let mut view = weights.view_mut(
        (upper_row, left_col), (nrows_view, ncols_view)
    );

    let pol_in_view = Polygon {
        vertices: vertices
            .iter()
            .map(|v| Vertex(v.0 - left_col as f64, v.1 - upper_row as f64))
            .collect(),
    };

    let l = vertices.len();

    for i in 0..l {
        let v1 = pol_in_view.vertices[i];
        let v2 = pol_in_view.vertices[(i + 1) % l];

        invert_from_edge(&mut view, v1, v2);
    }

    invert_from_edge(
        &mut view, pol_in_view.vertices[l - 1], pol_in_view.vertices[0],
    );

    Ok(weights)
}

fn get_corners(
    view: &MatrixViewMut<'_, f64, Dyn, Dyn>,
    nrows: usize,
    ncols: usize,
) -> (f64, f64, f64, f64) {
    (
        view[(0, 0)],
        view[(nrows - 1, 0)],
        view[(nrows - 1, ncols - 1)],
        view[(0, ncols - 1)],
    )
}

pub fn axis_aligned_aa(
    nrows: usize,
    ncols: usize,
    rectangle: &Rectangle,
) -> Result<DMatrix<f64>, AntiAliasingError> {
    if nrows < 3 || ncols < 3 {
        return Err(AntiAliasingError::MatrixTooSmall);
    }

    let left_col = rectangle.j_min.round() as usize;
    let right_col = rectangle.j_max.round() as usize;
    let upper_row = rectangle.i_min.round() as usize;
    let lower_row = rectangle.j_max.round() as usize;

    if right_col > ncols - 1 || lower_row > nrows - 1 {
        return Err(AntiAliasingError::OutOfBounds);
    }

    let nrows_view = lower_row - upper_row + 1;
    let ncols_view = right_col - left_col + 1;

    let mut weights = DMatrix::zeros(nrows, ncols);

    let mut view = weights.view_mut(
        (upper_row, left_col), (nrows_view, ncols_view)
    );

    view.fill(1.);

    if left_col == right_col {
        let (a, b, _, _) = get_corners(&view, nrows_view, ncols_view);

        let overlap = rectangle.j_max - rectangle.j_min;
        view.fill_column(left_col, overlap);

        view[(0, 0)] = a * overlap;
        view[(nrows_view - 1, 0)] = b * overlap;
    } else {
        let (a, b, c, d) = get_corners(&view, nrows_view, ncols_view);

        let left_overlap = 1. - (rectangle.j_min - 0.5).fract();
        view.fill_column(0, left_overlap);

        let right_overlap = (rectangle.j_max - 0.5) % 1.;
        view.fill_column(ncols_view - 1, right_overlap);

        view[(0, 0)] = a * left_overlap;
        view[(nrows_view - 1, 0)] = b * left_overlap;
        view[(nrows_view - 1, ncols_view - 1)] = c * right_overlap;
        view[(0, ncols_view - 1)] = d * right_overlap;
    }

    if lower_row == upper_row {
        let (a, _, _, d) = get_corners(&view, nrows_view, ncols_view);

        let overlap = rectangle.i_max - rectangle.i_min;
        view.fill_row(0, overlap);

        view[(0, 0)] = a * overlap;
        view[(0, ncols_view - 1)] = d * overlap;
    } else {
        let (a, b, c, d) = get_corners(&view, nrows_view, ncols_view);

        let upper_overlap = 1. - (rectangle.i_min - 0.5) % 1.;
        view.fill_row(0, upper_overlap);

        let lower_overlap = (rectangle.i_max - 0.5) % 1.;
        view.fill_row(nrows_view - 1, lower_overlap);

        view[(0, 0)] = a * upper_overlap;
        view[(nrows_view - 1, 0)] = b * lower_overlap;
        view[(nrows_view - 1, ncols_view - 1)] = c * lower_overlap;
        view[(0, ncols_view - 1)] = d * upper_overlap;
    }

    Ok(weights)
}

pub fn parallelogram_aa(
    nrows: usize,
    ncols: usize,
    parallelogram: &Parallelogram,
) -> Result<DMatrix<f64>, AntiAliasingError> {
    if nrows < 3 || ncols < 3 {
        return Err(AntiAliasingError::MatrixTooSmall);
    }

    let vertices = parallelogram.to_vertices();

    let (left_col, right_col, lower_row, upper_row) = get_bounding_box(&vertices)?;

    if right_col > ncols - 1 || lower_row > nrows - 1 {
        return Err(AntiAliasingError::OutOfBounds);
    }

    let nrows_view = lower_row - upper_row + 1;
    let ncols_view = right_col - left_col + 1;

    let mut weights = DMatrix::zeros(nrows, ncols);

    let mut view = weights.view_mut(
        (upper_row, left_col), (nrows_view, ncols_view)
    );

    view.fill(1.);

    let par_in_view = Parallelogram {
        center: (
            parallelogram.center.0 - left_col as f64,
            parallelogram.center.1 - upper_row as f64
        ),
        slope_1: parallelogram.slope_1,
        slope_2: parallelogram.slope_2,
        side_length_1: parallelogram.side_length_1,
        side_length_2: parallelogram.side_length_2,
    };

    let w = Vector2::new(1., parallelogram.slope_1).normalize();
    let h = Vector2::new(1., parallelogram.slope_2).normalize();

    let theta = w.angle(&h);

    let perp_dist_1 = 0.5 * par_in_view.side_length_2 * theta.sin();
    let perp_dist_2 = 0.5 * par_in_view.side_length_1 * theta.sin();

    for (dx, dy, w, slope) in [
        (w.x, w.y, perp_dist_1, par_in_view.slope_1),
        (h.x, h.y, perp_dist_2, par_in_view.slope_2),
    ] {
        let phi = slope.abs().atan();

        let c = -dy * par_in_view.center.0 + dx * par_in_view.center.1;

        let min_edge_dist = SQRT_TWO_HALVES * (PI_OVER_FOUR - phi).abs().sin();
        let max_edge_dist = SQRT_TWO_HALVES * (PI_OVER_FOUR + phi).sin();

        let crit_dist_1 = w - max_edge_dist;
        let crit_dist_2 = w - min_edge_dist;
        let crit_dist_3 = w + min_edge_dist;
        let crit_dist_4 = w + max_edge_dist;

        let sin_phi = dy.abs();
        let cos_phi = dx.abs();
        let tan_phi = slope.abs();

        let one_over_tan_phi = 1. / tan_phi;
        let one_over_sin_phi = 1. / sin_phi;
        let one_over_cos_phi = 1. / cos_phi;
        let one_over_2sincos = 1. / (2. * sin_phi * cos_phi);

        for col in 0..ncols_view {
            let x = col as f64;

            let mut acc = dy * x + c;
            let dacc = -dx;

            for row in 0..nrows_view {
                let dist = acc.abs();
                acc += dacc;

                let coverage = if dist < crit_dist_1 {
                    // fully inside
                    continue;
                } else if dist < crit_dist_2 {
                    // One corner of the pixel is outside
                    let d = max_edge_dist - (w - dist);

                    1. - d * d * one_over_2sincos
                } else if dist < crit_dist_3 {
                    // Two corners insde, two outside
                    let d = max_edge_dist - (dist - w);

                    let a = d * d * one_over_2sincos;

                    let b = if tan_phi > 1. {
                        0.5 * one_over_tan_phi * (d * one_over_cos_phi - 1.).powi(2)
                    } else {
                        0.5 * tan_phi * (d * one_over_sin_phi - 1.).powi(2)
                    };

                    a - b
                } else if dist < crit_dist_4 {
                    // One corner of the pixel is inside
                    let d = max_edge_dist - (dist - w);

                    d * d * one_over_2sincos
                } else {
                    // fully outside
                    0.
                };

                view[(row, col)] *= coverage;
            }
        }
    }

    Ok(weights)
}

pub fn agg_aa(
    nrows: usize,
    ncols: usize,
    polygon: &impl ConvexToVerticesCounterClockwise,
) -> Result<DMatrix<f64>, AntiAliasingError> {
    if nrows < 3 || ncols < 3 {
        return Err(AntiAliasingError::MatrixTooSmall);
    }

    let vertices = polygon.to_vertices();

    if vertices.len() < 3 {
        return Err(AntiAliasingError::InsufficientVertices);
    }

    let (left_col, right_col, lower_row, upper_row) = get_bounding_box(&vertices)?;

    if right_col > ncols - 1 || lower_row > nrows - 1 {
        return Err(AntiAliasingError::OutOfBounds);
    }

    let nrows_view = lower_row - upper_row + 1;
    let ncols_view = right_col - left_col + 1;

    let mut weights = DMatrix::zeros(nrows, ncols);

    let mut view = weights.view_mut(
        (upper_row, left_col), (nrows_view, ncols_view)
    );

    let vertices = vertices
        .iter()
        .map(|v| Vertex(v.0 - left_col as f64, v.1 - upper_row as f64))
        .collect::<Vec<Vertex>>();

    let mut path = PathStorage::new();
    path.move_to(vertices[0].0, vertices[0].1);
    for vertex in vertices.iter().skip(1) {
        path.line_to(vertex.0, vertex.1);
    }
    path.close_polygon(0);

    let mut buf = vec![0u8; ncols_view * nrows_view];

    unsafe {
        let mut rbuf = RenderingBuffer::new_with_buf(
            buf.as_mut_ptr(),
            ncols_view as u32,
            nrows_view as u32,
            ncols_view as i32,
        );

        let pixfmt = PixfmtGray8::new(&mut rbuf);
        let mut ren_base = RendererBase::new(pixfmt);
        ren_base.clear(&Gray8::new(0, 255));

        let mut ras = RasterizerScanlineAa::new();
        let mut sl = ScanlineU8::new();
        ras.add_path(&mut path, 0);
        render_scanlines_aa_solid(
            &mut ras, &mut sl, &mut ren_base, &Gray8::new(255, 255)
        );
    }

    for j in 0..ncols_view {
        for i in 0..nrows_view {
            view[(i, j)] = buf[i * ncols_view + j] as f64 / 255.;
        }
    }

    Ok(weights)
}

fn compute_intersection(
    a: &Vertex, b: &Vertex, v1: &Vertex, v2: &Vertex
) -> Vertex {
    let m1 = (a.1 - b.1) / (a.0 - b.0);
    let t1 = a.1 - m1 * a.0;

    let m2 = (v1.1 - v2.1) / (v1.0 - v2.0);
    let t2 = v1.1 - m2 * v1.0;

    let x0 = (t2 - t1) / (m1 - m2);
    let y0 = m1 * x0 + t1;

    Vertex(x0, y0)
}

fn inside(v: &Vertex, v1: &Vertex, v2: &Vertex) -> bool {
    let a = Vector2::from_vec(vec![v1.0 - v2.0, v1.1 - v2.1]);
    let b = Vector2::from_vec(vec![v1.0 - v.0, v1.1 - v.1]);

    a.x * b.y - a.y * b.x < 0.
}

/// Intersection by Sutherland-Hodgman clipping
pub fn intersect_convex_polygons(
    a: &impl ConvexToVerticesCounterClockwise,
    b: &impl ConvexToVerticesCounterClockwise,
) -> Result<Polygon, AntiAliasingError> {
    let va = a.to_vertices();
    let vb = b.to_vertices();

    if va.len() < 3 || vb.len() < 3 {
        return Err(AntiAliasingError::InsufficientVertices);
    }

    let la = va.len();

    let mut output = vb;

    for i in 0..la {
        let edge_start = va[i];
        let edge_end = va[(i + 1) % la];

        let input = output;
        output = Vec::new();

        for i in 0..input.len() {
            let curr = input[i];
            let prev = input[(i + input.len() - 1) % input.len()];

            let curr_inside = inside(&curr, &edge_start, &edge_end);
            let prev_inside = inside(&prev, &edge_start, &edge_end);

            if curr_inside {
                if !prev_inside {
                    output.push(compute_intersection(
                        &prev, &curr, &edge_start, &edge_end
                    ));
                }
                output.push(curr);
            } else if prev_inside {
                output.push(compute_intersection(
                    &prev, &curr, &edge_start, &edge_end
                ));
            }
        }
    }

    Ok(Polygon { vertices: output })
}
