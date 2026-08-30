use nalgebra::{DMatrix, Vector2};
use thiserror::Error;

use agg_rust::{
    color::Gray8, path_storage::PathStorage, pixfmt_gray::PixfmtGray8,
    rasterizer_scanline_aa::RasterizerScanlineAa, renderer_base::RendererBase,
    renderer_scanline::render_scanlines_aa_solid,
    rendering_buffer::RenderingBuffer, scanline_u::ScanlineU8,
};

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
pub struct Vertex {
    pub x: f64,
    pub y: f64,
}

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
            .map(|v| v.x)
            .sum::<f64>() / self.vertices.len() as f64;
        let y_0 = self.vertices
            .iter()
            .map(|v| v.y)
            .sum::<f64>() / self.vertices.len() as f64;

        let angles = self.vertices
            .iter()
            .map(|v| (v.y - y_0).atan2(v.x - x_0))
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
            Vertex { x: self.j_min, y: self.i_min },
            Vertex { x: self.j_min, y: self.i_max },
            Vertex { x: self.j_max, y: self.i_max },
            Vertex { x: self.j_max, y: self.i_min },
        ]
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Parallelogram {
    pub center_i: f64,
    pub center_j: f64,
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

        let x0 = Vector2::new(self.center_j, self.center_i);

        let v1 = x0 + (w - h);
        let v2 = x0 + (w + h);
        let v3 = x0 + (-w + h);
        let v4 = x0 + (-w - h);

        vec![
            Vertex { x: v1.x, y: v1.y }, Vertex { x: v2.x, y: v2.y },
            Vertex { x: v3.x, y: v3.y }, Vertex { x: v4.x, y: v4.y },
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

                Vertex {
                    x: self.center_j + self.radius_i * c * sin_phi
                        + self.radius_j * s * cos_phi,
                    y: self.center_i + self.radius_i * c * cos_phi
                        - self.radius_j * s * sin_phi,
                }
            })
            .collect()
    }
}

pub fn get_bounding_box(
    vertices: &[Vertex]
) -> Result<(usize, usize, usize, usize), AntiAliasingError> {
    let (Some(left_col), Some(right_col), Some(lower_row), Some(upper_row)) = (
        vertices.iter().map(|v| v.x).min_by(|x1, x2| x1.total_cmp(x2)),
        vertices.iter().map(|v| v.x).max_by(|x1, x2| x1.total_cmp(x2)),
        vertices.iter().map(|v| v.y).max_by(|y1, y2| y1.total_cmp(y2)),
        vertices.iter().map(|v| v.y).min_by(|y1, y2| y1.total_cmp(y2)),
    ) else {
        return Err(AntiAliasingError::InvalidVertices);
    };

    let left_col = left_col.round() as usize;
    let right_col = right_col.round() as usize;
    let lower_row = lower_row.round() as usize;
    let upper_row = upper_row.round() as usize;

    Ok((left_col, right_col, lower_row, upper_row))
}

pub fn agg_aa(
    nrows: usize,
    ncols: usize,
    polygon: &impl ConvexToVerticesCounterClockwise,
) -> Result<DMatrix<u8>, AntiAliasingError> {
    if nrows < 3 || ncols < 3 {
        return Err(AntiAliasingError::MatrixTooSmall);
    }

    let vertices = polygon.to_vertices();

    if vertices.len() < 3 {
        return Err(AntiAliasingError::InsufficientVertices);
    }

    let mut weights = DMatrix::zeros(nrows, ncols);

    // Transform our coords (pixel centers at integer coords)
    // to agg's coords (pixel boundaries at integer coords)
    let vertices = vertices
        .iter()
        .map(|v| Vertex { x: v.x + 0.5, y: v.y + 0.5 })
        .collect::<Vec<Vertex>>();

    if vertices.iter().any(|&v| {
        v.x < 0. || v.x > ncols as f64 || v.y < 0. || v.y > nrows as f64
    }) {
        return Err(AntiAliasingError::OutOfBounds);
    }

    // Coords are transposed because agg is row-major
    let mut path = PathStorage::new();
    path.move_to(vertices[0].y, vertices[0].x);
    for vertex in vertices.iter().skip(1) {
        path.line_to(vertex.y, vertex.x);
    }
    path.close_polygon(0);

    unsafe {
        let mut rbuf = RenderingBuffer::new_with_buf(
            weights.as_mut_ptr(),
            nrows as u32,
            ncols as u32,
            nrows as i32,
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

    Ok(weights)
}

fn compute_intersection(
    a: &Vertex, b: &Vertex, v1: &Vertex, v2: &Vertex
) -> Vertex {
    let m1 = (a.y - b.y) / (a.x - b.x);
    let t1 = a.y - m1 * a.x;

    let m2 = (v1.y - v2.y) / (v1.x - v2.x);
    let t2 = v1.y - m2 * v1.x;

    let x0 = (t2 - t1) / (m1 - m2);
    let y0 = m1 * x0 + t1;

    Vertex { x: x0, y: y0 }
}

fn inside(v: &Vertex, v1: &Vertex, v2: &Vertex) -> bool {
    let a = Vector2::from_vec(vec![v1.x - v2.x, v1.y - v2.y]);
    let b = Vector2::from_vec(vec![v1.x - v.x, v1.y - v.y]);

    a.x * b.y - a.y * b.x >= 0.
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

        for j in 0..input.len() {
            let curr = input[j];
            let prev = input[(j + input.len() - 1) % input.len()];

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
