#pragma once

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include "compositron/cdbs/vendored/agg/agg_path_storage.h"
#include "compositron/cdbs/vendored/agg/agg_pixfmt_gray.h"
#include "compositron/cdbs/vendored/agg/agg_rasterizer_scanline_aa.h"
#include "compositron/cdbs/vendored/agg/agg_renderer_base.h"
#include "compositron/cdbs/vendored/agg/agg_rendering_buffer.h"
#include "compositron/cdbs/vendored/agg/agg_scanline_u.h"
#include "compositron/cdbs/vendored/agg/agg_renderer_scanline.h"


namespace compositron::cdbs::anti_aliasing {


using Vertex = std::array<double, 2>;


class ConvexToVerticesCounterClockwise {
public:
    virtual std::vector<Vertex> to_vertices() = 0;
};


class Polygon : public ConvexToVerticesCounterClockwise {
public:
    std::vector<Vertex> vertices;

    Polygon(std::vector<Vertex> vertices) : vertices(vertices) {}

    std::vector<Vertex> to_vertices() override;
};


class Rectangle : public ConvexToVerticesCounterClockwise {
public:
    double i_min;
    double i_max;
    double j_min;
    double j_max;

    Rectangle(
        double i_min, double i_max, double j_min, double j_max
    ) : i_min(i_min), i_max(i_max), j_min(j_min), j_max(j_max) {}

    std::vector<Vertex> to_vertices() override;
};


class Parallelogram : public ConvexToVerticesCounterClockwise {
public:
    double center_i;
    double center_j;
    double slope_1;
    double slope_2;
    double side_length_1;
    double side_length_2;

    Parallelogram(
        double center_i,
        double center_j,
        double slope_1,
        double slope_2,
        double side_length_1,
        double side_length_2
    ) : center_i(center_i),
        center_j(center_j),
        slope_1(slope_1),
        slope_2(slope_2),
        side_length_1(side_length_1),
        side_length_2(side_length_2) {}

    std::vector<Vertex> to_vertices() override;
};


class Ellipse : public ConvexToVerticesCounterClockwise {
public:
    double center_i;
    double center_j;
    double radius_i;
    double radius_j;
    double phi;

    Ellipse(
        double center_i,
        double center_j,
        double radius_i,
        double radius_j,
        double phi
    ) : center_i(center_i),
        center_j(center_j),
        radius_i(radius_i),
        radius_j(radius_j),
        phi(phi) {}

    std::vector<Vertex> to_vertices() override;
};


std::array<size_t, 4> get_bounding_box(std::vector<Vertex>& vertices);


template <typename T>
Eigen::MatrixX<uint8_t> agg_aa(size_t nrows, size_t ncols, T polygon) {
    if (nrows < 3 || ncols < 3) {
        throw std::runtime_error("Input data is too small in size");
    }

    std::vector<Vertex> vertices = polygon.to_vertices();

    if (vertices.size() < 3) {
        throw std::runtime_error("Insufficient number of vertices");
    }

    Eigen::MatrixX<uint8_t> weights(nrows, ncols);
    weights.setZero();

    std::transform(
        vertices.begin(), vertices.end(), vertices.begin(),
        [](const auto& v) { return Vertex{v[0] + 0.5, v[1] + 0.5}; }
    );

    if (std::any_of(vertices.begin(), vertices.end(), [nrows, ncols](const auto& v) {
        return v[0] < 0 || v[0] > ncols || v[1] < 0 || v[1] > nrows;
    })) {
        throw std::runtime_error("Polygon vertices are out of matrix bounds");
    }

    agg::path_storage path;
    path.move_to(vertices[0][1], vertices[0][0]);
    for (const auto& vertex : vertices) {
        path.line_to(vertex[1], vertex[0]);
    }
    path.close_polygon();

    agg::rendering_buffer rbuf(
        weights.data(),
        static_cast<unsigned>(nrows),
        static_cast<unsigned>(ncols),
        static_cast<int>(nrows)
    );

    agg::pixfmt_gray8 pixfmt(rbuf);
    agg::renderer_base<decltype(pixfmt)> ren_base(pixfmt);
    ren_base.clear(agg::gray8(0, 255));

    agg::rasterizer_scanline_aa<> ras;
    agg::scanline_u8 sl;
    ras.add_path(path);
    agg::render_scanlines_aa_solid(ras, sl, ren_base, agg::gray8(255, 255));

    return weights;
}


Vertex compute_intersection(Vertex& a, Vertex& b, Vertex& v1, Vertex& v2);


bool inside(Vertex& v, Vertex& v1, Vertex& v2);


template <typename T, typename U>
Polygon intersect_convex_polygons(T a, U b) {
    std::vector<Vertex> va = a.to_vertices();
    std::vector<Vertex> vb = b.to_vertices();

    if (va.size() < 3 || vb.size() < 3) {
        throw std::runtime_error("Insufficient number of vertices");
    }

    int la = va.size();

    std::vector<Vertex> output = vb;

    for (auto i : std::views::iota(0, la)) {
        Vertex edge_start = va[i];
        Vertex edge_end = va[(i + 1) % la];

        std::vector<Vertex> input = output;
        output = std::vector<Vertex>();

        int lb = input.size();

        for (auto j : std::views::iota(0, lb)) {
            Vertex curr = input[i];
            Vertex prev = input[(j + lb - 1) % lb];

            bool curr_inside = inside(curr, edge_start, edge_end);
            bool prev_inside = inside(prev, edge_start, edge_end);

            if (curr_inside) {
                if (!prev_inside) {
                    output.push_back(compute_intersection(
                        prev, curr, edge_start, edge_end
                    ));
                }
                output.push_back(curr);
            } else if (prev_inside) {
                output.push_back(compute_intersection(
                    prev, curr, edge_start, edge_end
                ));
            }
        }
    }

    return Polygon(output);
}


} // namespace compositron::cdbs::anti_aliasing
