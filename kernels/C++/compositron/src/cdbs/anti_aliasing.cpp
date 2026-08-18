#include "compositron/cdbs/anti_aliasing.hpp"

#include <cmath>
#include <iterator>
#include <utility>
#include <ranges>

#include "compositron/constants.hpp"

namespace compositron::cdbs::anti_aliasing {


// class Polygon

std::vector<Vertex> Polygon::to_vertices() {
    auto [x_0, y_0] = std::ranges::fold_left(
        vertices, std::pair(0., 0.), [](auto acc, auto vertex) {
            return std::pair(acc.first + vertex[0], acc.second + vertex[1]);
        }
    );

    x_0 /= vertices.size();
    y_0 /= vertices.size();

    std::vector<double> angles;
    std::transform(
        vertices.begin(), vertices.end(), std::back_inserter(angles),
        [x_0, y_0](const auto& v) {
            return std::atan2(v[1] - y_0, v[0] - x_0);
        }
    );

    auto pairs = std::views::zip(vertices, angles);

    std::ranges::sort(pairs, [](const auto& a, const auto& b) {
        return std::get<1>(a) < std::get<1>(b);
    });

    return vertices;
}


// class Rectangle

std::vector<Vertex> Rectangle::to_vertices() {
    return {{j_min, i_min}, {j_min, i_max}, {j_max, i_max}, {j_max, i_min}};
}


// class Parallelogram

std::vector<Vertex> Parallelogram::to_vertices() {
    double half_width_1 = side_length_1 / 2;
    double half_width_2 = side_length_2 / 2;

    Eigen::Vector2d w(1, slope_1);
    Eigen::Vector2d h(1, slope_2);

    w = w.normalized() * half_width_1;
    h = h.normalized() * half_width_2;

    Eigen::Vector2d x0(center_j, center_i);

    auto v1 = x0 + (w + h);
    auto v2 = x0 + (w - h);
    auto v3 = x0 + (-w - h);
    auto v4 = x0 + (-w + h);

    return {{v1[0], v1[1]}, {v2[0], v2[1]}, {v3[0], v3[1]}, {v4[0], v4[1]}};
}


// class Ellipse

std::vector<Vertex> Ellipse::to_vertices() {
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    std::vector<Vertex> vertices;

    size_t n_vertices = 128;

    for (size_t i = 0; i < n_vertices; ++i) {
        double t = 2 * PI * i / n_vertices;
        double s = std::sin(t);
        double c = std::cos(t);

        vertices.push_back(Vertex{
            center_j + radius_i * c * sin_phi + radius_j * s * cos_phi,
            center_i + radius_i * c * cos_phi - radius_j * s * sin_phi
        });
    }

    return vertices;
}


std::array<size_t, 4> get_bounding_box(std::vector<Vertex>& vertices) {
    double left_col_x = std::ranges::min(
        vertices | std::views::transform([](const auto& v) { return v[0]; })
    );
    double right_col_x = std::ranges::max(
        vertices | std::views::transform([](const auto& v) { return v[0]; })
    );
    double lower_row_y = std::ranges::max(
        vertices | std::views::transform([](const auto& v) { return v[1]; })
    );
    double upper_row_y = std::ranges::min(
        vertices | std::views::transform([](const auto& v) { return v[1]; })
    );

    size_t left_col = std::max(std::round(left_col_x), 0.);
    size_t right_col = std::max(std::round(right_col_x), 0.);
    size_t lower_row = std::max(std::round(lower_row_y), 0.);
    size_t upper_row = std::max(std::round(upper_row_y), 0.);

    return {left_col, right_col, lower_row, upper_row};
}


Vertex compute_intersection(Vertex& a, Vertex& b, Vertex& v1, Vertex& v2) {
    double m1 = (a[1] - b[1]) / (a[0] - b[0]);
    double t1 = a[1] - m1 * a[0];

    double m2 = (v1[1] - v2[1]) / (v1[0] - v2[0]);
    double t2 = v1[1] - m2 * v1[0];

    double x0 = (t2 - t1) / (m1 - m2);
    double y0 = m1 * x0 + t1;

    return Vertex{x0, y0};
}


bool inside(Vertex& v, Vertex& v1, Vertex& v2) {
    Eigen::Vector2d a(v1[0] - v2[0], v1[1] - v2[1]);
    Eigen::Vector2d b(v1[0] - v[0], v1[1] - v[1]);

    return a[0] * b[1] - a[1] * b[0] < 0;
}


}
