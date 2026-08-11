#pragma once

#include <vector>


namespace compositron::cdbs::anti_aliasing {


using Vertex = double[2];

struct Polygon {
    std::vector<Vertex> vertices;
};


} // namespace compositron::cdbs::anti_aliasing
