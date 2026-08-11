#include "compositron/core/measurement.hpp"

#include <ostream>
#include <format>
#include <string>

#include "compositron/importers/SLOPE_n42_importer.hpp"


namespace compositron::core {


namespace measurement {


// class Shape

Shape::Shape(size_t dbs, size_t cdbs) : dbs(dbs), cdbs(cdbs) { }

std::ostream& operator<<(std::ostream& os, Shape const& shape) {
    return os << std::format(
        "Shape(DBS: {}, CDBS: {})", shape.dbs, shape.cdbs
    );
}


} // namespace measurement


// class Measurement

Measurement::Measurement(const std::string filename) {
    *this = importers::import_SLOPE_n42(filename);
}

measurement::Shape Measurement::shape() const {
    return {
        dbs.size(),
        cdbs.size(),
    };
}


} // namespace compositron::core
