#include "compositron/core/measurement.hpp"

#include <ostream>
#include <format>
#include <string>
#include <filesystem>

#include "compositron/importers/meps_dat_importer.hpp"
#include "compositron/importers/slope_n42_importer.hpp"


namespace compositron::core {


namespace measurement {


// class Shape

Shape::Shape(
    size_t dbs, size_t cdbs, size_t pals
) : dbs(dbs), cdbs(cdbs), pals(pals) { }

std::ostream& operator<<(std::ostream& os, Shape const& shape) {
    return os << std::format(
        "Shape(DBS: {}, CDBS: {}, PALS: {})", shape.dbs, shape.cdbs, shape.pals
    );
}


} // namespace measurement


// class Measurement

Measurement::Measurement(const std::string filepath) {
    std::filesystem::path path(filepath);

    if (path.extension() == ".n42") {
        *this = importers::import_SLOPE_n42(filepath);
    } else if (path.extension() == ".dat") {
        *this = importers::import_MePS_dat(filepath);
    }

    filename = path.filename();
    name = filename;
}

measurement::Shape Measurement::shape() const {
    return {
        dbs.size(),
        cdbs.size(),
        pals.size(),
    };
}


} // namespace compositron::core
