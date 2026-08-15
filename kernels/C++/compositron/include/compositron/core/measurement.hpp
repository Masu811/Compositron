#pragma once

#include <ostream>
#include <string>
#include <map>
#include <unordered_map>

#include "compositron/dbs/dbspectrum.hpp"
#include "compositron/cdbs/cdbspectrum.hpp"
#include "compositron/pals/palspectrum.hpp"


namespace compositron::core {


namespace measurement {


class Shape {
private:
    size_t dbs;
    size_t cdbs;
    size_t pals;

public:
    Shape() = delete;
    Shape(size_t dbs, size_t cdbs, size_t pals);

    friend std::ostream& operator<<(std::ostream& os, Shape const& shape);
};

std::ostream& operator<<(std::ostream& os, Shape const& shape);


} // namespace measurement


class Measurement {
public:
    std::string filename;
    std::string name;
    std::map<std::string, compositron::dbs::DBSpectrum> dbs;
    std::map<std::string, compositron::cdbs::CDBSpectrum> cdbs;
    std::map<std::string, compositron::pals::PALSpectrum> pals;
    std::unordered_map<std::string, std::string> metadata;

    Measurement() = default;

    Measurement(const std::string filename);

    measurement::Shape shape() const;
};


} // namespace compositron::core
