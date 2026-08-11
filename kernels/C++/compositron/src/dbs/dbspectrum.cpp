#include "compositron/dbs/dbspectrum.hpp"

#include <cmath>
#include <numeric>
#include <variant>

#include "compositron/core/utils.hpp"


namespace compositron::dbs {


using std::vector;
using std::string;
using std::tuple;


// class DBSpectrum

DBSpectrum::DBSpectrum(
    core::utils::Spectrum&& spectrum,
    core::utils::EnergyDetector detector
) : spectrum(spectrum), detector(detector) {
    counts = std::visit(
        [](const auto& spectrum) {
            return std::accumulate(
                spectrum.begin(),
                spectrum.end(),
                0,
                [](std::uint64_t acc, auto& x) { return acc + x; }
            );
        },
        this->spectrum.data
    );
    dcounts = std::sqrt(counts);
}


} // namespace compositron::dbs
