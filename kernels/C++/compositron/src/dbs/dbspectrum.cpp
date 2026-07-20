#include "compositron/dbs/dbspectrum.hpp"

namespace compositron::dbs {

using std::vector;
using std::string;
using std::tuple;

DBSpectrum::DBSpectrum(
    core::utils::SpectrumPtr spectrum,
    std::string detname,
    ecal_t ecal
) : spectrum(std::move(spectrum)), detname(std::move(detname)), ecal(ecal) {}

template<typename T>
DBSpectrum::DBSpectrum(vector<T>&& spectrum, string detname, ecal_t ecal)
: detname(std::move(detname)), ecal(ecal) {
    this->spectrum = core::utils::toMinTypeSpectrum(spectrum);
}

} // namespace compositron::dbs
