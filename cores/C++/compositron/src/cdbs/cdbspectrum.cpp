#include "compositron/cdbs/cdbspectrum.hpp"

namespace compositron::cdbs {

using std::vector;
using std::string;
using std::tuple;
using compositron::core::utils::Spectrum2DPtr;

CDBSpectrum::CDBSpectrum(
    Spectrum2DPtr spectrum,
    string detpair,
    ecal_t ecal
) : spectrum(std::move(spectrum)), detpair(std::move(detpair)), ecal(ecal) { }

template<typename T>
CDBSpectrum::CDBSpectrum(
    vector<T>&& spectrum,
    std::tuple<size_t, size_t> shape,
    string detpair,
    ecal_t ecal
) : detpair(std::move(detpair)), ecal(ecal) {
    this->spectrum = core::utils::toMinTypeSpectrum2D(spectrum);
}

} // namespace compositron::dbs
