#include "compositron/core/measurement.hpp"
#include "compositron/pals/model.hpp"

int main() {
    compositron::core::Measurement m(
        "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat"
    );

    compositron::pals::model::LifetimeModel lt_model{{
        {"lifetime_1", {180., 100., 5000.}},
        {"intensity_1", {0.79, 0., 1.}},
        {"lifetime_2", {350., 100., 5000.}},
        {"intensity_2", {0.2, 0., 1.}},
        {"lifetime_3", {2700., 100., 5000.}},
        {"intensity_3", {0.01, 0., 1.}},
        {"res_fwhm_1", {210., 50., 1000.}},
        {"res_intensity_1", {0.8, 0., 1.}},
        {"res_fwhm_2", {300., 50., 1000.}},
        {"res_intensity_2", {0.2, 0., 1.}},
        {"res_t0_2", {0., -200., 200.}},
    }};

    auto& s = m.pals.at("A");

    s.fit(lt_model, 16500, 18000);

    std::cout << s.fit_report() << std::endl;

    return 0;
}
