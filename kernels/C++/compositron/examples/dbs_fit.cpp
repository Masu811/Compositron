#include <iostream>

#include "compositron/core/measurement.hpp"
#include "compositron/dbs/dbspectrum.hpp"


int main() {
    compositron::core::Measurement m(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    );

    std::cout << m.shape() << std::endl;

    auto& d = m.dbs.at("OAB");

    d.default_analyze();

    auto s = d.lineshape_params["S"];
    auto w = d.lineshape_params["W"];
    auto v = d.lineshape_params["V/P"];
    auto p = d.lineshape_params["P/T"];

    std::cout << "  S:   " << s.val << " +/- " << s.err << std::endl;
    std::cout << "  W:   " << w.val << " +/- " << w.err << std::endl;
    std::cout << "  V/P: " << v.val << " +/- " << v.err << std::endl;
    std::cout << "  P/T: " << p.val << " +/- " << p.err << std::endl;

    return 0;
}
