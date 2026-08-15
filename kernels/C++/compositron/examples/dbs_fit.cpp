#include "compositron/core/measurement.hpp"
#include "compositron/dbs/dbspectrum.hpp"

#include <iostream>

int main() {
    compositron::core::Measurement m(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    );

    std::cout << m.shape() << std::endl;

    auto& s = m.dbs.at("OAB");

    s.default_analyze();

    std::cout << "S: "
              << s.lineshape_params["S"].val << " +/- "
              << s.lineshape_params["S"].err << std::endl;

    std::cout << "W: "
              << s.lineshape_params["W"].val << " +/- "
              << s.lineshape_params["W"].err << std::endl;

    std::cout << "V/P: "
              << s.lineshape_params["V/P"].val << " +/- "
              << s.lineshape_params["V/P"].err << std::endl;

    std::cout << "P/T: "
              << s.lineshape_params["P/T"].val << " +/- "
              << s.lineshape_params["P/T"].err << std::endl;

    return 0;
}
