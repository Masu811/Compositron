#include "compositron/cdbs/cdbspectrum.hpp"
#include "compositron/core/measurement.hpp"
#include "compositron/core/utils.hpp"
#include <chrono>
#include <fstream>

int main() {
    using compositron::core::utils::Unit;
    using compositron::core::utils::KEV;

    // Import

    compositron::core::Measurement m(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    );

    auto& c = m.cdbs.at("OAA x OAB");

    c.correct_ecal(compositron::core::utils::EcalCorrectionOrder::FIRST);

    // S parameter calculation

    compositron::cdbs::LineshapeParamDefinition ls_param {
        "S",
        {
            compositron::cdbs::DiagonalArea(
                Unit(2, KEV),
                Unit(1, KEV),
                Unit(0, KEV),
                Unit(0, KEV)
            )
        },
        {
            compositron::cdbs::DiagonalArea(
                Unit(10, KEV),
                Unit(1, KEV),
                Unit(0, KEV),
                Unit(0, KEV)
            )
        },
        false
    };

    auto start = std::chrono::steady_clock::now();
    c.calc_lineshape_param(ls_param);
    auto stop = std::chrono::steady_clock::now();

    auto s = c.lineshape_params.at("S");

    std::cout << "Calculated S paramter: " << s.val << " +/- " << s.err << std::endl;
    std::cout << "Calculation duration: "
              << std::chrono::duration<double, std::micro>(stop - start).count()
              << " µs" << std::endl;

    // Projection onto diagonal

    start = std::chrono::steady_clock::now();
    auto p = c.project_diagonal(
        Unit(0.1, KEV),
        Unit(2, KEV),
        Unit(10, KEV),
        false
    );
    stop = std::chrono::steady_clock::now();

    std::cout << "Projection duration: "
              << std::chrono::duration<double, std::micro>(stop - start).count()
              << " µs" << std::endl;

    // Write projection to file for inspection

    std::ofstream file("projection.csv");

    if (!file.is_open()) {
        std::cout << "Could not open output file" << std::endl;
        return 1;
    }

    for (const auto& x : p.spectrum) {
        file << x << "\n";
    }

    return 0;
}
