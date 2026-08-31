#include "compositron/cdbs/cdbspectrum.hpp"
#include "compositron/core/measurement.hpp"
#include "compositron/core/utils.hpp"
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>


void print_mean_and_std(std::string name, std::vector<double>& data) {
    const double mean =
        std::accumulate(data.begin(), data.end(), 0.0) / data.size();

    double sum_sq = 0.0;
    for (double x : data) {
        const double diff = x - mean;
        sum_sq += diff * diff;
    }

    const double stddev = std::sqrt(sum_sq / data.size());

    std::cout << name << ": " << mean << " ns +/- " << stddev << " ns\n";
}


int main() {
    using compositron::core::utils::Unit;
    using compositron::core::utils::KEV;

    size_t n_iter = 10000;

    std::vector<double> imports(n_iter);
    std::vector<double> corr(n_iter);
    std::vector<double> s_calc(n_iter);
    std::vector<double> proj_diag(n_iter);
    std::vector<double> proj_axes(n_iter);

    for (size_t i = 0; i < n_iter; ++i) {
        auto start = std::chrono::steady_clock::now();
        compositron::core::Measurement m(
            "../../../../testdata/depth-profile_Copper_0000.n42"
        );
        auto stop = std::chrono::steady_clock::now();
        imports[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        auto& c = m.cdbs.at("OAA x OAB");

        start = std::chrono::steady_clock::now();
        c.correct_ecal(compositron::core::utils::EcalCorrectionOrder::FIRST);
        stop = std::chrono::steady_clock::now();
        corr[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        compositron::cdbs::LineshapeParamDefinition ls_param {
            "S",
            {compositron::cdbs::DiagonalArea(Unit(2, KEV), Unit(1, KEV))},
            {compositron::cdbs::DiagonalArea(Unit(10, KEV), Unit(1, KEV))},
            false
        };

        start = std::chrono::steady_clock::now();
        c.calc_lineshape_param(ls_param);
        stop = std::chrono::steady_clock::now();
        s_calc[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        start = std::chrono::steady_clock::now();
        auto p = c.project_diagonal(
            Unit(0.1, KEV), Unit(2, KEV), Unit(100, KEV), false
        );
        stop = std::chrono::steady_clock::now();
        proj_diag[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        start = std::chrono::steady_clock::now();
        auto x = c.project_axes(compositron::cdbs::Axis::FIRST_DET);
        stop = std::chrono::steady_clock::now();
        proj_axes[i] = std::chrono::duration<double, std::micro>(stop - start).count();
    }

    print_mean_and_std("Import", imports);
    print_mean_and_std("Ecal Corr", corr);
    print_mean_and_std("S Param", s_calc);
    print_mean_and_std("Proj Diag", proj_diag);
    print_mean_and_std("Proj Axes", proj_axes);

    return 0;
}
