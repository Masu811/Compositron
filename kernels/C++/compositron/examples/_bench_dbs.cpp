#include "compositron/cdbs/cdbspectrum.hpp"
#include "compositron/core/measurement.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/dbs/dbspectrum.hpp"
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

    std::cout << name << ": " << mean << " µs +/- " << stddev << " µs\n";
}


int main() {
    using compositron::core::utils::Unit;
    using compositron::core::utils::KEV;
    using compositron::dbs::Area;

    size_t n_iter = 10000;

    std::vector<double> imports(n_iter);
    std::vector<double> corr(n_iter);
    std::vector<double> bg_sub(n_iter);
    std::vector<double> s_calc(n_iter);

    for (size_t i = 0; i < n_iter; ++i) {
        auto start = std::chrono::steady_clock::now();
        compositron::core::Measurement m(
            "../../../../testdata/depth-profile_Copper_0000.n42"
        );
        auto stop = std::chrono::steady_clock::now();
        imports[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        auto& d = m.dbs.at("OAB");

        start = std::chrono::steady_clock::now();
        d.correct_ecal(compositron::core::utils::EcalCorrectionOrder::FIRST);
        stop = std::chrono::steady_clock::now();
        corr[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        start = std::chrono::steady_clock::now();
        d.extract_peak(30, compositron::dbs::PeakModel::ERFLINEAR2GAUSS);
        stop = std::chrono::steady_clock::now();
        bg_sub[i] = std::chrono::duration<double, std::micro>(stop - start).count();

        compositron::dbs::LineshapeParamDefinition ls_param {
            "S",
            {Area(Unit(2., KEV))},
            {Area(Unit(20., KEV))},
            false
        };

        start = std::chrono::steady_clock::now();
        d.calc_lineshape_param(ls_param, compositron::dbs::BinIntegrationScheme::CONST);
        stop = std::chrono::steady_clock::now();
        s_calc[i] = std::chrono::duration<double, std::micro>(stop - start).count();
    }

    print_mean_and_std("Import", imports);
    print_mean_and_std("Ecal Corr", corr);
    print_mean_and_std("BG Sub", bg_sub);
    print_mean_and_std("S Param", s_calc);

    return 0;
}
