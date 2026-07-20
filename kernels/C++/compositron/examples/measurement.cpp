#include "compositron/core/measurement.hpp"

#include <iostream>

int main() {
    compositron::core::Measurement m(
        "../../../../testdata/depth-profile_Copper_0000.n42"
    );

    std::cout << m.shape() << std::endl;

    return 0;
}
