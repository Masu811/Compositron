#include <string>
#include <iostream>

#include <stacs/stacs.hpp>

int main() {
    std::string filename("../../../testdata/depth-profile_Copper_0000.n42");

    DopplerMeasurement m(filename);

    std::cout << m.shape() << "\n";

    return 0;
}
