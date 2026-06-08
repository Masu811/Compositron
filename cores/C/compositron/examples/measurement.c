#include "compositron/compositron.h"

#include <stdio.h>

int main() {
    Measurement *m = newMeasurement();

    if (m == NULL) return 1;

    SlopeN42ImportStatus status = import_SLOPE_n42(
        m, "../../../../testdata/depth-profile_Copper_0000.n42"
    );

    if (status != SUCCESS) {
        freeMeasurement(m);
        return status;
    }

    MeasurementShape shape = getMeasurementShape(m);

    fprintf(stderr, "Shape: dbs=%zu, cdbs=%zu\n", shape.dbs, shape.cdbs);

    freeMeasurement(m);

    return 0;
}
