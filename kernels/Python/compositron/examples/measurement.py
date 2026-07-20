import sys
sys.path.append("../src")
from compositron.core import Measurement

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

print(m.shape)
