from compositron.core.measurement import Measurement

m = Measurement(
    "../../../../testdata/depth-profile_Copper_0000.n42"
)

print(m.shape)

s = m.dbs["OAB"]

s.default_analyze()

print(s.lineshape_params)
