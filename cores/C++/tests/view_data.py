from stacs import DopplerMeasurement

m = DopplerMeasurement("../../../testdata/depth-profile_Copper_0000.n42", autocompute_singles=False)

print(m["OAB"].ecal)

m.show_singles()
