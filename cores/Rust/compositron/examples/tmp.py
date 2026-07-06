from pypalsfit import LifetimeMeasurement

lt_model = {
    "N": (1e6, 0),
    "background": (0, False),
    "t0": (500),
    "lifetime_1": (180, 100, 5000),
    "intensity_1": (0.79, 0, 1),
    "lifetime_2": (350, 100, 5000),
    "intensity_2": (0.2, 0, 1),
    "lifetime_3": (2700, 100, 5000),
    "intensity_3": (0.01, 0, 1),
    "res_sigma_1": (90, 30, 500),
    "res_intensity_1": (0.8, 0, 1),
    "res_t0_1": (0, False),
    "res_sigma_2": (126, 30, 500),
    "res_intensity_2": (0.2, 0, 1),
    "res_t0_2": (0, -500, 500)
}

s = LifetimeMeasurement(
    "../../../../testdata/1_W3Re6_2.00keV_301.0K_10.0Mio_PALS_01.05.26_08.42.22.dat",
    lt_model=lt_model,
    verbose=True,
    calibrate=True,
    fit_start_idx=16500,
    fit_end_idx=17000,
)
