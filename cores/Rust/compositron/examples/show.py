import matplotlib.pyplot as plt
import pandas as pd

num = pd.read_csv("./num.csv", header=None, index_col=False)
ana = pd.read_csv("./ana.csv", header=None, index_col=False)
py = pd.read_csv("./python.csv", header=None, index_col=False)

names = [
    "N",
    "t0",
    "lifetime_1",
    "intensity_1",
    "lifetime_2",
    "intensity_2",
    "lifetime_3",
    "intensity_3",
    "res_fwhm_1",
    "res_intensity_1",
    "res_fwhm_2",
    "res_intensity_2",
    "res_t0_2",
]

for i, name in enumerate(names):
    x = ana[i]
    y = num[i]
    z = py[i]

    fig, axs = plt.subplots(1, 2)

    fig.suptitle(name)

    axs[0].plot(x)
    axs[0].plot(y)

    axs[1].plot(z)

    plt.show()
