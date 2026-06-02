from pathlib import Path

import numpy as np
import numpy.typing as npt
from PIL import Image


def import_png(filename: str | Path) -> npt.NDArray[np.uint32]:
    img = np.asarray(Image.open(filename).convert("RGB"), dtype=np.uint32)

    spectrum = (img[..., 0] << 16) | (img[..., 1] << 8) | img[..., 2]

    return spectrum.astype(np.uint32)
