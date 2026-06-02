import numpy as np
import numpy.typing as npt


def to_min_type_spectrum(spectrum: npt.ArrayLike) -> npt.NDArray[np.unsignedinteger]:
    max = np.max(spectrum)

    if max <= np.iinfo(np.uint8).max:
        return np.asarray(spectrum, dtype=np.uint8)
    elif max <= np.iinfo(np.uint16).max:
        return np.asarray(spectrum, dtype=np.uint16)
    elif max <= np.iinfo(np.uint32).max:
        return np.asarray(spectrum, dtype=np.uint32)
    elif max <= np.iinfo(np.uint64).max:
        return np.asarray(spectrum, dtype=np.uint64)

    raise NotImplementedError(
        "Greater than 64 bit uint values are not supported"
    )
