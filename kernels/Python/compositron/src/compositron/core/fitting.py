from dataclasses import dataclass

import numpy as np


@dataclass
class SimpleFitParam:
    val: float
    err: float


@dataclass
class FitParam:
    val: float
    err: float = np.nan
    min: float = -np.inf
    max: float = np.inf
    vary: bool = True
    default_initialized: bool = False
