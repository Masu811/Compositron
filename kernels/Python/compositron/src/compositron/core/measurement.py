from pathlib import Path

from ..dbs import DBSpectrum
from ..cdbs import CDBSpectrum
from ..pals import PALSpectrum


class Shape:
    def __init__(self, dbs: int, cdbs: int, pals: int):
        self.dbs = dbs
        self.cdbs = cdbs
        self.pals = pals

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f"Shape(dbs={self.dbs}, cdbs={self.cdbs}, pals={self.pals})"


class Measurement:
    def __init__(self, filepath: str | None = None):
        self.filename: str | None = filepath
        self.name: str | None = filepath
        self.dbs: dict[str, DBSpectrum] = {}
        self.cdbs: dict[str, CDBSpectrum] = {}
        self.pals: dict[str, PALSpectrum] = {}
        self.metadata: dict[str, str] = {}

        if filepath is not None:
            if Path(filepath).suffix == ".n42":
                from ..importers.slope_n42_importer import import_slope_n42
                m = import_slope_n42(filepath)

            elif Path(filepath).suffix == ".dat":
                from ..importers.meps_dat_importer import import_meps_dat
                m = import_meps_dat(filepath)
            else:
                raise ValueError("Unsupported file format")

            self.dbs = m.dbs
            self.cdbs = m.cdbs
            self.pals = m.pals
            self.metadata = m.metadata

    @property
    def shape(self) -> Shape:
        return Shape(len(self.dbs), len(self.cdbs), len(self.pals))
