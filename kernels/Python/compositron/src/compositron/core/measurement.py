from ..dbs import DBSpectrum
from ..cdbs import CDBSpectrum


class Shape:
    def __init__(self, dbs: int, cdbs: int):
        self.dbs = dbs
        self.cdbs = cdbs

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f"Shape(dbs={self.dbs}, cdbs={self.cdbs})"


class Measurement:
    def __init__(self, filename: str | None = None):
        self.filename: str | None = filename
        self.name: str | None = filename
        self.dbs: dict[str, DBSpectrum] = {}
        self.cdbs: dict[str, CDBSpectrum] = {}
        self.metadata: dict[str, str] = {}

        if filename is not None:
            from ..importers.SLOPE_n42_importer import import_SLOPE_n42
            m = import_SLOPE_n42(filename)

            self.dbs = m.dbs
            self.cdbs = m.cdbs
            self.metadata = m.metadata

    @property
    def shape(self) -> Shape:
        return Shape(len(self.dbs), len(self.cdbs))
