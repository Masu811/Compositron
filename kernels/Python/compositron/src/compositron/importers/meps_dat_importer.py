from dataclasses import dataclass
from pathlib import Path

import numpy as np

from compositron.dbs.dbspectrum import DBSpectrum

from ..pals import PALSpectrum
from ..core.measurement import Measurement
from ..core.utils import EnergyDetector, LinearCalibration, Spectrum, TimingDetectorPair, to_min_type_spectrum


@dataclass
class MePSFileBundle:
    spectrum_file: Path
    param_file: Path
    start_phs_file: Path
    stop_phs_file: Path


@dataclass
class HeaderInfo:
    index: int
    sample: str
    energy: float
    temperature: float
    bin_width: float


def split_filepath_components(filepath: Path) -> list[str]:
    return filepath.name.split("_")


def get_all_filenames(filepath: Path) -> MePSFileBundle:
    filename_components = split_filepath_components(filepath)
    dir = filepath.parent

    critical_idx = len(filename_components) - 3

    match filename_components[critical_idx]:
        case "PALS":
            begin = "_".join(filename_components[:critical_idx])
            end = "_".join(filename_components[critical_idx+1:])
            return MePSFileBundle(
                dir / "_".join(filename_components),
                dir / "_".join([begin, "Measurement_Parameter", end]),
                dir / "_".join([begin, "PHS_Ch0", end]),
                dir / "_".join([begin, "PHS_Ch1", end]),
            )
        case "Parameter":
            begin = "_".join(filename_components[:critical_idx-1])
            end = "_".join(filename_components[critical_idx+1:])
            return MePSFileBundle(
                dir / "_".join([begin, "PALS", end]),
                dir / "_".join(filename_components),
                dir / "_".join([begin, "PHS_Ch0", end]),
                dir / "_".join([begin, "PHS_Ch1", end]),
            )
        case "PALS":
            begin = "_".join(filename_components[:critical_idx-1])
            end = "_".join(filename_components[critical_idx+1:])
            return MePSFileBundle(
                dir / "_".join([begin, "PALS", end]),
                dir / "_".join([begin, "Measurement_Parameter", end]),
                dir / "_".join(filename_components),
                dir / "_".join([begin, "PHS_Ch1", end]),
            )
        case "PALS":
            begin = "_".join(filename_components[:critical_idx-1])
            end = "_".join(filename_components[critical_idx+1:])
            return MePSFileBundle(
                dir / "_".join([begin, "PALS", end]),
                dir / "_".join([begin, "Measurement_Parameter", end]),
                dir / "_".join([begin, "PHS_Ch0", end]),
                dir / "_".join(filename_components),
            )
        case _:
            raise ValueError("Invalid filepath")


def parse_header_field(field: str, expected_unit: str) -> float:
    l = len(field)
    u = len(expected_unit)

    value = float(field[:l-u])
    unit = field[l-u:]

    if unit == expected_unit:
        return value

    raise ValueError(
        f"Unexpected unit {unit} in spectrum file header "
        f"(expected {expected_unit})"
    )


def get_header_info(header: str, m: Measurement) -> HeaderInfo:
    header_fields = [field.strip() for field in header.split("_")]

    header_size = len(header_fields)

    if header_size < 6:
        raise ValueError("Spectrum header has unexpected format")

    index = int(header_fields[0])
    bin_width = parse_header_field(header_fields[header_size-1], "ns")
    temperature = parse_header_field(header_fields[header_size-3], "K")
    energy = parse_header_field(header_fields[header_size-4], "keV")
    sample = "_".join(header_fields[1:header_size-5])

    m.metadata["Sample"] = sample
    m.metadata["PositronImplantationEnergy"] = f"{energy}"
    m.metadata["SampleTemperature"] = f"{temperature}"

    return HeaderInfo(index, sample, energy, temperature, bin_width)


def compare_filename_and_header(filepath: Path, header: HeaderInfo) -> None:
    filename_components = split_filepath_components(filepath)

    l = len(filename_components)

    if int(filename_components[0]) != header.index:
        # TODO: warn
        pass

    if parse_header_field(filename_components[l-6], "keV") != header.energy:
        # TODO: warn
        pass

    if parse_header_field(filename_components[l-5], "K") != header.temperature:
        # TODO: warn
        pass

    if "_".join(filename_components[1:l-7]) != header.sample:
        # TODO: warn
        pass


def import_spectrum_data(spectrum_data: list[str]) -> Spectrum:
    return to_min_type_spectrum(
        np.array([int(x.strip()) for x in spectrum_data])
    )


def import_palspectrum(filepath: Path, m: Measurement) -> PALSpectrum:
    with open(filepath, "r") as file:
        lines = file.readlines()

    header = lines[0]
    spectrum_data = lines[1:]

    header_info = get_header_info(header, m)

    compare_filename_and_header(filepath, header_info)

    spectrum = import_spectrum_data(spectrum_data)

    detpair = TimingDetectorPair(
        name = "A",
        tcal = LinearCalibration(0, 1000 * header_info.bin_width),
        corrected_tcal = None,
        tres = None,
    )

    return PALSpectrum(spectrum, detpair)


def import_params(filepath: Path, m: Measurement):
    with open(filepath, "r") as file:
        lines = file.readlines()

    for line in lines:
        if len(line) < 30:
            continue

        key = line[:29].strip()
        value = line[29:].strip()

        m.metadata[key] = value


def import_ph_spectrum(filepath: Path) -> None | DBSpectrum:
    with open(filepath, "r") as file:
        spectrum_data = file.readlines()

    spectrum = import_spectrum_data(spectrum_data)

    detector = EnergyDetector(
        name = "",
        ecal = LinearCalibration(np.nan, np.nan),
        corrected_ecal = None,
        eres = None,
    )

    return DBSpectrum(spectrum, detector)


def import_meps_dat(filepath: str) -> Measurement:
    datafiles = get_all_filenames(Path(filepath))

    m = Measurement()

    p = import_palspectrum(datafiles.spectrum_file, m)

    m.pals["A"] = p

    import_params(datafiles.param_file, m)

    d1 = import_ph_spectrum(datafiles.start_phs_file)
    if d1 is not None:
        d1.detector.name = "Start"
        m.dbs["Start"] = d1

    d2 = import_ph_spectrum(datafiles.stop_phs_file)
    if d2 is not None:
        d2.detector.name = "Stop"
        m.dbs["Stop"] = d2

    return m
