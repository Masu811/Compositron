from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from pathlib import Path
import warnings

import numpy as np
import numpy.typing as npt

from ..dbs import DBSpectrum, Ecal as D_Ecal
from ..cdbs import CDBSpectrum, Ecal as C_Ecal
from ..core import Measurement
from .png_importer import import_png


class NodeNotFoundError(BaseException): pass
class AttributeNotFoundError(BaseException): pass
class DuplicateNameError(BaseException): pass
class DanglingReferenceError(BaseException): pass
class EmptyFieldError(BaseException): pass
class SpectrumParseError(BaseException): pass
class EcalParseError(BaseException): pass
class InvalidEcalFormatError(BaseException): pass


def add_version(
    node: Element[str],
    metadata: dict[str, str],
):
    version, component = None, None

    for child in node:
        match child.tag:
            case "RadInstrumentComponentName":
                if child.text is not None:
                    component = child.text
            case "RadInstrumentComponentVersion":
                if child.text is not None:
                    version = child.text

    if version is not None and component is not None:
        metadata[f"Instrument {component} version"] = version


def gobble_hardware(
    hardware_node: Element[str],
    hardware: dict[str, Element[str]],
):
    for node in hardware_node:
        tag = node.tag

        if tag == "RadInstrumentHardwareElement" and "id" in node.attrib:
            hardware[node.attrib["id"]] = node


def gobble_instrument_information(
    info: Element[str],
    metadata: dict[str, str],
    hardware: dict[str, Element[str]],
):
    for node in info:
        tag = node.tag

        match tag:
            case "RadInstrumentManufacturerName":
                value = node.text
                if value is not None:
                    metadata["Instrument Manufacturer Name"] = value
            case "RadInstrumentModelName":
                value = node.text
                if value is not None:
                    metadata["Instrument Model Name"] = value
            case "RadInstrumentVersion":
                add_version(node, metadata);
            case "RadInstrumentHardware":
                gobble_hardware(node, hardware);


def gobble_metadata(
    parent: Element[str],
    metadata: dict[str, str],
):
    for node in parent:
        value = node.text
        if value is not None:
            metadata[node.tag] = value


def check_exported(node: Element[str]):
    creator_name = node.text

    if creator_name is None or creator_name != "STACS":
        return

    warnings.warn(
        "Imported data has been created with STACS and may have been altered"
    )


def get_id(node: Element[str]):
    if not "id" in node.attrib:
        raise AttributeNotFoundError("id")

    return node.attrib["id"]


def sort_root_children(
    root: Element[str],
    measurements: list[Element[str]],
    metadata: dict[str, str],
    ecals: dict[str, Element[str]],
    detectors: dict[str, Element[str]],
    detpairs: dict[str, Element[str]],
    hardware: dict[str, Element[str]],
):
    for node in root:
        tag = node.tag

        match tag:
            case "RadMeasurement":
                measurements.append(node)
            case "EnergyCalibration":
                ecals[get_id(node)] = node
            case "RadDetectorInformation":
                detectors[get_id(node)] = node
            case "RadDetectorCoincidencePair":
                detpairs[get_id(node)] = node
            case "RadInstrumentInformation":
                gobble_instrument_information(node, metadata, hardware)
            case "RadInstrumentDataCreatorName":
                check_exported(node)


def import_hardware_readout(
    readout_node: Element[str],
    hardware: dict[str, Element[str]],
    metadata: dict[str, str],
):
    if not "radInstrumentHardwareElementReference" in readout_node.attrib:
        raise AttributeNotFoundError("radInstrumentHardwareElementReference")

    hardware_id = readout_node.attrib["radInstrumentHardwareElementReference"]

    if not hardware_id in hardware:
        raise DanglingReferenceError(hardware_id)

    hardware_node = hardware[hardware_id]

    hardware_name_node = hardware_node.find("RadInstrumentHardwareElementName")

    if hardware_name_node is None:
        raise NodeNotFoundError("RadInstrumentHardwareElementName")

    hardware_name = hardware_name_node.text

    if hardware_name is None:
        raise EmptyFieldError("RadInstrumentHardwareElementName")

    set_value_node = readout_node.find("Set")
    if set_value_node is None:
        raise NodeNotFoundError("Set")
    set_value = set_value_node.text
    if set_value is None:
        raise EmptyFieldError("Set")

    is_value_node = readout_node.find("Is")
    if is_value_node is None:
        raise NodeNotFoundError("Is")
    is_value = is_value_node.text
    if is_value is None:
        raise EmptyFieldError("Is")

    metadata[f"{hardware_name}:Set Value"] = set_value
    metadata[f"{hardware_name}:Is Value"] = is_value


def get_or_parse_detname(
    det_id: str,
    detectors: dict[str, Element[str]],
    parsed_detnames: dict[str, str],
) -> str:
    if det_id in parsed_detnames:
        return parsed_detnames[det_id]

    if not det_id in detectors:
        raise DanglingReferenceError(det_id)

    det_node = detectors[det_id]

    detname_node = det_node.find("RadDetectorName")

    if detname_node is None:
        raise NodeNotFoundError("RadDetectorName")

    detname = detname_node.text

    if detname is None:
        raise EmptyFieldError("RadDetectorName")

    detname = detname.strip()

    parsed_detnames[det_id] = detname

    return detname


def get_or_parse_ecal(
    ecal_id: str,
    ecals: dict[str, Element[str]],
    parsed_ecals: dict[str, D_Ecal],
) -> D_Ecal:
    if ecal_id in parsed_ecals:
        return parsed_ecals[ecal_id]

    if not ecal_id in ecals:
        raise DanglingReferenceError(ecal_id)

    ecal_node = ecals[ecal_id]

    ecal_value_node = ecal_node.find("CoefficientValues")

    if ecal_value_node is None:
        raise NodeNotFoundError("CoefficientValues")

    values = ecal_value_node.text

    if values is None:
        raise EmptyFieldError("EnergyCalibration CoefficientValues")

    coeffs = np.fromstring(values, sep=" ", dtype=np.float64)

    if len(coeffs) < 2:
        raise InvalidEcalFormatError()

    ecal = (coeffs[0], coeffs[1])

    parsed_ecals[ecal_id] = ecal

    return ecal


def parse_dbspectrum(spectrum_node: Element[str]) -> npt.NDArray[np.uint32]:
    channel_data_node = spectrum_node.find("ChannelData")

    if channel_data_node is None:
        raise NodeNotFoundError("ChannelData")

    channel_data = channel_data_node.text

    if channel_data is None:
        raise EmptyFieldError("ChannelData")

    spectrum = np.fromstring(channel_data, sep=" ", dtype=np.uint32)

    return spectrum


def import_dbspectrum(
    spectrum_node: Element[str],
    detectors: dict[str, Element[str]],
    ecals: dict[str, Element[str]],
    parsed_detnames: dict[str, str],
    parsed_ecals: dict[str, D_Ecal],
    m: Measurement,
):
    if not "radDetectorInformationReference" in spectrum_node.attrib:
        raise AttributeNotFoundError("radDetectorInformationReference")

    det_id = spectrum_node.attrib["radDetectorInformationReference"]

    detname = get_or_parse_detname(det_id, detectors, parsed_detnames)

    if detname in m.dbs:
        raise DuplicateNameError(detname)

    if not "energyCalibrationReference" in spectrum_node.attrib:
        raise AttributeNotFoundError("energyCalibrationReference")

    ecal_id  = spectrum_node.attrib["energyCalibrationReference"]

    ecal = get_or_parse_ecal(ecal_id, ecals, parsed_ecals)

    spectrum = parse_dbspectrum(spectrum_node)

    m.dbs[detname] = DBSpectrum(spectrum, detname, ecal)


def parse_detpair(detpair_node: Element[str]) -> str:
    name_node = detpair_node.find("RadDetectorName")

    if name_node is None:
        raise NodeNotFoundError("RadDetectorName")

    if name_node.text is None:
        raise EmptyFieldError("RadDetectorName")

    return name_node.text


def get_ecal(
    detpair_node: Element[str],
    detectors: dict[str, Element[str]],
    parsed_detnames: dict[str, str],
    m: Measurement,
    det_node_name: str,
) -> D_Ecal:
    det_node = detpair_node.find(det_node_name)

    if det_node is None:
        raise NodeNotFoundError(det_node_name)

    if not "radDetectorInformationReference" in det_node.attrib:
        raise AttributeNotFoundError("radDetectorInformationReference")

    det_id = det_node.attrib["radDetectorInformationReference"]

    detname = get_or_parse_detname(det_id, detectors, parsed_detnames)

    if not detname in m.dbs:
        raise DanglingReferenceError(detname)

    return m.dbs[detname].ecal


def get_or_parse_c_ecal(
    detpair_node: Element[str],
    detectors: dict[str, Element[str]],
    parsed_detnames: dict[str, str],
    m: Measurement,
) -> C_Ecal:
    return (
        get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector1Name",
        ), get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector1Name",
        ),
    )


def parse_cdbspectrum(
    spectrum_node: Element[str],
    path: Path
) -> npt.NDArray[np.uint32]:
    if spectrum_node.text is None:
        raise EmptyFieldError("CDB Spectrum node does not contain any data")

    png_filename = spectrum_node.text

    directory = path.parent

    return import_png(directory / png_filename)


def import_cdbspectrum(
    spectrum_node: Element[str],
    detpairs: dict[str, Element[str]],
    detectors: dict[str, Element[str]],
    parsed_detnames: dict[str, str],
    m: Measurement,
    path: Path,
):
    if not "radDetectorInformationReference" in spectrum_node.attrib:
        raise AttributeNotFoundError("radDetectorInformationReference")

    detpair_id = spectrum_node.attrib["radDetectorInformationReference"]

    if not detpair_id in detpairs:
        raise DanglingReferenceError(detpair_id)

    detpair_node = detpairs[detpair_id]

    detpair = parse_detpair(detpair_node)

    if detpair in m.cdbs:
        raise DuplicateNameError(detpair)

    ecal = get_or_parse_c_ecal(detpair_node, detectors, parsed_detnames, m)

    spectrum = parse_cdbspectrum(spectrum_node, path)

    m.cdbs[detpair] = CDBSpectrum(spectrum, detpair, ecal)


def sort_meas_children(
    m: Measurement,
    meas_node: Element[str],
    ecals: dict[str, Element[str]],
    detectors: dict[str, Element[str]],
    detpairs: dict[str, Element[str]],
    hardware: dict[str, Element[str]],
    path: Path,
):
    parsed_ecals: dict[str, D_Ecal] = {}
    parsed_detnames: dict[str, str] = {}
    cdbs_nodes: list[Element[str]] = []

    for node in meas_node:
        match node.tag:
            case "StartDateTime":
                if node.text is not None:
                    m.metadata["StartDateTime"] = node.text
            case "RealTimeDuration":
                if node.text is not None:
                    m.metadata["RealTimeDuration"] = node.text
            case "Readout":
                import_hardware_readout(node, hardware, m.metadata)
            case "Spectrum":
                id = get_id(node)

                if "Coinc" in id:
                    cdbs_nodes.append(node)
                else:
                    import_dbspectrum(
                        node,
                        detectors,
                        ecals,
                        parsed_detnames,
                        parsed_ecals,
                        m,
                    )
            case _:
                if "Metadata" in node.tag:
                    gobble_metadata(node, m.metadata)
                elif node.text is not None:
                    m.metadata[node.tag] = node.text

    for node in cdbs_nodes:
        import_cdbspectrum(node, detpairs, detectors, parsed_detnames, m, path)


def import_m(root: Element[str], path: Path) -> Measurement:
    m = Measurement()

    ecals: dict[str, Element[str]] = {}
    detectors: dict[str, Element[str]] = {}
    detpairs: dict[str, Element[str]] = {}
    hardware: dict[str, Element[str]] = {}
    measurements: list[Element[str]] = []

    sort_root_children(
        root, measurements, m.metadata, ecals, detectors, detpairs, hardware
    )

    if len(measurements) > 1:
        raise NotImplementedError("Multiple RadMeasurement nodes in one file")
    elif len(measurements) == 0:
        raise NodeNotFoundError("RadMeasurement")

    meas_node = measurements[0]

    sort_meas_children(
        m, meas_node, ecals, detectors, detpairs, hardware, path
    )

    return m


def import_SLOPE_n42(filename: str) -> Measurement:
    path = Path(filename)

    doc = ElementTree.parse(filename)

    root = doc.getroot()

    return import_m(root, path)
