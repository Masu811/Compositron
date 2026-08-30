#include "compositron/importers/slope_n42_importer.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <format>
#include <stdexcept>
#include <iostream>
#include <utility>

#include "compositron/core/measurement.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/dbs/dbspectrum.hpp"
#include "compositron/cdbs/cdbspectrum.hpp"
#include "compositron/importers/png_importer.hpp"
#include "compositron/importers/vendored/rapidxml.hpp"
#include "compositron/importers/vendored/rapidxml_utils.hpp"


namespace {


namespace xml = rapidxml;
namespace fs = std::filesystem;

using namespace compositron;

using std::string;
using std::vector;
using std::unordered_map;

using compositron::core::Measurement;


void add_version(
    xml::xml_node<>* readout_node,
    unordered_map<string, string>& metadata
) {
    auto component_node = readout_node->first_node(
        "RadInstrumentComponentName"
    );

    if (component_node == nullptr) return;

    auto version_node = readout_node->first_node(
        "RadInstrumentComponentVersion"
    );

    if (version_node == nullptr) return;

    string c = component_node->value();
    string v = version_node->value();

    metadata[std::format("Instrument {} version", c)] = v;
}


string get_id(xml::xml_node<>* node, string name) {
    auto id_attr = node->first_attribute("id");

    if (id_attr == nullptr) {
        throw std::runtime_error(std::format(
            "{} node does not have an id", name
        ));
    }

    return id_attr->value();
}


void gobble_hardware(
    xml::xml_node<>* hardware_node,
    unordered_map<string, xml::xml_node<>*>& hardware
) {
    for (
        auto node = hardware_node->first_node();
        node != nullptr;
        node = node->next_sibling()
    ) {
        if (node->type() != xml::node_type::node_element) continue;

        string tag = node->name();

        if (tag != "RadInstrumentHardwareElement") continue;

        string id = get_id(node, "Hardware element");

        hardware[id] = node;
    }
}


void gobble_instrument_information(
    xml::xml_node<>* info,
    unordered_map<string, string>& metadata,
    unordered_map<string, xml::xml_node<>*>& hardware
) {
    for (
        auto node = info->first_node();
        node != nullptr;
        node = node->next_sibling()
    ) {
        if (node->type() != xml::node_type::node_element) continue;

        string tag = node->name();

        if (tag == "RadInstrumentManufacturerName") {
            metadata["Instrument Manufacturer Name"] = node->value();
        } else if (tag == "RadInstrumentModelName") {
            metadata["Instrument Model Name"] = node->value();
        } else if (tag == "RadInstrumentVersion") {
            add_version(node, metadata);
        } else if (tag == "RadInstrumentHardware") {
            gobble_hardware(node, hardware);
        }
    }
}


void gobble_metadata(
    xml::xml_node<>* parent,
    unordered_map<string, string>& metadata
) {
    for (
        auto node = parent->first_node();
        node != nullptr;
        node = node->next_sibling()
    ) {
        if (node->type() != xml::node_type::node_element) continue;
        metadata[node->name()] = node->value();
    }
}


void check_exported(xml::xml_node<>* creator_node) {
    string creator_name = creator_node->value();

    if (creator_name != "STACS") return;

    std::cerr << "Warning: Imported data has been created with STACS "
              << "and may have been altered\n";
}


void sort_root_children(
    xml::xml_node<>* root,
    vector<xml::xml_node<>*>& measurements,
    unordered_map<string, string>& metadata,
    unordered_map<string, xml::xml_node<>*>& ecals,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, xml::xml_node<>*>& detpairs,
    unordered_map<string, xml::xml_node<>*>& hardware
) {
    for (
        auto node = root->first_node();
        node != nullptr;
        node = node->next_sibling()
    ) {
        if (node->type() != xml::node_type::node_element) continue;

        string tag = node->name();

        if (tag == "RadMeasurement") {
            measurements.push_back(node);
        } else if (tag == "EnergyCalibration") {
            string id = get_id(node, "Energy calibration");
            ecals[id] = node;
        } else if (tag == "RadDetectorInformation") {
            string id = get_id(node, "Detector information");
            detectors[id] = node;
        } else if (tag == "RadDetectorCoincidencePair") {
            string id = get_id(node, "Detector information");
            detpairs[id] = node;
        } else if (tag == "RadInstrumentInformation") {
            gobble_instrument_information(node, metadata, hardware);
        } else if (tag == "RadInstrumentDataCreatorName") {
            check_exported(node);
        }
    }
}


void import_hardware_readout(
    xml::xml_node<>* readout_node,
    unordered_map<string, xml::xml_node<>*>& hardware,
    unordered_map<string, string>& metadata
) {
    auto hardware_id_attr = readout_node->first_attribute(
        "radInstrumentHardwareElementReference"
    );
    if (hardware_id_attr == nullptr) {
        throw std::runtime_error(
            "Hardware readout node is not referencing a hardware element node"
        );
    }
    string hardware_id = hardware_id_attr->value();
    if (hardware.find(hardware_id) == hardware.end()) {
        throw std::runtime_error(std::format(
            "Hardware readout is referencing hardware element {} but no "
            "such element is known", hardware_id
        ));
    }

    auto hardware_node = hardware.at(hardware_id);
    auto hardware_name_node = hardware_node->first_node(
        "RadInstrumentHardwareElementName"
    );
    if (hardware_name_node == nullptr) {
        throw std::runtime_error(
            "Hardware element node does not contain a name"
        );
    }
    string hardware_name = hardware_name_node->value();

    auto set_value_node = readout_node->first_node("Set");
    auto is_value_node = readout_node->first_node("Is");

    if (set_value_node == nullptr) {
        throw std::runtime_error(
            "Hardware readout node does not contain Set-value"
        );
    }
    if (is_value_node == nullptr) {
        throw std::runtime_error(
            "Hardware readout node does not contain Is-value"
        );
    }

    string set_value = set_value_node->value();
    string is_value = is_value_node->value();
    metadata[std::format("{}:Set Value", hardware_name)] = set_value;
    metadata[std::format("{}:Is Value", hardware_name)] = is_value;
}


string get_or_parse_detname(
    string det_id,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, string> parsed_detnames
) {
    if (parsed_detnames.find(det_id) != parsed_detnames.end()) {
        return parsed_detnames.at(det_id);
    }

    if (detectors.find(det_id) == detectors.end()) {
        throw std::runtime_error(std::format(
            "Spectrum is referencing detector information node '{}' "
            "but no such detector is known", det_id
        ));
    }

    auto det_node = detectors.at(det_id);

    auto detname_node = det_node->first_node("RadDetectorName");

    if (detname_node == nullptr) {
        throw std::runtime_error(
            "Detector information does not contain detector name"
        );
    }

    string detname = detname_node->value();

    parsed_detnames[det_id] = detname;

    return detname;
}


core::utils::LinearCalibration get_or_parse_ecal(
    string ecal_id,
    unordered_map<string, xml::xml_node<>*>& ecals,
    unordered_map<string, core::utils::LinearCalibration> parsed_ecals
) {
    if (parsed_ecals.find(ecal_id) != parsed_ecals.end()) {
        return parsed_ecals.at(ecal_id);
    }

    if (ecals.find(ecal_id) == ecals.end()) {
        throw std::runtime_error(std::format(
            "Spectrum is referencing energy calibration node '{}' "
            "but no such calibration is known", ecal_id
        ));
    }

    auto ecal_node = ecals.at(ecal_id);

    auto ecal_value_node = ecal_node->first_node("CoefficientValues");

    if (ecal_value_node == nullptr) {
        throw std::runtime_error(
            "Energy calibration does not contain coefficient values"
        );
    }

    string ecal_values = string(ecal_value_node->value());

    std::istringstream ss(ecal_values);
    double a, b;

    if (!(ss >> a)) {
        throw std::runtime_error("Invalid energy calibration format");
    }

    if (!(ss >> b)) {
        throw std::runtime_error("Invalid energy calibration format");
    }

    core::utils::LinearCalibration ecal = {a, b};

    parsed_ecals[ecal_id] = ecal;

    return ecal;
}


core::utils::Spectrum parse_dbspectrum(xml::xml_node<>* spectrum_node) {
    auto channel_data_node = spectrum_node->first_node("ChannelData");

    string channel_data(channel_data_node->value());

    if (channel_data == "") {
        throw std::runtime_error("Spectrum does not contain any data");
    }

    std::istringstream ss(channel_data);
    uint32_t bin;
    vector<uint32_t> spectrum;

    int length = std::ranges::count(channel_data, ' ') + 1;
    spectrum.reserve(length);

    while (ss >> bin) {
        spectrum.push_back(bin);
    }

    return core::utils::Spectrum(spectrum);
}


void import_dbspectrum(
    xml::xml_node<>* spectrum_node,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, xml::xml_node<>*>& ecals,
    unordered_map<string, string> parsed_detnames,
    unordered_map<string, core::utils::LinearCalibration> parsed_ecals,
    Measurement& m
) {
    auto det_id = spectrum_node->first_attribute(
        "radDetectorInformationReference"
    );

    if (det_id == nullptr) {
        throw std::runtime_error(
            "Could not find attribute 'radDetectorInformationReference' "
            "in node 'Spectrum'"
        );
    }

    string detname = get_or_parse_detname(
        det_id->value(), detectors, parsed_detnames
    );

    if (m.dbs.contains(detname)) {
        throw std::runtime_error(std::format(
            "Detectors of multiple spectra have the same name {}", detname
        ));
    }

    auto ecal_id = spectrum_node->first_attribute("energyCalibrationReference");

    if (ecal_id == nullptr) {
        throw std::runtime_error(
            "Could not find attribute 'energyCalibrationReference' "
            "in node 'Spectrum'"
        );
    }

    core::utils::LinearCalibration ecal = get_or_parse_ecal(
        ecal_id->value(), ecals, parsed_ecals
    );

    core::utils::Spectrum spectrum = parse_dbspectrum(spectrum_node);

    core::utils::EnergyDetector detector = {
        detname, ecal, std::nullopt, std::numeric_limits<double>::quiet_NaN()
    };

    m.dbs.emplace(detname, dbs::DBSpectrum(std::move(spectrum), detector));
}


string parse_detpair(xml::xml_node<>* detpair_node) {
    auto name_node = detpair_node->first_node("RadDetectorName");

    if (name_node == nullptr) {
        throw std::runtime_error(
            "Detector node does not contain a detector name"
        );
    }

    return name_node->value();
}


size_t parse_window(xml::xml_node<>* params, string param) {
    auto param_node = params->first_node(param.data());

    if (param_node == nullptr) {
        throw std::runtime_error(std::format(
            "Could not find node '{}' in the expected place", param
        ));
    }

    return static_cast<size_t>(std::stol(param_node->value()));
}


std::array<std::array<size_t, 2>, 2> get_window(xml::xml_node<>* detpair_node) {
    auto params = detpair_node->first_node("CoincidenceParameters");

    if (params == nullptr) {
        throw std::runtime_error(
            "Could not find node 'CoincidenceParameters' in the expected place"
        );
    }

    size_t offset_x = parse_window(params, "OffsetX");
    size_t offset_y = parse_window(params, "OffsetY");
    size_t dim_x = parse_window(params, "DimensionX");
    size_t dim_y = parse_window(params, "DimensionY");

    return {{{offset_x, offset_x + dim_x}, {offset_y, offset_y + dim_y}}};
}


core::utils::EnergyDetector get_detector(
    xml::xml_node<>* detpair_node,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, string> parsed_detnames,
    Measurement& m,
    string det_node_name
) {
    auto det_node = detpair_node->first_node(det_node_name.data());

    if (det_node == nullptr) {
        throw std::runtime_error(std::format(
            "Could not find node '{}' in the expected place", det_node_name
        ));
    }

    auto det_id = det_node->first_attribute(
        "radDetectorInformationReference"
    );

    if (det_id == nullptr) {
        throw std::runtime_error(std::format(
            "Could not find attribute 'radDetectorInformationReference' "
            "in node '{}'", det_node_name
        ));
    }

    string detname = get_or_parse_detname(
        det_id->value(), detectors, parsed_detnames
    );

    if (m.dbs.find(detname) == m.dbs.end()) {
        throw std::runtime_error(std::format(
            "Coincidence detector pair is referencing detector {} "
            "but no such detector is known", detname
        ));
    }

    return m.dbs.at(detname).detector;
}


std::array<core::utils::EnergyDetector, 2> get_or_parse_c_dets(
    xml::xml_node<>* detpair_node,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, string> parsed_detnames,
    Measurement& m
) {
    auto det1 = get_detector(
        detpair_node, detectors, parsed_detnames, m, "RadDetector1Name"
    );
    auto det2 = get_detector(
        detpair_node, detectors, parsed_detnames, m, "RadDetector2Name"
    );

    auto window = get_window(detpair_node);

    det1.ecal.offset += window[0][0] * det1.ecal.scale;
    det2.ecal.offset += window[1][0] * det2.ecal.scale;

    return {det1, det2};
}


core::utils::Spectrum2D parse_cdbspectrum(
    xml::xml_node<>* spectrum_node,
    fs::path& path
) {
    string png_filename = spectrum_node->value();

    if (png_filename == "") {
        throw std::runtime_error(
            "CDBSpectrum node does not contain PNG filename"
        );
    }

    auto directory = path.parent_path();

    importers::Array2D spectrum = importers::import_png(
        (directory / png_filename).string()
    );

    return core::utils::Spectrum2D(
        spectrum.data,
        spectrum.height,
        spectrum.width,
        core::utils::StorageOrder::ROWMAJ
    );
}


void import_cdbspectrum(
    xml::xml_node<>* spectrum_node,
    unordered_map<string, xml::xml_node<>*>& detpairs,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, string> parsed_detnames,
    Measurement& m,
    fs::path& path
) {
    auto detpair_id = spectrum_node->first_attribute(
        "radDetectorInformationReference"
    );

    if (detpair_id == nullptr) {
        throw std::runtime_error("Spectrum node is missing detector reference");
    }

    if (!detpairs.contains(detpair_id->value())) {
        throw std::runtime_error(std::format(
            "Spectrum node is referencing detector with id '{}', but no "
            "such detector is known", detpair_id->value()
        ));
    }

    auto detpair_node = detpairs.at(detpair_id->value());

    string detpair_name = parse_detpair(detpair_node);

    if (m.cdbs.contains(detpair_name)) {
        throw std::runtime_error(std::format(
            "Detectors of multiple spectra have the same name {}", detpair_name
        ));
    }

    std::array<core::utils::EnergyDetector, 2> dets = get_or_parse_c_dets(
        detpair_node, detectors, parsed_detnames, m
    );

    core::utils::Spectrum2D spectrum = parse_cdbspectrum(
        spectrum_node, path
    );

    core::utils::EnergyDetectorPair detpair = {
        detpair_name,
        dets[0],
        dets[1],
        std::numeric_limits<double>::quiet_NaN()
    };

    m.cdbs.emplace(detpair_name, cdbs::CDBSpectrum(std::move(spectrum), detpair));
}


void sort_meas_children(
    Measurement& m,
    xml::xml_node<>* meas_node,
    unordered_map<string, xml::xml_node<>*>& ecals,
    unordered_map<string, xml::xml_node<>*>& detectors,
    unordered_map<string, xml::xml_node<>*>& detpairs,
    unordered_map<string, xml::xml_node<>*>& hardware,
    fs::path& path
) {
    unordered_map<string, core::utils::LinearCalibration> parsed_ecals;
    unordered_map<string, string> parsed_detnames;
    vector<xml::xml_node<>*> cdbs_nodes;

    for (
        auto node = meas_node->first_node();
        node != nullptr;
        node = node->next_sibling()
    ) {
        string tag = node->name();

        if (tag == "StartDateTime") {
            m.metadata["StartDateTime"] = node->value();
        } else if (tag == "RealTimeDuration") {
            m.metadata["RealTimeDuration"] = node->value();
        } else if (tag == "Readout") {
            import_hardware_readout(node, hardware, m.metadata);
        } else if (tag == "Spectrum") {
            string id = get_id(node, "Spectrum");
            if (id.contains("Coinc")) {
                cdbs_nodes.push_back(node);
            } else {
                import_dbspectrum(
                    node, detectors, ecals, parsed_detnames, parsed_ecals, m
                );
            }
        } else if (tag.contains("Metadata")) {
            gobble_metadata(node, m.metadata);
        } else if (string(node->value()) != "") {
            m.metadata[tag] = node->value();
        }
    }

    for (auto node : cdbs_nodes) {
        import_cdbspectrum(
            node, detpairs, detectors, parsed_detnames, m, path
        );
    }
}


Measurement import_m(xml::xml_node<>* root, fs::path& path) {
    Measurement m{};

    m.filename = path.string();
    m.name = m.filename;

    unordered_map<string, xml::xml_node<>*> ecals;
    unordered_map<string, xml::xml_node<>*> detectors;
    unordered_map<string, xml::xml_node<>*> detpairs;
    unordered_map<string, xml::xml_node<>*> hardware;
    vector<xml::xml_node<>*> measurements;

    sort_root_children(
        root, measurements, m.metadata, ecals, detectors, detpairs, hardware
    );

    if (measurements.size() > 1) {
        throw std::runtime_error("Multiple RadMeasurement nodes in one file");
    } else if (measurements.size() == 0) {
        throw std::runtime_error(
            "Could not find required node RadMeasurement in the expected place"
        );
    }

    auto meas_node = measurements[0];

    sort_meas_children(
        m, meas_node, ecals, detectors, detpairs, hardware, path
    );

    return m;
}


} // anonymous namespace


namespace compositron::importers {


Measurement import_SLOPE_n42(const std::string& filepath) {
    fs::path path = fs::path(filepath);
    fs::path directory = path.parent_path();

    xml::file<> file(filepath.c_str());

    xml::xml_document<> doc;
    doc.parse<0>(file.data());

    xml::xml_node<>* root = doc.first_node();

    return import_m(root, path);
}


} // namespace compositron::importers
