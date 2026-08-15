#include "compositron/importers/meps_dat_importer.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <optional>

#include "compositron/core/measurement.hpp"
#include "compositron/core/utils.hpp"
#include "compositron/dbs/dbspectrum.hpp"
#include "compositron/pals/palspectrum.hpp"


namespace compositron::importers {


namespace {


struct MePSFileBundle {
    std::string spectrum_file;
    std::string param_file;
    std::string start_phs_file;
    std::string stop_phs_file;
};


struct HeaderInfo {
    int index;
    std::string sample;
    double energy;
    double temperature;
    double bin_width;
};


std::vector<std::string> split_filename_components(std::string filename) {
    return filename
        | std::views::split('_')
        | std::ranges::to<std::vector<std::string>>();
}


MePSFileBundle get_all_filenames(std::string filepath) {
    auto path = std::filesystem::path(filepath);

    std::string filename = path.filename();
    std::filesystem::path dir = path.parent_path();

    auto filename_components = split_filename_components(filename);

    int n_components = filename_components.size();

    int critical_idx = n_components - 3;

    if (filename_components[critical_idx] == "PALS") {
        std::string begin, end;
        for (const auto i : std::views::iota(0, critical_idx)) {
            begin += filename_components[i] + "_";
        }
        for (const auto i : std::views::iota(critical_idx + 1, n_components)) {
            end += "_" + filename_components[i];
        }
        return MePSFileBundle {
            filepath,
            dir / (begin + "Measurement_Parameter" + end),
            dir / (begin + "PHS_Ch0" + end),
            dir / (begin + "PHS_Ch1" + end),
        };
    } else if (filename_components[critical_idx] == "Parameter") {
        std::string begin, end;
        for (const auto i : std::views::iota(0, critical_idx - 1)) {
            begin += filename_components[i] + "_";
        }
        for (const auto i : std::views::iota(critical_idx + 1, n_components)) {
            end += "_" + filename_components[i];
        }
        return MePSFileBundle {
            dir / (begin + "PALS" + end),
            filepath,
            dir / (begin + "PHS_Ch0" + end),
            dir / (begin + "PHS_Ch1" + end),
        };
    } else if (filename_components[critical_idx] == "Ch0") {
        std::string begin, end;
        for (const auto i : std::views::iota(0, critical_idx - 1)) {
            begin += filename_components[i] + "_";
        }
        for (const auto i : std::views::iota(critical_idx + 1, n_components)) {
            end += "_" + filename_components[i];
        }
        return MePSFileBundle {
            dir / (begin + "PALS" + end),
            dir / (begin + "Measurement_Parameter" + end),
            filepath,
            dir / (begin + "PHS_Ch1" + end),
        };
    } else if (filename_components[critical_idx] == "Ch1") {
        std::string begin, end;
        for (const auto i : std::views::iota(0, critical_idx - 1)) {
            begin += filename_components[i] + "_";
        }
        for (const auto i : std::views::iota(critical_idx + 1, n_components)) {
            end += "_" + filename_components[i];
        }
        return MePSFileBundle {
            dir / (begin + "PALS" + end),
            dir / (begin + "Measurement_Parameter" + end),
            dir / (begin + "PHS_Ch0" + end),
            filepath,
        };
    }

    throw std::runtime_error("Invalid filename");
}


core::utils::Spectrum import_spectrum_data(std::ifstream& spectrum_data) {
    std::vector<uint32_t> spectrum;
    uint32_t value;

    while (spectrum_data >> value) {
        spectrum.push_back(value);
    }

    return core::utils::Spectrum(spectrum);
}


double parse_header_field(std::string& field, std::string expected_unit) {
    int l = field.size();
    int u = expected_unit.size();

    double value;

    try {
        value = std::stod(field);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Could not parse spectrum file header");
    }

    if (field.compare(l-u, u, expected_unit) == 0) {
        return value;
    }

    throw std::runtime_error(
        std::format(
            "Unexpected unit {} in spectrum file header (expected {})",
            field.substr(l-u, u), expected_unit
        )
    );
}


HeaderInfo get_header_info(std::string header, core::Measurement& m) {
    auto header_fields = header
        | std::views::split('_')
        | std::views::transform([](auto field) {
            auto str = std::string(field.begin(), field.end());
            core::utils::trim(str);
            return str;
        })
        | std::ranges::to<std::vector<std::string>>();

    int header_size = header_fields.size();

    if (header_size < 6) {
        throw std::runtime_error("Spectrum header has unexpected format");
    }

    std::string sample;

    for (auto i : std::views::iota(1, header_size - 5)) {
        sample += header_fields[i] + "_";
    }

    if (!sample.empty()) {
        sample.pop_back();
    }

    double energy = parse_header_field(header_fields[header_size - 4], "keV");
    double temperature = parse_header_field(header_fields[header_size - 3], "K");
    double bin_width = parse_header_field(header_fields[header_size - 1], "ns");

    HeaderInfo header_info {
        std::stoi(header_fields[0]),
        sample,
        energy,
        temperature,
        bin_width
    };

    m.metadata["Sample"] = sample;
    m.metadata["PositronImplantationEnergy"] = std::format("{}", energy);
    m.metadata["SampleTemperature"] = std::format("{}", temperature);

    return header_info;
}


void compare_filename_and_header(std::string filepath, HeaderInfo header) {
    std::string filename = std::filesystem::path(filepath).filename();

    auto filename_components = split_filename_components(filename);

    int l = filename_components.size();

    if (l < 6) {
        // TODO: warn
        return;
    }

    if (std::stoi(filename_components[0]) != header.index) {
        // TODO: warn
    }

    std::string sample;

    for (auto i : std::views::iota(1, l - 5)) {
        sample += filename_components[i] + "_";
    }

    if (!sample.empty()) {
        sample.pop_back();
    }

    if (sample != header.sample) {
        // TODO: warn
    }

    double energy = parse_header_field(filename_components[l - 6], "keV");

    if (energy != header.energy) {
        // TODO: warn
    }

    double temperature = parse_header_field(filename_components[l - 5], "K");

    if (temperature != header.temperature) {
        // TODO: warn
    }
}


void import_palspectrum(std::string filepath, core::Measurement& m) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::string header;
    std::getline(file, header);

    auto header_info = get_header_info(header, m);

    compare_filename_and_header(filepath, header_info);

    auto spectrum = import_spectrum_data(file);

    core::utils::TimingDetectorPair detpair {
        "A",
        core::utils::LinearCalibration {
            0, 1000 * header_info.bin_width
        },
        std::nullopt,
        NAN
    };

    m.pals.emplace("A", pals::PALSpectrum(std::move(spectrum), detpair));
}


void import_params(std::string filepath, core::Measurement& m) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        // TODO: warn
        return;
    }

    std::string line;

    while (std::getline(file, line)) {
        core::utils::trim(line);

        if (line.size() < 30) continue;

        std::string key = line.substr(0, 29);
        std::string value = line.substr(29);

        core::utils::trim(key);
        core::utils::trim(value);

        m.metadata[key] = value;
    }
}


std::optional<dbs::DBSpectrum> import_ph_spectrum(std::string filepath) {
    std::ifstream file(filepath);

    if (!file.is_open()) return std::nullopt;

    auto spectrum = import_spectrum_data(file);

    core::utils::EnergyDetector detector {
        "",
        core::utils::LinearCalibration { NAN, NAN },
        std::nullopt,
        NAN
    };

    return dbs::DBSpectrum(std::move(spectrum), detector);
}


} // anonymous namespace


core::Measurement import_MePS_dat(const std::string& filepath) {
    auto datafiles = get_all_filenames(filepath);

    core::Measurement m{};

    m.name = std::filesystem::path(filepath).filename();

    import_palspectrum(datafiles.spectrum_file, m);

    import_params(datafiles.param_file, m);

    std::optional<dbs::DBSpectrum> d1 = import_ph_spectrum(datafiles.start_phs_file);

    if (d1.has_value()) {
        auto& d = d1.value();
        d.detector.name = "Start";
        m.dbs.emplace("Start", std::move(d));
    }

    std::optional<dbs::DBSpectrum> d2 = import_ph_spectrum(datafiles.stop_phs_file);

    if (d2.has_value()) {
        auto& d = d2.value();
        d.detector.name = "Stop";
        m.dbs.emplace("Stop", std::move(d));
    }

    return m;
}


} // namespace compositron::importers
