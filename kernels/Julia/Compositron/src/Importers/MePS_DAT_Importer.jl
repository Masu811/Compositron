using ..Core
using ...Utils
using ...PALS
using ...DBS

export import_meps_dat


struct MePSFileBundle
    spectrum_file::String
    param_file::String
    start_phs_file::String
    stop_phs_file::String
end


struct HeaderInfo
    index::UInt32
    sample::String
    energy::Float64
    temperature::Float64
    bin_width::Float64
end


function split_filepath_components(filepath::String)::Vector{String}
    split(basename(filepath), "_")
end


function get_all_filenames(filepath::String)::MePSFileBundle
    filename_components = split_filepath_components(filepath)
    dir = dirname(filepath)

    critical_idx = length(filename_components) - 2

    if filename_components[critical_idx] == "PALS"
        start = join(filename_components[1:critical_idx-1], "_")
        stop = join(filename_components[critical_idx+1:end], "_")
        return MePSFileBundle(
            filepath,
            joinpath(dir, join([start, "Measurement_Parameter", stop], "_")),
            joinpath(dir, join([start, "PHS_Ch0", stop], "_")),
            joinpath(dir, join([start, "PHS_Ch1", stop], "_")),
        )
    elseif filename_components[critical_idx] == "Parameter"
        start = join(filename_components[1:critical_idx-2], "_")
        stop = join(filename_components[critical_idx+1:end], "_")
        return MePSFileBundle(
            joinpath(dir, join([start, "PALS", stop], "_")),
            filepath,
            joinpath(dir, join([start, "PHS_Ch0", stop], "_")),
            joinpath(dir, join([start, "PHS_Ch1", stop], "_")),
        )
    elseif filename_components[critical_idx] == "Ch0"
        start = join(filename_components[1:critical_idx-2], "_")
        stop = join(filename_components[critical_idx+1:end], "_")
        return MePSFileBundle(
            joinpath(dir, join([start, "PALS", stop], "_")),
            joinpath(dir, join([start, "Measurement_Parameter", stop], "_")),
            filepath,
            joinpath(dir, join([start, "PHS_Ch1", stop], "_")),
        )
    elseif filename_components[critical_idx] == "Ch1"
        start = join(filename_components[1:critical_idx-2], "_")
        stop = join(filename_components[critical_idx+1:end], "_")
        return MePSFileBundle(
            joinpath(dir, join([start, "PALS", stop], "_")),
            joinpath(dir, join([start, "Measurement_Parameter", stop], "_")),
            joinpath(dir, join([start, "PHS_Ch0", stop], "_")),
            filepath,
        )
    end

    throw(ErrorException("Invalid Filename"))
end


function parse_header_field(field::AbstractString, expected_unit::String)::Float64
    l = length(field)
    u = length(expected_unit)

    value = parse(Float64, field[1:l-u])
    unit = field[l-u+1:end]

    if unit != expected_unit
        throw(ErrorException(
            "Unexpected unit $unit in spectrum file header "
            * "(expected $expected_unit)"
        ))
    end

    value
end


function get_header_info!(header::String, m::Measurement)::HeaderInfo
    header_fields = header |> x -> split(x, "_") .|> strip

    header_size = length(header_fields)

    if header_size < 6
        throw(ErrorException("Spectrum header has unexpected format"))
    end

    header_info = HeaderInfo(
        parse(UInt32, header_fields[1]),
        join(header_fields[2:header_size-4], "_"),
        parse_header_field(header_fields[header_size-3], "keV"),
        parse_header_field(header_fields[header_size-2], "K"),
        parse_header_field(header_fields[header_size], "ns"),
    )

    m.metadata["Sample"] = header_info.sample
    m.metadata["PositronImplantationEnergy"] = "$(header_info.energy)"
    m.metadata["SampleTemperature"] = "$(header_info.temperature)"

    header_info
end


function compare_filename_and_header(filepath::String, header::HeaderInfo)
    filename = basename(filepath)

    filename_components = split_filepath_components(filename)

    l = length(filename_components)

    if parse(UInt32, filename_components[1]) != header.index
        # TODO: warn
    end

    if parse_header_field(filename_components[l-5], "keV") != header.energy
        # TODO: warn
    end

    if parse_header_field(filename_components[l-4], "K") != header.temperature
        # TODO: warn
    end

    if join(filename_components[2:l-7], "_") != header.sample
        # TODO: warn
    end
end


function import_spectrum_data(channel_data::Vector{String})::Vector{<:Unsigned}
    min_spectrum(channel_data .|> strip .|> x -> parse(UInt32, x))
end


function import_palspectrum!(filepath::String, m::Measurement)::PALSpectrum
    spectrum_file_lines = readlines(filepath)

    header = spectrum_file_lines[1]

    channel_data = spectrum_file_lines[2:end]

    header_info = get_header_info!(header, m)

    compare_filename_and_header(filepath, header_info)

    spectrum = import_spectrum_data(channel_data)

    detpair = TimingDetectorPair(
        "A",
        LinearCalibration(0., 1000 * header_info.bin_width),
        nothing,
        nothing,
    )

    m.pals["A"] = PALSpectrum(spectrum, detpair)
end


function import_params!(filepath::String, m::Measurement)
    if !isfile(filepath)
        # TODO: warn
        return
    end

    for line in readlines(filepath)
        length(line) < 30 && continue

        key = strip(line[1:29])
        value = strip(line[30:end])

        m.metadata[key] = value
    end
end


function import_ph_spectrum(filepath::String)::Union{Nothing, DBSpectrum}
    if !isfile(filepath)
        # TODO: warn
        println("File $filepath not found")
        return
    end

    spectrum_data = readlines(filepath)

    spectrum = import_spectrum_data(spectrum_data)

    detector = EnergyDetector(
        "",
        LinearCalibration(NaN, NaN),
        nothing,
        nothing,
    )

    DBSpectrum(spectrum, detector)
end


function import_meps_dat(filepath::String)::Measurement
    datafiles = get_all_filenames(filepath)

    m = Measurement()

    import_palspectrum!(datafiles.spectrum_file, m)

    import_params!(datafiles.param_file, m)

    d1 = import_ph_spectrum(datafiles.start_phs_file)

    if !isnothing(d1)
        d1.detector.name = "Start"
        m.dbs["Start"] = d1
    end

    d2 = import_ph_spectrum(datafiles.stop_phs_file)

    if !isnothing(d2)
        d2.detector.name = "Stop"
        m.dbs["Stop"] = d2
    end

    m
end
