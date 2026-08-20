using DataStructures: OrderedDict

import LightXML as xml

using ...DBS: DBSpectrum
using ...CDBS: CDBSpectrum
using ...Utils: min_spectrum, min_spectrum2d, LinearCalibration, EnergyDetector, EnergyDetectorPair
using ..Core: Measurement

function add_version!(
    node::xml.XMLElement,
    metadata::Dict{String, String},
)
    version = nothing; component = nothing

    component_node = xml.find_element(node, "RadInstrumentComponentName")
    if isnothing(component_node)
        return
    end

    version_node = xml.find_element(node, "RadInstrumentComponentVersion")
    if isnothing(version_node)
        return
    end

    component = xml.content(component_node)
    version = xml.content(version_node)

    metadata["Instrument $component version"] = version
end

function gobble_hardware!(
    hardware_node::xml.XMLElement,
    hardware::Dict{String, xml.XMLElement},
)
    for node in xml.child_elements(hardware_node)
        id = xml.attribute(node, "id")
        isnothing(id) && continue

        if xml.name(node) == "RadInstrumentHardwareElement"
            hardware[id] = node
        end
    end
end

function gobble_instrument_information!(
    info::xml.XMLElement,
    metadata::Dict{String, String},
    hardware::Dict{String, xml.XMLElement},
)
    for node in xml.child_elements(info)
        tag = xml.name(node)

        if tag == "RadInstrumentManufacturerName"
            metadata["Instrument Manufacturer Name"] = xml.content(node)
        elseif tag == "RadInstrumentModelName"
            metadata["Instrument Model Name"] = xml.content(node)
        elseif tag == "RadInstrumentVersion"
            add_version!(node, metadata)
        elseif tag == "RadInstrumentHardware"
            gobble_hardware!(node, hardware)
        end
    end
end

function gobble_metadata!(
    parent::xml.XMLElement, metadata::Dict{String, String}
)
    for node in xml.child_elements(parent)
        metadata[xml.name(node)] = xml.content(node)
    end
end

function check_exported(creator_node::xml.XMLElement)
    xml.content(creator_node) != "STACS" && return

    println(
        stderr,
        "Imported data has been created with STACS and may have been altered"
    )
end

function sort_root_children!(
    root::xml.XMLElement,
    measurements::Vector{xml.XMLElement},
    metadata::Dict{String, String},
    ecals::Dict{String, xml.XMLElement},
    detectors::Dict{String, xml.XMLElement},
    detpairs::Dict{String, xml.XMLElement},
    hardware::Dict{String, xml.XMLElement},
)
    for node in xml.child_elements(root)
        tag = xml.name(node)

        if tag == "RadMeasurement"
            push!(measurements, node)
        elseif tag == "EnergyCalibration"
            id = xml.attribute(node, "id")
            isnothing("id") && continue
            ecals[id] = node
        elseif tag == "RadDetectorInformation"
            id = xml.attribute(node, "id")
            isnothing("id") && continue
            detectors[id] = node
        elseif tag == "RadDetectorCoincidencePair"
            id = xml.attribute(node, "id")
            isnothing("id") && continue
            detpairs[id] = node
        elseif tag == "RadInstrumentInformation"
            gobble_instrument_information!(node, metadata, hardware)
        elseif tag == "RadInstrumentDataCreatorName"
            check_exported(node)
        end
    end
end

function import_hardware_readout!(
    readout_node::xml.XMLElement,
    hardware::Dict{String, xml.XMLElement},
    metadata::Dict{String, String},
)
    hardware_id = xml.attribute(
        readout_node, "radInstrumentHardwareElementReference"
    )

    if isnothing(hardware_id)
        throw(ErrorException(
            "Hardware readout node does not reference a hardware element node"
        ))
    end

    hardware_node = nothing

    try
        hardware_node = hardware[hardware_id]
    catch
        throw(ErrorException(
            "Hardware readout node references a hardware element "*
            "$(hardware_id), but no such element could be found"
        ))
    end

    hardware_name_node = xml.find_element(
        hardware_node, "RadInstrumentHardwareElementName"
    )

    if isnothing(hardware_name_node) || xml.content(hardware_name_node) == ""
        throw(ErrorException("Hardware element node does not contain a name"))
    end

    hardware_name = xml.content(hardware_name_node)

    set_value_node = xml.find_element(readout_node, "Set")
    if isnothing(set_value_node)
        throw(ErrorException(
            "Hardware readout node does not contain Set-value"
        ))
    end
    set_value = xml.content(set_value_node)

    is_value_node = xml.find_element(readout_node, "Is")
    if isnothing(set_value_node)
        throw(ErrorException("Hardware readout node does not contain Is-value"))
    end
    is_value = xml.content(is_value_node)

    metadata["$(hardware_name):Set Value"] = set_value
    metadata["$(hardware_name):Is Value"] = is_value
end

function get_or_parse_detname!(
    det_id::String,
    detectors::Dict{String, xml.XMLElement},
    parsed_detnames::Dict{String, String},
)::String
    det_id in keys(parsed_detnames) && return parsed_detnames[det_id]

    det_node = nothing

    try
        det_node = detectors[det_id]
    catch
        throw(ErrorException(
            "Could not find node for detector reference '$(det_id)'"
        ))
    end

    detname_node = xml.find_element(det_node, "RadDetectorName")

    if isnothing(detname_node)
        throw(ErrorException(
            "Detector at reference '$(det_id)' does not appear to have "*
            "a name specified"
        ))
    end

    detname = xml.content(detname_node)

    parsed_detnames[det_id] = detname

    detname
end

function get_or_parse_ecal!(
    ecal_id::String,
    ecals::Dict{String, xml.XMLElement},
    parsed_ecals::Dict{String, LinearCalibration},
)::LinearCalibration
    ecal_id in keys(parsed_ecals) && return parsed_ecals[ecal_id]

    ecal_node = nothing

    try
        ecal_node = ecals[ecal_id]
    catch
        throw(ErrorException(
            "Could not find node for energy calibration reference '$(ecal_id)'"
        ))
    end

    ecal_value_node = xml.find_element(ecal_node, "CoefficientValues")

    if isnothing(ecal_value_node) || xml.content(ecal_value_node) == ""
        throw(ErrorException(
            "Energy calibration node does not appear to contain coefficients"
        ))
    end

    values = (
        xml.content(ecal_value_node)
        |> data -> split(data, " ")
        .|> x -> parse(Float64, x)
    )

    if length(values) < 2
        throw(ErrorException(
            "Energy calibration coefficients appear to be in the wrong format"
        ))
    end

    ecal = LinearCalibration(values[1], values[2])

    # Ecal is given for 0-based indices, Julia uses 1-based indices
    ecal.offset -= ecal.scale

    parsed_ecals[ecal_id] = ecal

    ecal
end

function parse_spectrum(spectrum_node::xml.XMLElement)::Vector{<:Unsigned}
    channel_data_node = xml.find_element(spectrum_node, "ChannelData")
    if isnothing(channel_data_node) || xml.content(channel_data_node) == ""
        throw(ErrorException("Spectrum Node does not contain channel data"))
    end

    channel_data = xml.content(channel_data_node)

    spectrum = (
        channel_data
        |> data -> split(data, " ")
        .|> x -> parse(UInt, x)
    )

    min_spectrum(spectrum)
end

function import_dbspectrum!(
    spectrum_node::xml.XMLElement,
    detectors::Dict{String, xml.XMLElement},
    ecals::Dict{String, xml.XMLElement},
    parsed_detnames::Dict{String, String},
    parsed_ecals::Dict{String, LinearCalibration},
    m::Measurement
)
    det_id = xml.attribute(spectrum_node, "radDetectorInformationReference")
    if isnothing(det_id)
        throw(ErrorException(
            "Could not find attribute 'radDetectorInformationReference' "*
            "in node 'Spectrum'"
        ))
    end
    detname = get_or_parse_detname!(det_id, detectors, parsed_detnames)

    if detname in keys(m.dbs)
        throw(ErrorException(
            "Detectors of multiple spectra have the same name '$(detname)'"
        ))
    end

    ecal_id = xml.attribute(spectrum_node, "energyCalibrationReference")
    if isnothing(ecal_id)
        throw(ErrorException("Spectrum does not reference energy calibration"))
    end

    ecal = get_or_parse_ecal!(ecal_id, ecals, parsed_ecals)

    spectrum = parse_spectrum(spectrum_node)

    detector = EnergyDetector(detname, ecal, nothing, nothing)

    m.dbs[detname] = DBSpectrum(spectrum, detector)
end

function parse_detpair(detpair_node::xml.XMLElement)::String
    name_node = xml.find_element(detpair_node, "RadDetectorName")
    if isnothing(name_node) || xml.content(name_node) == ""
        throw(ErrorException(
            "Detector node does not contain a detector name"
        ))
    end
    xml.content(name_node)
end

function get_detector!(
    det_element_name::String,
    detpair_node::xml.XMLElement,
    detectors::Dict{String, xml.XMLElement},
    parsed_detnames::Dict{String, String},
    m::Measurement
)
    det_node = xml.find_element(detpair_node, det_element_name)
    if isnothing(det_node)
        throw(ErrorException(
            "Could not find node '$(det_element_name)' in the expected place"
        ))
    end
    det_id = xml.attribute(det_node, "radDetectorInformationReference")
    if isnothing(det_id)
        throw(ErrorException(
            "Could not find attribute 'radDetectorInformationReference' "*
            "in node '$(det_element_name)'"
        ))
    end
    detname = get_or_parse_detname!(det_id, detectors, parsed_detnames)
    try
        s = m.dbs[detname]
        return s.detector
    catch
        throw(ErrorException(
            "Coincidence detector pair is referencing detector $(detname), "*
            "but no such detector is known"
        ))
    end
end

function get_or_parse_c_dets!(
    detpair_node::xml.XMLElement,
    detectors::Dict{String, xml.XMLElement},
    parsed_detnames::Dict{String, String},
    m::Measurement
)::Tuple{EnergyDetector, EnergyDetector}
    (
        get_detector!(
            "RadDetector1Name", detpair_node, detectors, parsed_detnames, m,
        ), get_detector!(
            "RadDetector2Name", detpair_node, detectors, parsed_detnames, m,
        ),
    )
end

function parse_cdbspectrum(
    spectrum_node::xml.XMLElement,
    path::String,
)::Matrix{<:Unsigned}
    png_filename = xml.content(spectrum_node)

    if png_filename == ""
        throw(ErrorException(
            "CDBSpectrum node does not reference a PNG file"
        ))
    end

    directory = dirname(path)

    import_png(joinpath(directory, png_filename))
end

function import_cdbspectrum!(
    spectrum_node::xml.XMLElement,
    detectors::Dict{String, xml.XMLElement},
    detpairs::Dict{String, xml.XMLElement},
    parsed_detnames::Dict{String, String},
    m::Measurement,
    path::String,
)
    detpair_id = xml.attribute(spectrum_node, "radDetectorInformationReference")

    if isnothing(detpair_id)
        throw(ErrorException("Spectrum Node is missing detector reference"))
    end

    if !(detpair_id in keys(detpairs))
        throw(ErrorException(
            "A Spectrum node is referencing detector with id $(detpair_id), "*
            "but no such detector is known"
        ))
    end

    detpair_node = detpairs[detpair_id]
    detpair_name = parse_detpair(detpair_node)

    if detpair_name in keys(m.cdbs)
        throw(ErrorException(
            "Detectors of multiple spectra have the same name $detpair_name"
        ))
    end

    dets = get_or_parse_c_dets!(detpair_node, detectors, parsed_detnames, m)

    spectrum = parse_cdbspectrum(spectrum_node, path)

    detpair = EnergyDetectorPair(detpair_name, dets[1], dets[2], nothing)

    m.cdbs[detpair_name] = CDBSpectrum(spectrum, detpair)
end

function sort_meas_children!(
    m::Measurement,
    meas_node::xml.XMLElement,
    ecals::Dict{String, xml.XMLElement},
    detectors::Dict{String, xml.XMLElement},
    detpairs::Dict{String, xml.XMLElement},
    hardware::Dict{String, xml.XMLElement},
    filename::String,
)
    parsed_ecals = Dict{String, LinearCalibration}()
    parsed_detnames = Dict{String, String}()

    cdbs_nodes = Vector{xml.XMLElement}()

    for node in xml.child_elements(meas_node)
        tag = xml.name(node)

        if tag == "StartDateTime"
            m.metadata["StartDateTime"] = xml.content(node)
        elseif tag == "RealTimeDuration"
            m.metadata["RealTimeDuration"] = xml.content(node)
        elseif tag == "Readout"
            import_hardware_readout!(node, hardware, m.metadata)
        elseif tag == "Spectrum"
            id = xml.attribute(node, "id")
            isnothing(id) && continue
            if contains(id, "Coinc")
                push!(cdbs_nodes, node)
            else
                import_dbspectrum!(
                    node, detectors, ecals, parsed_detnames, parsed_ecals, m
                )
            end
        elseif contains(tag, "Metadata")
            gobble_metadata!(node, m.metadata)
        elseif xml.content(node) != ""
            m.metadata[tag] = xml.content(node)
        end
    end

    for node in cdbs_nodes
        import_cdbspectrum!(
            node, detectors, detpairs, parsed_detnames, m, filename,
        )
    end
end

function import_SLOPE_n42(filename::String)::Measurement
    doc = xml.parse_file(filename)

    root = xml.root(doc)

    m = Measurement()

    ecals = Dict{String, xml.XMLElement}()
    detectors = Dict{String, xml.XMLElement}()
    detpairs = Dict{String, xml.XMLElement}()
    hardware = Dict{String, xml.XMLElement}()
    measurements = Vector{xml.XMLElement}()

    sort_root_children!(
        root, measurements, m.metadata, ecals, detectors, detpairs, hardware,
    )

    if length(measurements) > 1
        throw(ErrorException("Multiple RadMeasurement nodes in one file"))
    elseif length(measurements) == 0
        throw(ErrorException("No RadMeasurement found in file"))
    end

    meas_node = measurements[1]

    sort_meas_children!(
        m, meas_node, ecals, detectors, detpairs, hardware, filename,
    )

    xml.free(doc)

    m
end
