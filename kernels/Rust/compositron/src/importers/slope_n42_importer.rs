use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::num::{ParseIntError, ParseFloatError};
use std::path::{Path, PathBuf};

use thiserror::Error;
use roxmltree;
use roxmltree::Node;

use crate::dbs::DBSpectrum;
use crate::cdbs::CDBSpectrum;
use crate::core::Measurement;
use crate::core::utils::{EnergyDetector, EnergyDetectorPair, LinearCalibration, Spectrum, Spectrum2D};
use crate::importers::png_importer;
use crate::spectrum2d_match;


#[derive(Debug, Error)]
pub enum N42FormatError {
    #[error("Could not find required node {node} in the expected place")]
    NodeNotFoundError {
        node: String,
    },

    #[error("Could not find required attribute {attr} in node {node}")]
    AttributeNotFoundError {
        attr: String,
        node: String,
    },

    #[error("Detectors of multiple spectra have the same name {name}")]
    DuplicateNameError {
        name: String,
    },

    #[error(
        "A spectrum is referencing a {node} node with id {value} for its \
        {repr}, but no such node could be found"
    )]
    DanglingReferenceError {
        node: String,
        attr: String,
        repr: String,
        value: String,
    },

    #[error("Node {node} node does not contain any data")]
    EmptyFieldError {
        node: String,
    },

    #[error("Spectrum data has invalid format")]
    SpectrumParseError {
        #[from]
        inner: ParseIntError,
    },

    #[error("Energy calibration has invalid format")]
    EcalParseError {
        #[from]
        inner: ParseFloatError,
    },

    #[error("Energy calibration has invalid format")]
    InvalidEcalFormatError,

    #[error("Error while importing PNG file")]
    PNGImportError{
        #[from]
        source: png_importer::ImportError,
    },

    #[error("N42 format not yet supported: {detail}")]
    NotImplementedError {
        detail: String,
    },
}


#[derive(Debug, Error)]
pub enum ImportError {
    #[error("Error while opening N42 file {path}")]
    FileSystemError {
        #[source]
        inner: std::io::Error,
        path: PathBuf,
    },

    #[error("Error while reading N42 file {path}")]
    XMLParseError {
        #[source]
        inner: roxmltree::Error,
        path: PathBuf,
    },

    #[error("Error while importing N42 file {path}")]
    InvalidN42FormatError {
        #[source]
        inner: N42FormatError,
        path: PathBuf,
    },
}


#[derive(Debug)]
enum N42Ref {
    Ecal,
    Det,
    Detpair,
}


impl N42Ref {
    fn node(&self) -> &str {
        match self {
            N42Ref::Ecal => "EnergyCalibration",
            N42Ref::Det => "RadDetectorInformation",
            N42Ref::Detpair => "RadDetectorCoincidencePair",
        }
    }

    fn attr(&self) -> &str {
        match self {
            N42Ref::Ecal => "energyCalibrationReference",
            N42Ref::Det => "radDetectorInformationReference",
            N42Ref::Detpair => "radDetectorInformationReference",
        }
    }

    fn repr(&self) -> &str {
        match self {
            N42Ref::Ecal => "energy calibration",
            N42Ref::Det => "detector name",
            N42Ref::Detpair => "detector pair name",
        }
    }
}


fn find_node<'a, 'input>(
    parent: &Node<'a, 'input>, tag: &str
) -> Result<Node<'a, 'input>, N42FormatError> {
    match parent.children().find(|&e| e.has_tag_name(tag)) {
        Some(node) => Ok(node),
        None => Err(N42FormatError::NodeNotFoundError {
            node: tag.to_string(),
        }),
    }
}


fn find_attr<'a, 'input>(
    node: &Node<'a, 'input>, attr: &str
) -> Result<&'a str, N42FormatError> {
    match node.attribute(attr) {
        Some(value) => Ok(value),
        None => Err(N42FormatError::AttributeNotFoundError {
            attr: attr.into(),
            node: node.tag_name().name().into(),
        }),
    }
}


fn add(key: &str, node: &Node, metadata: &mut HashMap<String, String>) {
    if let Some(text) = node.text() {
        metadata.insert(key.into(), text.into());
    }
}


fn add_version(node: &Node, metadata: &mut HashMap<String, String>) {
    let mut version: Option<&str> = None;
    let mut component : Option<&str> = None;

    for child in node.children() {
        if !node.is_element() { continue; }

        let child_tag = child.tag_name().name();

        match child_tag {
            "RadInstrumentComponentName" => {
                if let Some(text) = child.text() {
                    component = Some(text);
                }
            },
            "RadInstrumentComponentVersion" => {
                if let Some(text) = child.text() {
                    version = Some(text);
                }
            },
            _ => {}
        }
    }

    if let Some(v) = version && let Some(c) = component {
        metadata.insert(
            format!("Instrument {c} version"), v.into()
        );
    }
}


fn gobble_hardware<'a, 'input>(
    hardware_node: &Node<'a, 'input>,
    hardware: &mut HashMap<String, Node<'a, 'input>>
) {
    for node in hardware_node.children() {
        if !node.is_element() { continue; }

        let tag = node.tag_name().name();

        if tag == "RadInstrumentHardwareElement"
            && let Ok(id) = find_attr(&node, "id")
        {
            hardware.insert(id.into(), node);
        }
    }
}


fn gobble_instrument_information<'a, 'input>(
    info: &Node<'a, 'input>,
    metadata: &mut HashMap<String, String>,
    hardware: &mut HashMap<String, Node<'a, 'input>>,
) {
    for node in info.children() {
        if !node.is_element() { continue; }

        let tag = node.tag_name().name();

        match tag {
            "RadInstrumentManufacturerName" => {
                add("Instrument Manufacturer Name", &node, metadata);
            },
            "RadInstrumentModelName" => {
                add("Instrument Model Name", &node, metadata);
            },
            "RadInstrumentVersion" => {
                add_version(&node, metadata);
            },
            "RadInstrumentHardware" => {
                gobble_hardware(&node, hardware);
            },
            _ => {}
        }
    }
}


fn gobble_metadata(parent: &Node, metadata: &mut HashMap<String, String>) {
    for node in parent.children().filter(|&node| node.is_element()) {
        if let Some(text) = node.text() {
            metadata.insert(
                node.tag_name().name().to_string(), text.to_string()
            );
        }
    }
}


fn check_exported(creator_node: &Node) {
    let Some(creator_name) = creator_node.text() else { return; };

    if creator_name != "STACS" { return; }

    eprintln!(
        "Imported data has been created with STACS and may have been altered"
    );
}


fn sort_root_children<'a, 'input>(
    root: &Node<'a, 'input>,
    measurements: &mut Vec<Node<'a, 'input>>,
    metadata: &mut HashMap<String, String>,
    ecals: &mut HashMap<String, Node<'a, 'input>>,
    detectors: &mut HashMap<String, Node<'a, 'input>>,
    detpairs: &mut HashMap<String, Node<'a, 'input>>,
    hardware: &mut HashMap<String, Node<'a, 'input>>,
) {
    for node in root.children() {
        if !node.is_element() { continue; }

        let tag = node.tag_name().name();

        match tag {
            "RadMeasurement" => {
                measurements.push(node);
            },
            "EnergyCalibration" => {
                if let Ok(id) = find_attr(&node, "id") {
                    ecals.insert(id.into(), node);
                };
            },
            "RadDetectorInformation" => {
                if let Ok(id) = find_attr(&node, "id") {
                    detectors.insert(id.into(), node);
                };
            },
            "RadDetectorCoincidencePair" => {
                if let Ok(id) = find_attr(&node, "id") {
                    detpairs.insert(id.into(), node);
                };
            },
            "RadInstrumentInformation" => {
                gobble_instrument_information(&node, metadata, hardware);
            },
            "RadInstrumentDataCreatorName" => {
                check_exported(&node);
            },
            _ => {}
        }
    }
}


fn import_hardware_readout(
    readout_node: &Node,
    hardware: &HashMap<String, Node>,
    metadata: &mut HashMap<String, String>,
) -> Result<(), N42FormatError> {
    let hardware_id = find_attr(
        readout_node, "radInstrumentHardwareElementReference"
    )?;
    let Some(hardware_node) = hardware.get(hardware_id) else {
        return Err(N42FormatError::DanglingReferenceError {
            node: "Readout".into(),
            attr: "radInstrumentHardwareElementReference".into(),
            repr: "hardware readout".into(),
            value: hardware_id.into(),
        });
    };

    let hardware_name_node = find_node(
        hardware_node, "RadInstrumentHardwareElementName",
    )?;
    let Some(hardware_name) = hardware_name_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: "RadInstrumentHardwareElementName".into(),
        });
    };

    let set_value_node = find_node(readout_node, "Set")?;
    let Some(set_value) = set_value_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: "Set".into(),
        });
    };

    let is_value_node = find_node(readout_node, "Is")?;
    let Some(is_value) = is_value_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: "Is".into(),
        });
    };

    metadata.insert(
        format!("{hardware_name}:Set Value"), set_value.into()
    );
    metadata.insert(
        format!("{hardware_name}:Is Value"), is_value.into()
    );

    Ok(())
}


fn get_or_parse_detname(
    det_id: &str,
    detectors: &HashMap<String, Node>,
    parsed_detnames: &mut HashMap<String, String>,
) -> Result<String, N42FormatError> {
    if let Some(detname) = parsed_detnames.get(det_id) {
        return Ok(detname.clone());
    }

    let Some(det_node) = detectors.get(det_id) else {
        return Err(N42FormatError::DanglingReferenceError {
            node: N42Ref::Det.node().into(),
            attr: N42Ref::Det.attr().into(),
            repr: N42Ref::Det.repr().into(),
            value: det_id.into(),
        });
    };

    let detname_node = find_node(det_node, "RadDetectorName")?;

    let detname = match detname_node.text() {
        Some(text) => text,
        None => return Err(N42FormatError::EmptyFieldError {
            node: det_node.tag_name().name().into(),
        }),
    };

    parsed_detnames.insert(det_id.into(), detname.into());

    Ok(detname.into())
}


fn get_or_parse_ecal(
    ecal_id: &str,
    ecals: &HashMap<String, Node>,
    parsed_ecals: &mut HashMap<String, LinearCalibration>,
) -> Result<LinearCalibration, N42FormatError> {
    if let Some(ecal) = parsed_ecals.get(ecal_id) {
        return Ok(ecal.clone());
    }

    let Some(ecal_node) = ecals.get(ecal_id) else {
        return Err(N42FormatError::DanglingReferenceError {
            node: N42Ref::Ecal.node().into(),
            attr: N42Ref::Ecal.attr().into(),
            repr: N42Ref::Ecal.repr().into(),
            value: ecal_id.into(),
        });
    };

    let ecal_value_node = find_node(&ecal_node, "CoefficientValues")?;

    let Some(values) = ecal_value_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: ecal_value_node.tag_name().name().into(),
        });
    };

    let values = values.split_whitespace();

    let mut ecal: Vec<f64> = Vec::new();

    for x in values {
        let x = x.parse::<f64>()?;
        ecal.push(x);
    }

    if ecal.len() < 2 {
        return Err(N42FormatError::InvalidEcalFormatError);
    }

    let ecal = LinearCalibration { offset: ecal[0], scale: ecal[1] };

    parsed_ecals.insert(ecal_id.into(), ecal);

    Ok(ecal)
}


fn parse_spectrum(spectrum_node: &Node) -> Result<Spectrum, N42FormatError> {
    let channel_data_node = find_node(spectrum_node, "ChannelData")?;

    let Some(channel_data) = channel_data_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: channel_data_node.tag_name().name().into(),
        });
    };

    let channel_data = channel_data.split_whitespace();

    let mut spectrum: Vec<u32> = Vec::new();

    for x in channel_data {
        let x = x.parse::<u32>()?;
        spectrum.push(x);
    }

    Ok(Spectrum::from(spectrum))
}


fn import_dbspectrum(
    spectrum_node: &Node,
    detectors: &HashMap<String, Node>,
    ecals: &HashMap<String, Node>,
    parsed_detnames: &mut HashMap<String, String>,
    parsed_ecals: &mut HashMap<String, LinearCalibration>,
    m: &mut Measurement,
) -> Result<(), N42FormatError> {
    let det_id = find_attr(spectrum_node, N42Ref::Det.attr())?;
    let detname = get_or_parse_detname(det_id, detectors, parsed_detnames)?;

    if m.dbs.contains_key(&detname) {
        return Err(N42FormatError::DuplicateNameError {
            name: detname.into(),
        });
    }

    let ecal_id = find_attr(spectrum_node, N42Ref::Ecal.attr())?;
    let ecal = get_or_parse_ecal(ecal_id, ecals, parsed_ecals)?;

    let spectrum = parse_spectrum(spectrum_node)?;

    let detector = EnergyDetector {
        name: detname.clone(),
        ecal,
        corrected_ecal: None,
        eres: None,
    };

    m.dbs.insert(detname, DBSpectrum::new(spectrum, detector));

    Ok(())
}


fn parse_detpair(detpair_node: &Node) -> Result<String, N42FormatError> {
    let name_node = find_node(detpair_node, "RadDetectorName")?;

    match name_node.text() {
        Some(text) => Ok(text.to_string()),
        None => Err(N42FormatError::EmptyFieldError {
            node: name_node.tag_name().name().into(),
        }),
    }
}


fn parse_window(params: &Node, param: &str) -> Result<usize, N42FormatError> {
    let param_node = find_node(&params, param)?;

    let Some(param_text) = param_node.text() else {
        return Err(N42FormatError::EmptyFieldError { node: param.into() });
    };

    Ok(param_text.parse::<usize>()?)
}


fn get_window(
    detpair_node: &Node
) -> Result<((usize, usize), (usize, usize)), N42FormatError> {
    let params = find_node(detpair_node, "CoincidenceParameters")?;

    let offset_x = parse_window(&params, "OffsetX")?;
    let offset_y = parse_window(&params, "OffsetY")?;
    let dim_x = parse_window(&params, "DimensionX")?;
    let dim_y = parse_window(&params, "DimensionY")?;

    Ok(((offset_x, offset_x + dim_x), (offset_y, offset_y + dim_y)))
}


fn get_ecal(
    detpair_node: &Node,
    detectors: &HashMap<String, Node>,
    parsed_detnames: &mut HashMap<String, String>,
    m: &Measurement,
    det_node_name: &str,
) -> Result<LinearCalibration, N42FormatError> {
    let det_node = find_node(detpair_node, det_node_name)?;
    let det_id = find_attr(&det_node, N42Ref::Detpair.attr())?;
    let detname = get_or_parse_detname(det_id, detectors, parsed_detnames)?;

    match m.dbs.get(&detname) {
        Some(spectrum) => Ok(spectrum.detector.ecal),
        None => {
            Err(N42FormatError::DanglingReferenceError {
                node: "RadDetector1Name".into(),
                attr: N42Ref::Detpair.attr().into(),
                repr: "detector name".into(),
                value: detname,
            })
        }
    }
}


fn get_or_parse_c_ecal(
    detpair_node: &Node,
    detectors: &HashMap<String, Node>,
    parsed_detnames: &mut HashMap<String, String>,
    m: &Measurement,
) -> Result<(LinearCalibration, LinearCalibration), N42FormatError> {
    Ok((
        get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector1Name"
        )?,
        get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector2Name"
        )?,
    ))
}


fn parse_cdbspectrum(
    spectrum_node: &Node,
    path: &Path
) -> Result<Spectrum2D, N42FormatError> {
    let Some(png_filename) = spectrum_node.text() else {
        return Err(N42FormatError::EmptyFieldError {
            node: spectrum_node.tag_name().name().into(),
        });
    };

    let directory = match path.parent() {
        Some(dir) => dir,
        None => Path::new(""),
    };

    let spectrum = png_importer::import_png(&directory.join(png_filename))?;

    Ok(Spectrum2D::from(spectrum))
}


fn import_cdbspectrum(
    spectrum_node: &Node,
    detpairs: &HashMap<String, Node>,
    detectors: &HashMap<String, Node>,
    parsed_detnames: &mut HashMap<String, String>,
    m: &mut Measurement,
    path: &Path,
) -> Result<(), N42FormatError> {
    let detpair_id = find_attr(
        &spectrum_node, N42Ref::Detpair.attr()
    )?;

    let Some(detpair_node) = detpairs.get(detpair_id) else {
        return Err(N42FormatError::DanglingReferenceError {
            node: N42Ref::Detpair.node().into(),
            attr: N42Ref::Detpair.attr().into(),
            repr: N42Ref::Detpair.repr().into(),
            value: detpair_id.into(),
        });
    };

    let detpair_name = parse_detpair(&detpair_node)?;

    if m.cdbs.contains_key(&detpair_name) {
        return Err(N42FormatError::DuplicateNameError {
            name: detpair_name.into(),
        });
    }

    let mut ecal = get_or_parse_c_ecal(
        detpair_node, detectors, parsed_detnames, m,
    )?;

    let window = get_window(&detpair_node)?;

    ecal.0.offset += window.0.0 as f64 * ecal.0.scale;
    ecal.1.offset += window.1.0 as f64 * ecal.1.scale;

    let det_1_node = find_node(detpair_node, "RadDetector1Name")?;
    let det_1_id = find_attr(&det_1_node, N42Ref::Detpair.attr())?;
    let det_1_name = get_or_parse_detname(det_1_id, detectors, parsed_detnames)?;

    let det_2_node = find_node(detpair_node, "RadDetector2Name")?;
    let det_2_id = find_attr(&det_2_node, N42Ref::Detpair.attr())?;
    let det_2_name = get_or_parse_detname(det_2_id, detectors, parsed_detnames)?;

    let spectrum = parse_cdbspectrum(spectrum_node, path)?;

    let first_det = EnergyDetector {
        name: det_1_name,
        ecal: ecal.0,
        corrected_ecal: None,
        eres: None,
    };

    let second_det = EnergyDetector {
        name: det_2_name,
        ecal: ecal.1,
        corrected_ecal: None,
        eres: None,
    };

    let detpair = EnergyDetectorPair {
        name: detpair_name.clone(),
        first_det,
        second_det,
        eres: None,
    };

    m.cdbs.insert(detpair_name, CDBSpectrum::new(spectrum, detpair));

    Ok(())
}


fn sort_meas_children<'a, 'input>(
    m: &mut Measurement,
    meas_node: &Node,
    ecals: &HashMap<String, Node>,
    detectors: &HashMap<String, Node>,
    detpairs: &HashMap<String, Node>,
    hardware: &HashMap<String, Node>,
    path: &Path,
) -> Result<(), N42FormatError> {

    let mut parsed_ecals: HashMap<String, LinearCalibration> = HashMap::new();
    let mut parsed_detnames: HashMap<String, String> = HashMap::new();
    let mut cdbs_nodes: Vec<Node> = Vec::new();

    for node in meas_node.children() {
        if !node.is_element() { continue; }

        let tag = node.tag_name().name();

        match tag {
            "StartDateTime" => {
                add("StartDateTime", &node, &mut m.metadata);
            },
            "RealTimeDuration" => {
                add("RealTimeDuration", &node, &mut m.metadata);
            },
            "Readout" => import_hardware_readout(
                &node, hardware, &mut m.metadata
            )?,
            "Spectrum" => {
                let Some(id) = node.attribute("id") else { continue; };

                if id.contains("Coinc") {
                    cdbs_nodes.push(node);
                } else {
                    import_dbspectrum(
                        &node,
                        detectors,
                        ecals,
                        &mut parsed_detnames,
                        &mut parsed_ecals,
                        m,
                    )?;
                }
            },
            _ => {
                if tag.contains("Metadata") {
                    gobble_metadata(&node, &mut m.metadata);
                } else if let Some(text) = node.text() {
                    m.metadata.insert(tag.into(), text.into());
                }
            }
        }
    }

    for node in cdbs_nodes.iter() {
        import_cdbspectrum(
            node, detpairs, detectors, &mut parsed_detnames, m, path,
        )?;
    }

    Ok(())
}


fn import_m(
    root: &Node, path: &Path
) -> Result<Measurement, N42FormatError> {
    let mut m = Measurement::new();

    let mut ecals: HashMap<String, Node> = HashMap::new();
    let mut detectors: HashMap<String, Node> = HashMap::new();
    let mut detpairs: HashMap<String, Node> = HashMap::new();
    let mut hardware: HashMap<String, Node> = HashMap::new();
    let mut measurements: Vec<Node> = Vec::new();

    sort_root_children(
        &root,
        &mut measurements,
        &mut m.metadata,
        &mut ecals,
        &mut detectors,
        &mut detpairs,
        &mut hardware,
    );

    if measurements.len() > 1 {
        return Err(N42FormatError::NotImplementedError {
            detail: "Multiple RadMeasurement nodes in one file".into(),
        })
    } else if measurements.len() == 0 {
        return Err(N42FormatError::NodeNotFoundError {
            node: "RadMeasurement".into(),
        });
    }

    let meas_node = measurements[0];

    sort_meas_children(
        &mut m, &meas_node, &ecals, &detectors, &detpairs, &hardware, path
    )?;

    Ok(m)
}


pub fn import_n42(filename: &str) -> Result<Measurement, ImportError> {
    let path = Path::new(filename);

    let mut file = File::open(&path).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: filename.into(),
        }
    )?;

    let mut s = String::new();

    file.read_to_string(&mut s).map_err(
        |err| ImportError::FileSystemError {
            inner: err,
            path: filename.into(),
        }
    )?;

    let doc = roxmltree::Document::parse(&s).map_err(
        |err| ImportError::XMLParseError {
            inner: err,
            path: filename.into(),
        }
    )?;

    let root = doc.root_element();

    import_m(&root, &path).map_err(
        |err| ImportError::InvalidN42FormatError {
            inner: err,
            path: filename.into(),
        }
    )
}
