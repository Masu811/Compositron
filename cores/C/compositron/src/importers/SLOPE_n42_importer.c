#include "compositron/importers/SLOPE_n42_importer.h"

#include "compositron/core/utils.h"
#include "compositron/dbs/dbspectrum.h"
#include "compositron/cdbs/cdbspectrum.h"
#include "compositron/core/measurement.h"
#include "compositron/importers/png_importer.h"

#include <libxml2/libxml/parser.h>
#include <libxml2/libxml/tree.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


static xmlNodePtr findElementChild(xmlNodePtr parent, const char *name) {
    for (
        xmlNodePtr node = xmlFirstElementChild(parent);
        node != NULL;
        node = xmlNextElementSibling(node)
    ) {
        if (node->name != NULL && strcmp((char*)node->name, name) == 0) {
            return node;
        }
    }
    return NULL;
}

static char *make_version_key(char *component_name) {
    char *component_version_template = "Instrument %s version";

    size_t templ_size = strlen(component_version_template);
    size_t name_size = strlen(component_name);

    char *new = (char*)malloc((templ_size + name_size + 1) * sizeof(char));

    if (new == NULL) {
        fprintf(stderr, "Could not allocate memory for metadata entry\n");
        return NULL;
    }

    sprintf(new, component_version_template, component_name);

    return new;
}


static SlopeN42ImportStatus add_version(
    xmlNodePtr readout_node, Hashmap *metadata
) {
    SlopeN42ImportStatus status = SUCCESS;

    char *component_name = NULL, *version = NULL;
    char *key = NULL;

    xmlNodePtr component_node = findElementChild(
        readout_node, "RadInstrumentComponentName"
    );

    xmlNodePtr version_node = findElementChild(
        readout_node, "RadInstrumentComponentVersion"
    );

    if (component_node == NULL || version_node == NULL) {
        fprintf(stderr, "");
        status = NODE_NOT_FOUND_ERROR;
        return status;
    }

    component_name = (char*)xmlNodeGetContent(component_node);
    version = (char*)xmlNodeGetContent(version_node);

    if (component_name == NULL || version == NULL) {
        status = EMPTY_FIELD_ERROR;
        goto free;
    }

    key = make_version_key(component_name);

    if (key == NULL) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    if (!hashmap_insert(metadata, key, version)) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

free:

    if (component_name != NULL) free(component_name);
    if (status != SUCCESS) {
        if (key != NULL) free(key);
        if (version != NULL) free(version);
    }

    return status;
}


static SlopeN42ImportStatus gobble_hardware(
    xmlNodePtr hardware_node, Hashmap *hardware
) {
    SlopeN42ImportStatus status = SUCCESS;

    for (
        xmlNodePtr node = xmlFirstElementChild(hardware_node);
        node != NULL && status == SUCCESS;
        node = xmlNextElementSibling(node)
    ) {
        char *tag = (char*)node->name;

        if (tag == NULL || strcmp(tag, "RadInstrumentHardwareElement") != 0) {
            continue;
        }

        char *id = (char*)xmlGetProp(node, (xmlChar*)"id");

        if (id == NULL) {
            fprintf(stderr, "Hardware element node is missing an id\n");
            return ATTRIBUTE_NOT_FOUND_ERROR;
        }

        if (!hashmap_insert(hardware, id, node)) {
            free(id);
            return MEMORY_ALLOCATION_ERROR;
        }
    }

    return status;
}


static SlopeN42ImportStatus gobble_instrument_info_node(
    xmlNodePtr node, Hashmap *metadata, const char *key
) {
    char* content = (char*)xmlNodeGetContent(node);

    if (content == NULL) {
        return EMPTY_FIELD_ERROR;
    }

    size_t new_len = strlen(key);

    char *new_key = (char*)malloc((new_len + 1) * sizeof(char));

    if (new_key == NULL) {
        fprintf(stderr, "Could not allocate memory for metadata entry\n");
        return MEMORY_ALLOCATION_ERROR;
    }

    strcpy(new_key, key);

    new_key[new_len] = '\0';

    if (!hashmap_insert(metadata, new_key, content)) {
        free(content);
        return MEMORY_ALLOCATION_ERROR;
    }

    return SUCCESS;
}


static SlopeN42ImportStatus gobble_instrument_information(
    xmlNodePtr info, Hashmap *metadata, Hashmap *hardware
) {
    SlopeN42ImportStatus status = SUCCESS;

    for (
        xmlNodePtr node = xmlFirstElementChild(info);
        node != NULL && status == SUCCESS;
        node = xmlNextElementSibling(node)
    ) {
        char *tag = (char*)node->name;

        if (strcmp(tag, "RadInstrumentManufacturerName") == 0) {
            status = gobble_instrument_info_node(
                node, metadata, "Instrument Manufacturer Name"
            );
        } else if (strcmp(tag, "RadInstrumentModelName") == 0) {
            status = gobble_instrument_info_node(
                node, metadata, "Instrument Model Name"
            );
        } else if (strcmp(tag, "RadInstrumentVersion") == 0) {
            status = add_version(node, metadata);
        } else if (strcmp(tag, "RadInstrumentHardware") == 0) {
            status = gobble_hardware(node, hardware);
        }
    }

    return status;
}


static SlopeN42ImportStatus gobble_generic_text_node(
    xmlNodePtr node, Hashmap *metadata
) {
    char *content = (char*)xmlNodeGetContent(node);

    if (content == NULL) return SUCCESS;

    content = trim(content);

    if (strlen(content) == 0) {
        free(content);
        return SUCCESS;
    }

    char *tag_cpy = strdup((char*)node->name);

    if (tag_cpy == NULL) {
        fprintf(
            stderr, "Could not allocate memory for metadata entry\n"
        );
        free(content);
        return MEMORY_ALLOCATION_ERROR;
    }

    if (!hashmap_insert(metadata, tag_cpy, content)) {
        free(tag_cpy);
        free(content);
        return MEMORY_ALLOCATION_ERROR;
    }

    return SUCCESS;
}


static SlopeN42ImportStatus gobble_metadata(
    xmlNodePtr parent, Hashmap *metadata
) {
    SlopeN42ImportStatus status = SUCCESS;

    for (
        xmlNodePtr node = xmlFirstElementChild(parent);
        node != NULL && status == SUCCESS;
        node = xmlNextElementSibling(node)
    ) {
        status = gobble_generic_text_node(node, metadata);
    }

    return status;
}


static void check_exported(xmlNodePtr creator_node) {
    char *creator_name = (char*)xmlNodeGetContent(creator_node);

    if (creator_name == NULL) return;

    if (strcmp(creator_name, "STACS") == 0) {
        fprintf(
            stderr,
            "Warning: Imported data has been created with STACS and may have "
            "been altered\n"
        );
    }

    free(creator_name);
}


static SlopeN42ImportStatus gobble_meas_node(
    xmlNodePtr node, xmlNodePtr *measurements
) {
    if (*measurements != NULL) {
        fprintf(
            stderr,
            "Multiple RadMeasurements per file are presently not supported\n"
        );
        return NOT_IMPLEMENTED_ERROR;
    }

    *measurements = node;

    return SUCCESS;
}


static SlopeN42ImportStatus gobble_xml_node(
    xmlNodePtr node,
    Hashmap *map,
    const char *node_description
) {
    char *id = (char*)xmlGetProp(node, (xmlChar*)"id");

    if (id == NULL) {
        fprintf(stderr, "%s node is missing an id\n", node_description);
        return ATTRIBUTE_NOT_FOUND_ERROR;
    }

    if (!hashmap_insert(map, id, node)) {
        free(id);
        return MEMORY_ALLOCATION_ERROR;
    }

    return SUCCESS;
}


static SlopeN42ImportStatus sort_root_children(
    xmlNodePtr root,
    xmlNodePtr *measurements,
    Hashmap *metadata,
    Hashmap *ecals,
    Hashmap *detectors,
    Hashmap *detpairs,
    Hashmap *hardware
) {
    SlopeN42ImportStatus status = SUCCESS;

    for (
        xmlNodePtr node = xmlFirstElementChild(root);
        node != NULL && status == SUCCESS;
        node = xmlNextElementSibling(node)
    ) {
        char *tag = (char*)node->name;

        if (tag == NULL) continue;

        if (strcmp(tag, "RadMeasurement") == 0) {
            status = gobble_meas_node(node, measurements);
        } else if (strcmp(tag, "EnergyCalibration") == 0) {
            status = gobble_xml_node(node, ecals, "Energy Calibration");
        } else if (strcmp(tag, "RadDetectorInformation") == 0) {
            status = gobble_xml_node(node, detectors, "Detector Information");
        } else if (strcmp(tag, "RadDetectorCoincidencePair") == 0) {
            status = gobble_xml_node(node, detpairs, "Detector Information");
        } else if (strcmp(tag, "RadInstrumentInformation") == 0) {
            status = gobble_instrument_information(node, metadata, hardware);
        } else if (strcmp(tag, "RadInstrumentDataCreatorName") == 0) {
            check_exported(node);
        }
    }

    return status;
}


static SlopeN42ImportStatus import_hardware_readout(
    xmlNodePtr readout_node,
    Hashmap *hardware,
    Hashmap *metadata
) {
    SlopeN42ImportStatus status = SUCCESS;

    char *hardware_id = NULL, *hardware_name = NULL;
    char *set_value = NULL, *is_value = NULL;
    char *is_key = NULL, *set_key = NULL;

    hardware_id = (char*)xmlGetProp(
        readout_node, (xmlChar*)"radInstrumentHardwareElementReference"
    );

    if (hardware_id == NULL) {
        fprintf(
            stderr,
            "Hardware readout node is not referencing a hardware element node\n"
        );
        return ATTRIBUTE_NOT_FOUND_ERROR;
    }

    xmlNodePtr hardware_node = hashmap_get(hardware, hardware_id);

    if (hardware_node == NULL) {
        fprintf(
            stderr,
            "Hardware readout is referencing hardware element '%s', but no "
            "such element is known\n",
            hardware_id
        );
        status = DANGLING_REFERENCE_ERROR;
        goto free;
    }

    xmlNodePtr hardware_name_node = findElementChild(
        hardware_node, "RadInstrumentHardwareElementName"
    );

    if (hardware_name_node == NULL) {
        fprintf(
            stderr,
            "Hardware element node '%s' does not contain a name\n",
            hardware_id
        );
        status = NODE_NOT_FOUND_ERROR;
        goto free;
    }

    hardware_name = (char*)xmlNodeGetContent(hardware_name_node);

    if (hardware_name == NULL) {
        fprintf(
            stderr,
            "Hardware element node '%s' does not contain a name\n",
            hardware_id
        );
        status = EMPTY_FIELD_ERROR;
        goto free;
    }

    xmlNodePtr set_value_node = findElementChild(readout_node, "Set");
    xmlNodePtr is_value_node = findElementChild(readout_node, "Is");

    if (set_value_node == NULL || is_value_node == NULL) {
        fprintf(
            stderr,
            "Hardware readout node '%s' does not contain Set- or Is-value\n",
            hardware_name
        );
        status = NODE_NOT_FOUND_ERROR;
        goto free;
    }

    set_value = (char*)xmlNodeGetContent(set_value_node);
    is_value = (char*)xmlNodeGetContent(is_value_node);

    if (set_value == NULL || is_value == NULL) {
        fprintf(
            stderr,
            "Hardware readout node '%s' does not contain Set- or Is-value\n",
            hardware_name
        );
        status = EMPTY_FIELD_ERROR;
        goto free;
    }

    char *set_appendix = ":Set Value";
    char *is_appendix = ":Is Value";

    size_t name_size = strlen(hardware_name);
    size_t set_appendix_size = strlen(set_appendix);
    size_t is_appendix_size = strlen(is_appendix);

    set_key = (char*)malloc(name_size + set_appendix_size + 1);
    is_key = (char*)malloc(name_size + is_appendix_size + 1);

    if (set_key == NULL || is_key == NULL) {
        fprintf(stderr, "Could not allocate hardware readout keys\n");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    sprintf(set_key, "%s%s", hardware_name, set_appendix);
    sprintf(is_key, "%s%s", hardware_name, is_appendix);

    set_key[name_size + set_appendix_size] = '\0';
    is_key[name_size + is_appendix_size] = '\0';

    if (!hashmap_insert(metadata, set_key, set_value)) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    if (!hashmap_insert(metadata, is_key, is_value)) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

free:

    if (hardware_id != NULL) free(hardware_id);
    if (hardware_name != NULL) free(hardware_name);
    if (status != SUCCESS) {
        if (set_value != NULL) free(set_value);
        if (is_value != NULL) free(is_value);
        if (set_key != NULL) free(set_key);
        if (is_key != NULL) free(is_key);
    }

    return status;
}


static char *get_or_parse_detname(
    char *det_id,
    Hashmap *detectors,
    Hashmap *parsed_detnames
) {
    char *detname = hashmap_get(parsed_detnames, det_id);

    if (detname != NULL) return detname;

    xmlNodePtr det_node = hashmap_get(detectors, det_id);

    if (det_node == NULL) {
        fprintf(
            stderr,
            "Spectrum is referencing detector information node '%s', "
            "but no such detector is known\n",
            det_id
        );
        return NULL;
    }

    xmlNodePtr name_node = findElementChild(det_node, "RadDetectorName");

    if (name_node == NULL) {
        fprintf(
            stderr, "Detector information does not contain detector name\n"
        );
        return NULL;
    }

    detname = (char*)xmlNodeGetContent(name_node);

    if (detname == NULL) {
        fprintf(
            stderr, "Detector information does not contain detector name\n"
        );
        return NULL;
    }

    if (!hashmap_insert(parsed_detnames, det_id, detname)) {
        free(detname);
        return NULL;
    }

    return detname;
}


static d_ecal_t *get_or_parse_ecal(
    char* ecal_id,
    Hashmap *ecals,
    Hashmap *parsed_ecals
) {
    d_ecal_t *ecal_ptr = hashmap_get(parsed_ecals, ecal_id);
    if (ecal_ptr != NULL) return ecal_ptr;

    d_ecal_t *ecal = NULL;
    char *ecal_values = NULL;

    xmlNodePtr ecal_node = hashmap_get(ecals, ecal_id);

    if (ecal_node == NULL) {
        fprintf(
            stderr,
            "Spectrum is referencing energy calibration node '%s', "
            "but no such calibration is known\n",
            ecal_id
        );
        return NULL;
    }

    xmlNodePtr ecal_value_node = findElementChild(
        ecal_node, "CoefficientValues"
    );

    if (ecal_value_node == NULL) {
        fprintf(
            stderr, "Energy calibration does not contain coefficient values\n"
        );
        return NULL;
    }

    ecal_values = (char*)xmlNodeGetContent(ecal_value_node);

    if (ecal_values == NULL) {
        fprintf(
            stderr, "Energy calibration does not contain coefficient values\n"
        );
        goto free;
    }

    ecal = (d_ecal_t*)calloc(1, sizeof(d_ecal_t));
    if (ecal == NULL) return NULL;

    char *ptr = ecal_values, *endptr;

    ecal->c1 = strtod(ecal_values, &endptr);
    ecal->c2 = strtod(endptr, NULL);

    if (ecal->c2 == 0) {
        fprintf(stderr, "Invalid energy calibration format\n");
        free(ecal);
        ecal = NULL;
        goto free;
    }

    if (!hashmap_insert(parsed_ecals, ecal_id, ecal)) {
        free(ecal);
        ecal = NULL;
    }

free:

    if (ecal_values != NULL) free(ecal_values);
    if (ecal == NULL && ecal_id != NULL) free(ecal_id);

    return ecal;
}


static size_t count_whitespace(char *str) {
    size_t count = 0;
    while (*str != '\0') {
        if (*str == ' ') count++;
        str++;
    }
    return count;
}


static SlopeN42ImportStatus parse_channel_data(
    char *channel_data,
    uint64_t *arr
) {
    uint64_t number;
    char *ptr = channel_data, *endptr;
    size_t i = 0;

    errno = 0;

    while (*ptr != '\0') {
        number = strtoull(ptr, &endptr, 10);
        if (endptr == ptr || errno != 0) break;
        arr[i++] = number;
        ptr = endptr;
    }

    if (i == 0 || errno != 0) {
        fprintf(stderr, "Spectrum has invalid format\n");
        return SPECTRUM_PARSE_ERROR;
    }

    return SUCCESS;
}


static SlopeN42ImportStatus parse_dbspectrum(
    xmlNodePtr spectrum_node,
    Spectrum *spectrum
) {
    SlopeN42ImportStatus status = SUCCESS;

    char *channel_data = NULL;
    uint64_t *arr = NULL;

    xmlNodePtr channel_data_node = findElementChild(
        spectrum_node, "ChannelData"
    );

    if (channel_data_node == NULL) {
        fprintf(stderr, "Spectrum node does not contain ChannelData node\n");
        return NODE_NOT_FOUND_ERROR;
    }

    channel_data = (char*)xmlNodeGetContent(channel_data_node);

    if (channel_data == NULL) {
        fprintf(stderr, "Spectrum node does not contain any data\n");
        return EMPTY_FIELD_ERROR;
    }

    channel_data = trim(channel_data);

    if (strlen(channel_data) == 0) {
        fprintf(stderr, "Spectrum node does not contain any data\n");
        status = EMPTY_FIELD_ERROR;
        goto free;
    }

    size_t size = count_whitespace(channel_data) + 1;

    arr = (uint64_t*)malloc(size * sizeof(uint64_t));

    if (arr == NULL) {
        fprintf(stderr, "Could not allocate array for Spectrum\n");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    status = parse_channel_data(channel_data, arr);

    spectrum->data = arr;
    spectrum->size = size;

    toMinTypeSpectrum(spectrum);

free:

    if (channel_data != NULL) free(channel_data);
    if (status != SUCCESS && arr != NULL) free(arr);

    return status;
}


static SlopeN42ImportStatus import_dbspectrum(
    xmlNodePtr spectrum_node,
    Hashmap *detectors,
    Hashmap *ecals,
    Hashmap *parsed_detnames,
    Hashmap *parsed_ecals,
    Measurement *m
) {
    SlopeN42ImportStatus status = SUCCESS;
    char *det_id = NULL, *detname = NULL, *ecal_id = NULL;
    char *owned_detname_1 = NULL, *owned_detname_2 = NULL;
    d_ecal_t *ecal = NULL;
    DBSpectrum *d = NULL;

    det_id = (char*)xmlGetProp(
        spectrum_node, (xmlChar*)"radDetectorInformationReference"
    );

    if (det_id == NULL) {
        fprintf(
            stderr,
            "Could not find attribute 'radDetectorInformationReference' "
            "in node 'Spectrum'\n"
        );
        return ATTRIBUTE_NOT_FOUND_ERROR;
    }

    detname = get_or_parse_detname(det_id, detectors, parsed_detnames);

    if (detname == NULL) {
        free(det_id);
        status = DANGLING_REFERENCE_ERROR;
        goto free;
    }

    if (hashmap_get(m->dbs, detname) != NULL) {
        fprintf(
            stderr,
            "Detectors of multiple spectra have the same name '%s'\n",
            detname
        );
        status = DUPLICATE_NAME_ERROR;
        goto free;
    }

    owned_detname_1 = strdup(detname);
    owned_detname_2 = strdup(detname);

    if (owned_detname_1 == NULL || owned_detname_2 == NULL) {
        fprintf(stderr, "Could not copy detname\n");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    ecal_id = (char*)xmlGetProp(
        spectrum_node, (xmlChar*)"energyCalibrationReference"
    );

    if (ecal_id == NULL) {
        fprintf(
            stderr,
            "Could not find attribute 'energyCalibrationReference' "
            "in node 'Spectrum'\n"
        );
        status = ATTRIBUTE_NOT_FOUND_ERROR;
        goto free;
    }

    ecal = get_or_parse_ecal(ecal_id, ecals, parsed_ecals);

    if (ecal == NULL) {
        status = INVALID_ECAL_OR_MALLOC_ERROR;
        goto free;
    }

    Spectrum spectrum;

    status = parse_dbspectrum(spectrum_node, &spectrum);

    if (status != SUCCESS) {
        goto free;
    }

    d = newDBSpectrum();

    if (d == NULL) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    d->detname = owned_detname_1;
    d->spectrum = spectrum;
    d->ecal = *ecal;

    if (!hashmap_insert(m->dbs, owned_detname_2, d)) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

free:

    if (status != SUCCESS) {
        if (d == NULL) {
            if (owned_detname_1 != NULL) free(owned_detname_1);
        } else {
            freeDBSpectrum(d);
        }
        if (owned_detname_2 != NULL) free(owned_detname_2);
    }

    return status;
}


static char *parse_detpair(xmlNodePtr detpair_node) {
    xmlNodePtr name_node = findElementChild(detpair_node, "RadDetectorName");

    if (name_node == NULL) {
        fprintf(stderr, "Detector pair node does not contain a pair name\n");
        return NULL;
    }

    char *detpair = (char*)xmlNodeGetContent(name_node);

    if (detpair == NULL) {
        fprintf(stderr, "Detector pair is missing name\n");
        return NULL;
    }

    return detpair;
}


static d_ecal_t get_ecal(
    xmlNodePtr detpair_node,
    Hashmap *detectors,
    Hashmap *parsed_detnames,
    Measurement *m,
    const char *det_node_name
) {
    d_ecal_t ecal = {0};

    xmlNodePtr det_node = findElementChild(detpair_node, det_node_name);

    if (det_node == NULL) {
        fprintf(
            stderr,
            "Could not find node '%s' in the expected place\n",
            det_node_name
        );
        return ecal;
    }

    char *det_id = (char*)xmlGetProp(
        det_node, (xmlChar*)"radDetectorInformationReference"
    );

    if (det_id == NULL) {
        fprintf(
            stderr,
            "Could not find attribute 'radDetectorInformationReference' "
            "in node %s\n",
            det_node_name
        );
        goto free;
    }

    char *detname = get_or_parse_detname(det_id, detectors, parsed_detnames);

    if (detname == NULL) {
        fprintf(
            stderr,
            "Detector pair is referencing detector '%s', but no such "
            "detector is known\n",
            det_id
        );
        goto free;
    }

    DBSpectrum *d = hashmap_get(m->dbs, detname);

    if (d == NULL) {
        fprintf(
            stderr,
            "Detector pair is referencing detector '%s', but no such "
            "detector is known\n",
            detname
        );
        goto free;
    }

    ecal = d->ecal;

free:

    if (det_id != NULL) free(det_id);

    return ecal;
}


static c_ecal_t get_or_parse_c_ecal(
    xmlNodePtr detpair_node,
    Hashmap *detectors,
    Hashmap *parsed_detnames,
    Measurement *m
) {
    return (c_ecal_t){
        get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector1Name"
        ),
        get_ecal(
            detpair_node, detectors, parsed_detnames, m, "RadDetector2Name"
        ),
    };
}


static char *build_png_filepath(
    const char* png_filename,
    const char* path
) {
    #ifdef _WIN32
        const char *dir_end = strrchr(path, '\\');
    #else
        const char *dir_end = strrchr(path, '/');
    #endif

    size_t dir_len = dir_end == NULL ? 0 : dir_end - path + 1;
    size_t file_len = strlen(png_filename);

    char *png_filepath = (char*)calloc(dir_len + file_len + 1, sizeof(char));

    if (png_filepath == NULL) {
        fprintf(stderr, "Could not allocate memory for png file path\n");
        return NULL;
    }

    memcpy(png_filepath, path, dir_len);
    memcpy(png_filepath + dir_len, png_filename, file_len);

    return png_filepath;
}


static Spectrum2D parse_cdbspectrum(
    xmlNodePtr spectrum_node,
    const char* path
) {
    Spectrum2D spectrum = {0};
    char *png_filename = NULL, *png_filepath = NULL;

    png_filename = (char*)xmlNodeGetContent(spectrum_node);

    if (png_filename == NULL) {
        fprintf(stderr, "CDB Spectrum node does not contain data\n");
        goto free;
    }

    png_filename = trim(png_filename);

    if (strlen(png_filename) == 0) {
        fprintf(stderr, "CDB Spectrum node does not contain PNG file name\n");
        goto free;
    }

    png_filepath = build_png_filepath(png_filename, path);

    if (png_filepath == NULL) goto free;

    Image png_content = import_png(png_filepath);

    if (png_content.data == NULL) goto free;

    spectrum.data = png_content.data;
    spectrum.width = png_content.width;
    spectrum.height = png_content.height;

    toMinTypeSpectrum2D(&spectrum);

free:

    if (png_filepath != NULL) free(png_filepath);
    if (png_filename != NULL) xmlFree(png_filename);

    return spectrum;
}


static xmlNodePtr get_detpair_node(
    xmlNodePtr spectrum_node,
    Hashmap *detpairs,
    SlopeN42ImportStatus *status
) {
    char *detpair_id = (char*)xmlGetProp(
        spectrum_node, (xmlChar*)"radDetectorInformationReference"
    );

    if (detpair_id == NULL) {
        fprintf(stderr, "Spectrum node is missing detector reference\n");
        *status = ATTRIBUTE_NOT_FOUND_ERROR;
        return NULL;
    }

    xmlNodePtr detpair_node = hashmap_get(detpairs, detpair_id);

    if (detpair_node == NULL) {
        fprintf(
            stderr,
            "Spectrum node is referencing detector pair with id '%s', but no "
            "such detector pair is known\n",
            detpair_id
        );
        *status = DANGLING_REFERENCE_ERROR;
        goto free;
    }

free:

    free(detpair_id);

    if (*status == SUCCESS) {
        return detpair_node;
    } else {
        return NULL;
    }
}


static SlopeN42ImportStatus import_cdbspectrum(
    xmlNodePtr spectrum_node,
    Hashmap *detpairs,
    Hashmap *detectors,
    Hashmap *parsed_detnames,
    Measurement *m,
    const char* path
) {
    SlopeN42ImportStatus status = SUCCESS;
    char *detpair = NULL;
    char *owned_detpair_1 = NULL, *owned_detpair_2 = NULL;
    Spectrum2D spectrum = {0};
    CDBSpectrum *c = NULL;

    xmlNodePtr detpair_node = get_detpair_node(
        spectrum_node, detpairs, &status
    );

    if (status != SUCCESS) goto free;

    detpair = parse_detpair(detpair_node);

    if (detpair == NULL) {
        status = NODE_NOT_FOUND_ERROR;
        goto free;
    }

    owned_detpair_1 = strdup(detpair);
    owned_detpair_2 = strdup(detpair);

    if (owned_detpair_1 == NULL || owned_detpair_2 == NULL) {
        fprintf(stderr, "Could not copy detname\n");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    if (hashmap_get(m->cdbs, detpair) != NULL) {
        fprintf(
            stderr,
            "Detectors of multiple spectra have the same name '%s'\n",
            detpair
        );
        status = DUPLICATE_NAME_ERROR;
        goto free;
    }

    c_ecal_t ecal = get_or_parse_c_ecal(
        detpair_node, detectors,parsed_detnames, m
    );

    if (ecal.ecal_1.c2 == 0 || ecal.ecal_2.c2 == 0) {
        status = INVALID_ECAL_OR_MALLOC_ERROR;
        goto free;
    }

    spectrum = parse_cdbspectrum(spectrum_node, path);

    if (spectrum.data == NULL) {
        status = PNG_IMPORT_ERROR;
        goto free;
    }

    c = newCDBSpectrum();

    if (c == NULL) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    c->detpair = owned_detpair_1,
    c->spectrum = spectrum;
    c->ecal = ecal;

    if (!hashmap_insert(m->cdbs, owned_detpair_2, c)) {
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

free:

    if (status != SUCCESS) {
        if (c == NULL) {
            if (owned_detpair_1 != NULL) free(owned_detpair_1);
        } else {
            freeCDBSpectrum(c);
        }
        if (owned_detpair_2 != NULL) free(owned_detpair_2);
    }
    if (detpair != NULL) free(detpair);

    return status;
}


static SlopeN42ImportStatus gobble_spectrum(
    xmlNodePtr spectrum_node,
    Hashmap *detectors,
    Hashmap *ecals,
    Hashmap *parsed_detnames,
    Hashmap *parsed_ecals,
    Measurement *m,
    Vector *cdbs_nodes
) {
    SlopeN42ImportStatus status = SUCCESS;

    char *id = (char*)xmlGetProp(spectrum_node, (xmlChar*)"id");

    if (id == NULL) {
        fprintf(stderr, "Spectrum node does not have an id\n");
        return ATTRIBUTE_NOT_FOUND_ERROR;
    }

    if (strstr(id, "Coinc")) {
        vector_push(cdbs_nodes, &spectrum_node);
    } else {
        status = import_dbspectrum(
            spectrum_node, detectors, ecals, parsed_detnames, parsed_ecals, m
        );
    }

    free(id);

    return status;
}


static SlopeN42ImportStatus gobble_meas_child(
    xmlNodePtr node,
    Measurement *m,
    const xmlNodePtr meas_node,
    Hashmap *ecals,
    Hashmap *detectors,
    Hashmap *hardware,
    Hashmap *parsed_detnames,
    Hashmap *parsed_ecals,
    Vector *cdbs_nodes
) {
    SlopeN42ImportStatus status = SUCCESS;

    char *tag = (char*)node->name;

    if (strcmp(tag, "Spectrum") == 0) {
        status = gobble_spectrum(
            node, detectors, ecals, parsed_detnames, parsed_ecals, m, cdbs_nodes
        );
    } else if (strcmp(tag, "Readout") == 0) {
        import_hardware_readout(node, hardware, m->metadata);
    } else if (strstr(tag, "Metadata")) {
        gobble_metadata(node, m->metadata);
    } else {
        status = gobble_generic_text_node(node, m->metadata);
    }

    return status;
}


static SlopeN42ImportStatus sort_meas_children(
    Measurement *m,
    const xmlNodePtr meas_node,
    Hashmap *ecals,
    Hashmap *detectors,
    Hashmap *detpairs,
    Hashmap *hardware,
    const char *path
) {
    SlopeN42ImportStatus status = SUCCESS;

    Hashmap *parsed_ecals = hashmap_new_with_deallocator(0, free);
    Hashmap *parsed_detnames = hashmap_new_with_deallocator(0, free);
    Vector *cdbs_nodes = vector_new(20, sizeof(xmlNodePtr));

    if (parsed_ecals == NULL || parsed_detnames == NULL || cdbs_nodes == NULL) {
        perror("Could not allocate xml node lookup buffers");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    for (
        xmlNodePtr node = xmlFirstElementChild(meas_node);
        node != NULL;
        node = xmlNextElementSibling(node)
    ) {
        status = gobble_meas_child(
            node,
            m,
            meas_node,
            ecals,
            detectors,
            hardware,
            parsed_detnames,
            parsed_ecals,
            cdbs_nodes
        );
        if (status != SUCCESS) goto free;
    }

    for (size_t i = 0; i < cdbs_nodes->size; ++i) {
        xmlNodePtr *node = vector_get(cdbs_nodes, i);
        if (node == NULL) continue;

        status = import_cdbspectrum(
            *node,
            detpairs,
            detectors,
            parsed_detnames,
            m,
            path
        );
    }

free:

    if (parsed_detnames != NULL) hashmap_free(parsed_detnames);
    if (parsed_ecals != NULL) hashmap_free(parsed_ecals);
    if (cdbs_nodes != NULL) vector_free(cdbs_nodes);

    return status;
}


static SlopeN42ImportStatus save_filename(
    Measurement *m,
    const char *filename
) {
    size_t filename_size = strlen(filename);

    if (filename_size > 200) {
        fprintf(stderr, "File path too long\n");
        return INVALID_FILEPATH_ERROR;
    }

    m->filename = (char*)calloc(filename_size + 1, sizeof(char));
    m->name = (char*)calloc(filename_size + 1, sizeof(char));

    if (m->filename == NULL || m->name == NULL) {
        perror("Could not save file path");
        if (m->filename != NULL) free(m->filename);
        if (m->name != NULL) free(m->name);
        return MEMORY_ALLOCATION_ERROR;
    }

    strcpy(m->filename, filename);
    strcpy(m->name, filename);

    return SUCCESS;
}

static SlopeN42ImportStatus import_m(
    Measurement *m,
    xmlNodePtr root,
    const char *filename
) {
    SlopeN42ImportStatus status = save_filename(m, filename);

    if (status != SUCCESS) return status;

    Hashmap *ecals = hashmap_new(0);
    Hashmap *detectors = hashmap_new(0);
    Hashmap *detpairs = hashmap_new(0);
    Hashmap *hardware = hashmap_new(0);
    xmlNodePtr measurements = NULL;

    if (
        ecals == NULL || detectors == NULL ||
        detpairs == NULL || hardware == NULL
    ) {
        perror("Could not allocate xml node lookup buffers");
        status = MEMORY_ALLOCATION_ERROR;
        goto free;
    }

    status = sort_root_children(
        root, &measurements, m->metadata, ecals, detectors, detpairs, hardware
    );

    if (status != SUCCESS) goto free;

    if (measurements == NULL) {
        fprintf(
            stderr,
            "Could not find required node 'RadMeasurement' in the "
            "expected place\n"
        );
        status = NODE_NOT_FOUND_ERROR;
        goto free;
    }

    status = sort_meas_children(
        m, measurements, ecals, detectors, detpairs, hardware, filename
    );

free:

    if (ecals != NULL) hashmap_free(ecals);
    if (detectors != NULL) hashmap_free(detectors);
    if (detpairs != NULL) hashmap_free(detpairs);
    if (hardware != NULL) hashmap_free(hardware);

    return status;
}


SlopeN42ImportStatus import_SLOPE_n42(Measurement *m, const char *filename) {
    if (m == NULL) return NOT_INITIALIZED_ERROR;

    xmlDocPtr doc = xmlParseFile(filename);

    if (doc == NULL) {
        perror("Could not open .n42 file");
        return FILE_SYSTEM_ERROR;
    }

    xmlNodePtr root = xmlDocGetRootElement(doc);

    SlopeN42ImportStatus status = import_m(m, root, filename);

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return status;
}
