use std::fs::File;
use std::io::{Seek, Read};
use std::path::Path;

use thiserror::Error;

use crate::core::utils::LinearCalibration;
use crate::waveforms::waveform::{PeakPolarity, Trace, WaveformArray};


#[derive(Debug, Error)]
pub enum TRCImportError {
    #[error("I/O Error")]
    IOError {
        #[from]
        inner: std::io::Error,
    },

    #[error("Invalid UTF-8")]
    UTF8Error {
        #[from]
        inner: std::string::FromUtf8Error,
    },
}


#[derive(Debug, Clone, Copy)]
pub enum CommType {
    Byte,
    Word,
}


#[derive(Debug)]
pub enum DataArray {
    Byte(Vec<i8>),
    Word(Vec<i16>),
}


#[derive(Debug)]
pub struct TrcFile {
    pub header_size: usize,
    pub byte_order: ByteOrder,
    pub descriptor: WaveDesc1,
    pub usertext: String,
    pub trigtime: Vec<(f64, f64)>,
    pub ristime: Vec<f64>,
    pub data_array_1: DataArray,
    pub data_array_2: DataArray,
}


#[derive(Debug)]
pub struct WaveDesc1 {
    pub descriptor_name: String,  // Will contain "WAVEDESC"
    pub template_name: String,  // Name of the template
    pub comm_type: i16,  // Data format: 0 = byte; 1 = word
    pub comm_order: i16,  // Data order: 0 Hi First; 1 = Lo First
    // The following variables specify the lengths of all blocks of which the entire waveform (as it is currently being read) is composed. If a block length is zero, that block is (currently) not present. Blocks and arrays that are present will be found in the same order as shown below.
    pub wave_desc_length: i32,  // Length in bytes of WAVEDESC1 block
    pub user_text_length: i32,  // Length in bytes of USERTEXT block
    pub res_desc1: i32,  // Reserved
    pub trig_time_array: i32,  // Length in bytes of TRIGTIME array
    pub ris_time_array: i32,  // Length in bytes of RISTIME array
    pub res_array_1: i32,  // Reserved
    pub wave_array_1: i32,  // Length in bytes of first simple data array
    pub wave_array_2: i32,  // Length in bytes of second simple data array
    pub res_array_2: i32,  // Reserved
    pub res_array_3: i32,  // Reserved
    // The following variables identify the instrument.
    pub instrument_name: String,
    pub instrument_number: i32,
    pub trace_label: String,  // Waveform identifier
    pub reserved_data_count: i32,  // Reserved
    // The following variables describe the waveform type and the time at which the waveform was generated.
    pub wave_array_count: i32,  // Number of data points in a data array. If there are two arrays (e.g., FFT or Extrema waveform), this number applies to each array separately.
    pub points_per_screen: i32,  // Nominal number of data points on the screen
    pub first_valid: i32,  // Number of points to skip before forst good point. FIRST_VALID_POINT = 0 for normal waveforms.
    pub last_valid: i32,  // Index of last good data point in record before padding (blanking) was started. LAST_VALID_POINT = WAVE_ARRAY_COUNT - 1 except for aborted Sequence and Roll Mode acquisitions.
    pub first_point: i32,  // Indicates the data offset relative to the beginning of the trace buffer
    pub sparsing_factor: i32,  // Indicates the sparsing into data block
    pub segment_no: i32,  // For Sequence waveforms, index of the segment
    pub subarray_count: i32,  // For Sequence waveforms, acquired segment count, between 0 and NOM_SUBARRAY_COUNT
    pub sweeps_per_acq: i32,  // For Average or Extrema waveforms, number of sweeps accumulated, else 1
    pub points_per_pair: i16,  // For Peak Dectect waveforms (which always include data points in DATA_ARRAY_1 and min/max pairs in DATA_ARRAY_2), the number of data points for each min/max pair
    pub pair_offset: i16,  // For Peak Dectect waveforms, the number of data points by which the first min/max pair in DATA_ARRAY_2 is offset relative to the first data value in DATA_ARRAY_1
    pub vertical_gain: f32,  // Total gain of waveform, units per lsb
    pub vertical_offset: f32,  // Total vertical offset of waveform. To get floating values from raw data: VERTICAL_GAIN * data - VERTICAL_OFFSET
    pub max_value: f32,  // Maximum allowed value; corresponds to the upper edge of the grid
    pub min_value: f32,  // Minimum allowed value; corresponds to the lower edge of the grid
    pub nominal_bits: i16,  // Intrinsic precision of the observation
    pub nom_subarray_counts: i16,  // For Sequence waveforms, nominal segment count, else 1
    pub horizontal_interval: f32,  // Sampling interval, the nominal time between successive points in the data
    pub horizontal_offset: f64,  // Trigger offset in time domain for zero'th sweep of trigger, measured as seconds from trigger to zero'th data point (i.e., actual trigger delay)
    pub pixel_offset: f64,  // Time from trigger to zero'th pixel of display segment in time domain, measured in seconds (i.e., nominal trigger delay)
    pub vert_unit: String,  // Vertical axis unit
    pub hor_unit: String,  // Horizontal axis unit
    pub horiz_uncertainty: f32,  // Uncertainty from one acquisition to the next, of the horizontal offset in seconds
    pub trigger_time_second: f64,
    pub trigger_time_minute: u8,
    pub trigger_time_hour: u8,
    pub trigger_time_day: u8,
    pub trigger_time_month: u8,
    pub trigger_time_year: i16,
    pub trigger_time_dummy: i16,
    pub acq_duration: f32,  // Duration of the acquisition (in sec) for multi-trigger waveforms (e.g., Sequence, RIS and Average)
    pub ca_record_type: i16,  // Type of waveform
    pub processing_done: i16,  // Indication of any processing done. 0 = no processing; 1 = fir filter; 2 = interpolated; 3 = sparsed; 4 = autoscaled; 5 = no result; 6 = rolling; 7 = cumulative.
    pub reserved5: i16,
    pub ris_sweeps: i16,  // For RIS acquisitions, number of sweeps from which waveform is calculated, else 1
    // **The information below is based on the legacy enumeration list. It may not be valid for newer oscilloscope models.
    pub time_base: i16,  // Enumerated time/div **
    pub vertical_coupling: i16,  // Enumerated channel coupling value
    pub probe_attenuation: f32,  // Probe attenuation value
    pub fixed_vertical_gain: i16,  // Enumerated vertical gain **
    pub band_width_limit: i16,
    pub vertical_vernier: f32,
    pub acq_vertical_offset: f32,
    pub wave_source: i16,  // Waveform source input
}


#[repr(usize)]
pub enum WaveDesc1FieldOffset {
    DescriptorName = 0,
    TemplateName = 16,
    CommType = 32,
    CommOrder = 34,
    WaveDescLength = 36,
    UserTextLength = 40,
    ResDesc1 = 44,
    TrigTimeArray = 48,
    RisTimeArray = 52,
    ResArray1 = 56,
    WaveArray1 = 60,
    WaveArray2 = 64,
    ResArray2 = 68,
    ResArray3 = 72,
    InstrumentName = 76,
    InstrumentNumber = 92,
    TraceLabel = 96,
    ReservedDataCount = 112,
    WaveArrayCount = 116,
    PointsPerScreen = 120,
    FirstValid = 124,
    LastValid = 128,
    FirstPoint = 132,
    SparsingFactor = 136,
    SegmentNo = 140,
    SubarrayCount = 144,
    SweepsPerAcq = 148,
    PointsPerPair = 152,
    PairOffset = 154,
    VerticalGain = 156,
    VerticalOffset = 160,
    MaxValue = 164,
    MinValue = 168,
    NominalBits = 172,
    NomSubarrayCounts = 174,
    HorizontalInterval = 176,
    HorizontalOffset = 180,
    PixelOffset = 188,
    VertUnit = 196,
    HorUnit = 244,
    HorizUncertainty = 292,
    TriggerTimeSeconds = 296,
    TriggerTimeMinutes = 304,
    TriggerTimeHours = 305,
    TriggerTimeDays = 306,
    TriggerTimeMonths = 307,
    TriggerTimeYear = 308,
    TriggerTimeDummy = 310,
    AcqDuration = 312,
    CaRecordType = 316,
    ProcessingDone = 318,
    Reserved5 = 320,
    RisSweeps = 322,
    TimeBase = 324,
    VerticalCoupling = 326,
    ProbeAttenuation = 328,
    FixedVerticalGain = 332,
    BandWidthLimit = 334,
    VerticalVernier = 336,
    AcqVerticalOffset = 340,
    WaveSource = 344,
}


#[derive(Debug, Clone, Copy)]
pub enum ByteOrder {
    LittleEndian,
    BigEndian,
}


trait FromBytes<const N: usize> {
    fn from_bytes(bytes: [u8; N], order: ByteOrder) -> Self;
}


macro_rules! impl_from_bytes {
    ($target_type:ident, $size:expr) => {
        impl FromBytes<$size> for $target_type {
            fn from_bytes(bytes: [u8; $size], order: ByteOrder) -> Self {
                match order {
                    ByteOrder::LittleEndian => $target_type::from_le_bytes(bytes),
                    ByteOrder::BigEndian => $target_type::from_be_bytes(bytes),
                }
            }
        }
    };
}


impl_from_bytes!(i16, 2);
impl_from_bytes!(i32, 4);
impl_from_bytes!(f32, 4);
impl_from_bytes!(f64, 8);
impl_from_bytes!(u8, 1);


fn read_num<const N: usize, T: FromBytes<N>>(
    file: &mut File,
    order: ByteOrder,
) -> Result<T, TRCImportError> {
    let mut buffer = [0u8; N];
    file.read_exact(&mut buffer)?;
    Ok(T::from_bytes(buffer, order))
}


fn read_str<const N: usize>(file: &mut File) -> Result<String, TRCImportError> {
    let mut buffer = [0u8; N];
    file.read_exact(&mut buffer)?;

    Ok(
        buffer
            .iter()
            .map(|&c| c as char)
            .filter(|&c| c.is_alphanumeric())
            .collect::<String>()
    )
}


fn get_header_size(file: &mut File) -> Result<usize, TRCImportError> {
    file.seek(std::io::SeekFrom::Start(1))?;

    let mut buffer = [0u8; 1];

    file.read_exact(&mut buffer)?;

    Ok((buffer[0] as char).to_digit(16).unwrap() as usize + 2)
}


fn get_comm_type(
    file: &mut File, offset: usize
) -> Result<CommType, TRCImportError> {
    file.seek(
        std::io::SeekFrom::Start(
            (offset + WaveDesc1FieldOffset::CommType as usize) as u64
        )
    )?;

    let mut buffer = [0u8; 2];
    file.read_exact(&mut buffer)?;

    if i16::from_le_bytes(buffer) == 0 {
        Ok(CommType::Byte)
    } else {
        Ok(CommType::Word)
    }
}


fn get_byte_order(
    file: &mut File, offset: usize
) -> Result<ByteOrder, TRCImportError> {
    file.seek(
        std::io::SeekFrom::Start(
            (offset + WaveDesc1FieldOffset::CommOrder as usize) as u64
        )
    )?;

    let mut buffer = [0u8; 2];
    file.read_exact(&mut buffer)?;

    if i16::from_le_bytes(buffer) == 0 {
        Ok(ByteOrder::BigEndian)
    } else {
        Ok(ByteOrder::LittleEndian)
    }
}


fn read_descriptor_block(
    file: &mut File,
    offset: usize,
    order: ByteOrder,
) -> Result<WaveDesc1, TRCImportError> {
    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    Ok(WaveDesc1 {
        descriptor_name: read_str::<16>(file)?,
        template_name: read_str::<16>(file)?,
        comm_type: read_num(file, order)?,
        comm_order: read_num(file, order)?,
        wave_desc_length: read_num(file, order)?,
        user_text_length: read_num(file, order)?,
        res_desc1: read_num(file, order)?,
        trig_time_array: read_num(file, order)?,
        ris_time_array: read_num(file, order)?,
        res_array_1: read_num(file, order)?,
        wave_array_1: read_num(file, order)?,
        wave_array_2: read_num(file, order)?,
        res_array_2: read_num(file, order)?,
        res_array_3: read_num(file, order)?,
        instrument_name: read_str::<16>(file)?,
        instrument_number: read_num(file, order)?,
        trace_label: read_str::<16>(file)?,
        reserved_data_count: read_num(file, order)?,
        wave_array_count: read_num(file, order)?,
        points_per_screen: read_num(file, order)?,
        first_valid: read_num(file, order)?,
        last_valid: read_num(file, order)?,
        first_point: read_num(file, order)?,
        sparsing_factor: read_num(file, order)?,
        segment_no: read_num(file, order)?,
        subarray_count: read_num(file, order)?,
        sweeps_per_acq: read_num(file, order)?,
        points_per_pair: read_num(file, order)?,
        pair_offset: read_num(file, order)?,
        vertical_gain: read_num(file, order)?,
        vertical_offset: read_num(file, order)?,
        max_value: read_num(file, order)?,
        min_value: read_num(file, order)?,
        nominal_bits: read_num(file, order)?,
        nom_subarray_counts: read_num(file, order)?,
        horizontal_interval: read_num(file, order)?,
        horizontal_offset: read_num(file, order)?,
        pixel_offset: read_num(file, order)?,
        vert_unit: read_str::<48>(file)?,
        hor_unit: read_str::<48>(file)?,
        horiz_uncertainty: read_num(file, order)?,
        trigger_time_second: read_num(file, order)?,
        trigger_time_minute: read_num(file, order)?,
        trigger_time_hour: read_num(file, order)?,
        trigger_time_day: read_num(file, order)?,
        trigger_time_month: read_num(file, order)?,
        trigger_time_year: read_num(file, order)?,
        trigger_time_dummy: read_num(file, order)?,
        acq_duration: read_num(file, order)?,
        ca_record_type: read_num(file, order)?,
        processing_done: read_num(file, order)?,
        reserved5: read_num(file, order)?,
        ris_sweeps: read_num(file, order)?,
        time_base: read_num(file, order)?,
        vertical_coupling: read_num(file, order)?,
        probe_attenuation: read_num(file, order)?,
        fixed_vertical_gain: read_num(file, order)?,
        band_width_limit: read_num(file, order)?,
        vertical_vernier: read_num(file, order)?,
        acq_vertical_offset: read_num(file, order)?,
        wave_source: read_num(file, order)?,
    })
}


fn read_usertext(
    file: &mut File,
    offset: usize,
    size: usize,
) -> Result<String, TRCImportError> {
    if size == 0 {
        return Ok(String::new());
    }

    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    let mut buffer = vec![0u8; size];
    file.read_exact(&mut buffer)?;

    let s = String::from_utf8(buffer)?;

    Ok(s)
}


fn read_trigtime(
    file: &mut File,
    offset: usize,
    size: usize,
    order: ByteOrder,
) -> Result<Vec<(f64, f64)>, TRCImportError> {
    let n_elements = size / 16;
    let mut result: Vec<(f64, f64)> = Vec::with_capacity(n_elements);

    if size == 0 {
        return Ok(result);
    }

    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    let mut buffer = [0u8; 8];

    match order {
        ByteOrder::LittleEndian => {
            for _ in 0..n_elements {
                file.read_exact(&mut buffer)?;
                let x = f64::from_le_bytes(buffer);
                file.read_exact(&mut buffer)?;
                let y = f64::from_le_bytes(buffer);
                result.push((x, y));
            }
        },
        ByteOrder::BigEndian => {
            for _ in 0..n_elements {
                file.read_exact(&mut buffer)?;
                let x = f64::from_be_bytes(buffer);
                file.read_exact(&mut buffer)?;
                let y = f64::from_be_bytes(buffer);
                result.push((x, y));
            }
        },
    }

    Ok(result)
}


fn read_ristime(
    file: &mut File,
    offset: usize,
    size: usize,
    order: ByteOrder
) -> Result<Vec<f64>, TRCImportError> {
    let n_elements = size / 8;

    let mut result: Vec<f64> = Vec::with_capacity(n_elements);

    if size == 0 {
        return Ok(result);
    }

    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    let mut buffer = [0u8; 8];

    match order {
        ByteOrder::LittleEndian => {
            for _ in 0..n_elements {
                file.read_exact(&mut buffer)?;
                result.push(f64::from_le_bytes(buffer));
            }
        },
        ByteOrder::BigEndian => {
            for _ in 0..n_elements {
                file.read_exact(&mut buffer)?;
                result.push(f64::from_be_bytes(buffer));
            }
        },
    }

    Ok(result)
}


fn read_data_array(
    file: &mut File,
    offset: usize,
    size: usize,
    comm_type: CommType,
    order: ByteOrder,
) -> Result<DataArray, TRCImportError> {
    file.seek(std::io::SeekFrom::Start(offset as u64))?;
    match comm_type {
        CommType::Byte => {
            let mut buffer = vec![0u8; size];
            file.read_exact(&mut buffer)?;

            let result: Vec<i8> = buffer
                    .into_iter()
                    .map(|b| b as i8)
                    .collect();

            Ok(DataArray::Byte(result))
        },
        CommType::Word => {
            let mut buffer = vec![0u8; size];
            file.read_exact(&mut buffer)?;

            let mut result = Vec::with_capacity(size / 2);

            match order {
                ByteOrder::LittleEndian => {
                    for chunk in buffer.chunks_exact(2) {
                        let val = i16::from_le_bytes([chunk[0], chunk[1]]);
                        result.push(val);
                    }
                }
                ByteOrder::BigEndian => {
                    for chunk in buffer.chunks_exact(2) {
                        let val = i16::from_be_bytes([chunk[0], chunk[1]]);
                        result.push(val);
                    }
                }
            }

            Ok(DataArray::Word(result))
        },
    }
}


pub fn import_lecroy_trc(
    path: &Path,
    peak_polarity: PeakPolarity,
) -> Result<WaveformArray, TRCImportError> {
    let mut file = File::open(path)?;

    let mut offset = 0;

    let header_size = get_header_size(&mut file)?;
    offset += header_size;

    let comm_type = get_comm_type(&mut file, offset)?;
    let byte_order = get_byte_order(&mut file, offset)?;
    let descriptor = read_descriptor_block(&mut file, offset, byte_order)?;
    offset += descriptor.wave_desc_length as usize;

    let usertext = read_usertext(
        &mut file,
        offset,
        descriptor.user_text_length as usize,
    )?;
    offset += descriptor.user_text_length as usize;

    let trigtime = read_trigtime(
        &mut file,
        offset,
        descriptor.trig_time_array as usize,
        byte_order
    )?;
    offset += descriptor.trig_time_array as usize;

    let ristime = read_ristime(
        &mut file,
        offset,
        descriptor.ris_time_array as usize,
        byte_order
    )?;
    offset += descriptor.ris_time_array as usize;

    let data_array_1 = read_data_array(
        &mut file,
        offset,
        descriptor.wave_array_1 as usize,
        comm_type,
        byte_order,
    )?;
    offset += descriptor.wave_array_1 as usize;

    let data_array_2 = read_data_array(
        &mut file,
        offset,
        descriptor.wave_array_2 as usize,
        comm_type,
        byte_order,
    )?;

    let trc_file = TrcFile {
        header_size,
        byte_order,
        descriptor,
        usertext,
        trigtime,
        ristime,
        data_array_1,
        data_array_2,
    };

    let tot_size = trc_file.descriptor.wave_array_1;
    let n_subarrays = trc_file.descriptor.subarray_count;
    let dsize = if trc_file.descriptor.comm_type == 0 { 1 } else { 2 };

    let subarray_size = tot_size / (n_subarrays * dsize);

    let mut traces: Vec<Trace> = Vec::with_capacity(n_subarrays as usize);

    match &trc_file.data_array_1 {
        DataArray::Byte(data) => {
            let size = data.len();
            for i in 0..n_subarrays {
                let start: usize = (i * subarray_size) as usize;
                let stop: usize = ((i + 1) * subarray_size) as usize;
                if stop > size {
                    break;
                }
                traces.push(Trace::I8(data[start..stop].to_vec()));
            }
        }
        DataArray::Word(data) => {
            let size = data.len();
            for i in 0..n_subarrays {
                let start: usize = (i * subarray_size) as usize;
                let stop: usize = ((i + 1) * subarray_size) as usize;
                if stop > size {
                    break;
                }
                traces.push(Trace::I16(data[start..stop].to_vec()));
            }
        }
    }

    let mut trigger_times: Vec<f64> = Vec::new();

    for i in 0..n_subarrays as usize {
        trigger_times.push(trc_file.trigtime[i].0 + trc_file.trigtime[i].1);
    }

    Ok(WaveformArray::new(
        traces,
        trc_file.descriptor.horizontal_interval as f64,
        trigger_times,
        LinearCalibration {
            scale: trc_file.descriptor.vertical_gain as f64,
            offset: trc_file.descriptor.vertical_offset as f64,
        },
        peak_polarity,
    ))
}
