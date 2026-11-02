use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::sync::Arc;

use thiserror::Error;

use crate::core::utils::{Detector, read_str_dyn};
use crate::eventlist::eventlist::{CorrelatedEventLists, EventList, EventListColumn};


#[derive(Debug, Error)]
pub enum BinImportError {
    #[error("I/O Error")]
    IOError {
        #[from]
        inner: std::io::Error,
    },

    #[error("Malformed Header")]
    MalformedHeader,

    #[error("Bin file has wrong signature (expected \"CDBSDATA\", got \"{got}\")")]
    WrongSignature {
        got: String,
    },
}


#[derive(Debug, Clone, Copy)]
pub struct Flags(u16);


impl Flags {
    pub fn sorted(&self) -> bool {
        self.0 & 0b001 == 1
    }

    pub fn overflow_fixed(&self) -> bool {
        self.0 & 0b010 == 1
    }

    pub fn excerpt(&self) -> bool {
        self.0 & 0b100 == 1
    }
}


#[derive(Clone, Copy)]
pub struct Signature([u8; 8]);


impl Signature {
    pub fn to_string(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(&self.0)
    }
}


impl Debug for Signature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}


#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Header {
    pub signature: Signature,
    pub flags: Flags,
    pub version: u16,
    pub compression: u32,
    pub timestamp: u64,
    pub board_count: u32,
    pub event_count: u64,
}


const HEADER_SIZE: usize = std::mem::size_of::<Header>();


impl Header {
    pub fn from_file(file: &mut impl Read) -> Result<Self, std::io::Error> {
        let mut buf = [0u8; HEADER_SIZE];

        file.read_exact(&mut buf)?;

        Ok(unsafe {
            std::ptr::read_unaligned(buf.as_ptr().cast::<Header>())
        })
    }
}


#[repr(C, packed)]
#[derive(Debug, Clone, Copy)]
pub struct VersionableStringInfo {
    pub timestamp: u64,  // ns since Unix epoch
    pub next_link: u64,
    pub last_link: u64,
    pub len: u32,
}


const VSTRINGINFO_SIZE: usize = std::mem::size_of::<VersionableStringInfo>();


impl VersionableStringInfo {
    pub fn from_file(file: &mut impl Read) -> Result<Self, std::io::Error> {
        let mut buf = [0u8; VSTRINGINFO_SIZE];

        file.read_exact(&mut buf)?;

        Ok(unsafe {
            std::ptr::read_unaligned(buf.as_ptr().cast::<VersionableStringInfo>())
        })
    }
}


#[derive(Debug, Clone)]
pub struct VersionableString {
    pub info: VersionableStringInfo,
    pub data: String,
}


impl VersionableString {
    pub fn from_file(
        file: &mut impl Read,
    ) -> Result<Self, std::io::Error> {
        let info = VersionableStringInfo::from_file(file)?;
        let data = read_str_dyn(file, info.len as usize)?;

        Ok(VersionableString { info, data })
    }
}


#[derive(Debug, Clone)]
pub struct CDBSBinFileInfo {
    pub header: Header,
    pub settings: Vec<VersionableString>,
    pub metadata: VersionableString,
    pub config: VersionableString,
    pub team: VersionableString,
}


impl CDBSBinFileInfo {
    fn from_file(
        file: &mut impl Read,
    ) -> Result<Self, BinImportError> {
        let header = Header::from_file(file)?;

        if header.signature.to_string() != "CDBSDATA" {
            return Err(BinImportError::WrongSignature {
                got: header.signature.to_string().into()
            });
        }

        let mut settings = Vec::new();

        for _ in 0..header.board_count {
            settings.push(VersionableString::from_file(file)?);
        }

        let metadata = VersionableString::from_file(file)?;
        let config = VersionableString::from_file(file)?;
        let team = VersionableString::from_file(file)?;

        Ok(CDBSBinFileInfo { header, settings, metadata, config, team })
    }
}


#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Entry {
    pub timestamp: u64,
    pub channel: u16,
    pub det_id: u8,
}


const ENTRY_SIZE: usize = std::mem::size_of::<Entry>();


impl Entry {
    pub fn from_file(file: &mut impl Read) -> Result<Self, std::io::Error> {
        let mut buf = [0u8; ENTRY_SIZE];

        file.read_exact(&mut buf)?;

        Ok(unsafe {
            std::ptr::read_unaligned(buf.as_ptr().cast::<Entry>())
        })
    }
}


#[derive(Debug, Clone)]
pub struct CDBSBinFile {
    pub info: CDBSBinFileInfo,
    pub eventlists: CorrelatedEventLists,
}


impl CDBSBinFile {
    pub fn from_file(filepath: &str) -> Result<Self, BinImportError> {
        let file = File::open(filepath)?;
        let mut reader = BufReader::new(file);

        let info = CDBSBinFileInfo::from_file(&mut reader)?;

        let event_count = info.header.event_count as usize;

        let mut det_ids = HashMap::new();

        for _ in 0..event_count {
            let event = Entry::from_file(&mut reader)?;

            if !det_ids.contains_key(&event.det_id) {
                det_ids.insert(event.det_id, (Vec::new(), Vec::new()));
            }

            let det = &mut det_ids.get_mut(&event.det_id).unwrap();

            det.0.push(event.timestamp);
            det.1.push(event.channel);
        }

        let mut elists = HashMap::new();

        for (det_id, elist_cols) in det_ids.into_iter() {
            let detname = format!("Det {det_id}");
            elists.insert(
                detname.clone(),
                EventList::try_new(
                    Detector {
                        name: detname,
                        ecal: None,
                        tcal: None,
                        eres: None,
                        tres: None,
                    },
                    HashMap::from([
                        (
                            "timestamp".into(),
                            EventListColumn::U64(Arc::from(elist_cols.0))
                        ),
                        (
                            "energy".into(),
                            EventListColumn::U16(Arc::from(elist_cols.1))
                        ),
                    ])
                ).unwrap()
            );
        }

        Ok(CDBSBinFile { info, eventlists: CorrelatedEventLists::new(elists) })
    }
}


impl Display for CDBSBinFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CDBSBinFile <{:?}, data: Vec of size {:?}>",
            self.info,
            self.eventlists,
        )
    }
}


#[derive(Debug)]
pub struct LazyCDBSBinFile {
    pub info: CDBSBinFileInfo,
    reader: BufReader<File>,
    data_block_offset: u64,
}


impl LazyCDBSBinFile {
    pub fn from_file(filepath: &str) -> Result<Self, BinImportError> {
        let file = File::open(filepath)?;
        let mut reader = BufReader::new(file);

        let info = CDBSBinFileInfo::from_file(&mut reader)?;

        let data_block_offset = reader.stream_position()?;

        Ok(LazyCDBSBinFile { info, reader, data_block_offset })
    }

    /// Fill buf, get how many entries were read
    pub fn load_buf(
        &mut self, buf: &mut [Entry], n_entries_offset: u64
    ) -> Result<usize, std::io::Error> {
        self.reader.seek(
            std::io::SeekFrom::Start(
                self.data_block_offset + n_entries_offset * ENTRY_SIZE as u64
            )
        )?;

        let count = buf.len();

        let bytes = unsafe {
            std::slice::from_raw_parts_mut(
                buf.as_mut_ptr().cast::<u8>(), count * ENTRY_SIZE,
            )
        };

        Ok(self.reader.read(bytes)? / ENTRY_SIZE)
    }
}


impl Display for LazyCDBSBinFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "LazyCDBSBinFile <{:?}, data_block_offset: {}>",
            self.info,
            self.data_block_offset,
        )
    }
}
