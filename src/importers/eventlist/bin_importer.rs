use std::collections::HashMap;
use std::io::Read;
use std::sync::Arc;

use thiserror::Error;

use crate::core::utils::ByteOrder::BigEndian;
use crate::core::utils::{ByteOrder, Detector};
use crate::eventlist::eventlist::{
    AnalysisError, EventList, EventListColumn, EventListColumnType,
    OwnedEventListColumn
};
use crate::owned_elistcol_match;


#[derive(Debug, Error)]
pub enum BinImportError {
    #[error("I/O Error")]
    IOError {
        #[from]
        inner: std::io::Error,
    },

    #[error("Not enough RAM to load bin file")]
    OOM {
        #[from]
        inner: std::collections::TryReserveError,
    },

    #[error("Could not assemble EventList")]
    ConstructionError {
        #[from]
        inner: AnalysisError,
    },
}


fn parse_rows(
    elist_cols: &mut [OwnedEventListColumn],
    data: &[u8],
    byte_order: ByteOrder
) -> Result<(), BinImportError> {
    let mut file_ptr = 0;

    if data.len() == 0 || elist_cols.len() == 0 {
        return Ok(());
    }

    let mut n_bytes_per_row = 0;

    for col in elist_cols.iter() {
        n_bytes_per_row += match col {
            OwnedEventListColumn::U8(_) => std::mem::size_of::<u8>(),
            OwnedEventListColumn::U16(_) => std::mem::size_of::<u16>(),
            OwnedEventListColumn::U32(_) => std::mem::size_of::<u32>(),
            OwnedEventListColumn::U64(_) => std::mem::size_of::<u64>(),
            OwnedEventListColumn::I8(_) => std::mem::size_of::<i8>(),
            OwnedEventListColumn::I16(_) => std::mem::size_of::<i16>(),
            OwnedEventListColumn::I32(_) => std::mem::size_of::<i32>(),
            OwnedEventListColumn::I64(_) => std::mem::size_of::<i64>(),
            OwnedEventListColumn::F32(_) => std::mem::size_of::<f32>(),
            OwnedEventListColumn::F64(_) => std::mem::size_of::<f64>(),
        };
    }

    let nrows = data.len() / n_bytes_per_row;

    for col in elist_cols.iter_mut() {
        owned_elistcol_match!(col, arr => arr.try_reserve_exact(nrows)?)
    }

    macro_rules! column_match {
        ($column:expr, $data:expr, $file_ptr:expr, $from_bytes:ident;
         $( $ctype:ident => $dtype:ty ),* $(,)?) => {
            match $column {
                $(
                    OwnedEventListColumn::$ctype(arr) => {
                        let Some(buffer) = $data[$file_ptr..]
                            .get(..std::mem::size_of::<$dtype>())
                            .and_then(|s| s.try_into().ok())
                        else { break; };

                        let num = <$dtype>::$from_bytes(buffer);
                        $file_ptr += buffer.len();
                        arr.push(num);
                    }
                )*
            }
        };
    }

    match byte_order {
        ByteOrder::LittleEndian => {
            for _ in 0..nrows {
                for col in elist_cols.iter_mut() {
                    column_match!(
                        col, data, file_ptr, from_le_bytes;
                        U8 => u8,
                        U16 => u16,
                        U32 => u32,
                        U64 => u64,
                        I8 => i8,
                        I16 => i16,
                        I32 => i32,
                        I64 => i64,
                        F32 => f32,
                        F64 => f64
                    );
                }
            }
        },
        BigEndian => {
            for _ in 0..nrows {
                for col in elist_cols.iter_mut() {
                    column_match!(
                        col, data, file_ptr, from_be_bytes;
                        U8 => u8,
                        U16 => u16,
                        U32 => u32,
                        U64 => u64,
                        I8 => i8,
                        I16 => i16,
                        I32 => i32,
                        I64 => i64,
                        F32 => f32,
                        F64 => f64
                    );
                }
            }
        }
    }

    Ok(())
}


pub fn import_bin_eventlist_row_major(
    filepath: &str,
    columns: Vec<(String, EventListColumnType)>,
    detector: Detector,
    byte_order: ByteOrder,
    data_block_offset_bytes: u64,
) -> Result<EventList, BinImportError> {
    let mut file = std::fs::File::open(filepath)?;
    let file_len = file.metadata()?.len();

    let mut data = Vec::new();
    data.try_reserve_exact(file_len as usize)?;

    file.read_to_end(&mut data)?;

    let data = &data[data_block_offset_bytes as usize..];

    let mut elist_cols: Vec<OwnedEventListColumn> = Vec::new();

    for (_, col_type) in &columns {
        elist_cols.push(match col_type {
            EventListColumnType::U8 => OwnedEventListColumn::U8(Vec::new()),
            EventListColumnType::U16 => OwnedEventListColumn::U16(Vec::new()),
            EventListColumnType::U32 => OwnedEventListColumn::U32(Vec::new()),
            EventListColumnType::U64 => OwnedEventListColumn::U64(Vec::new()),
            EventListColumnType::I8 => OwnedEventListColumn::I8(Vec::new()),
            EventListColumnType::I16 => OwnedEventListColumn::I16(Vec::new()),
            EventListColumnType::I32 => OwnedEventListColumn::I32(Vec::new()),
            EventListColumnType::I64 => OwnedEventListColumn::I64(Vec::new()),
            EventListColumnType::F32 => OwnedEventListColumn::F32(Vec::new()),
            EventListColumnType::F64 => OwnedEventListColumn::F64(Vec::new()),
        });
    }

    parse_rows(&mut elist_cols, &data, byte_order)?;

    let mut elist = HashMap::new();

    macro_rules! column_match {
        ($column:expr, $col_name:expr, $elist:expr;
         $( $ctype:ident ),* $(,)?) => {
            match $column {
                $(
                    OwnedEventListColumn::$ctype(arr) => {
                        $elist.insert(
                            $col_name.clone(),
                            EventListColumn::$ctype(Arc::from(arr)),
                        );
                    }
                )*
            }
        };
    }

    for ((col_name, _), col) in columns.iter().zip(elist_cols) {
        column_match!(
            col, col_name, elist;
            U8, U16, U32, U64, I8, I16, I32, I64, F32, F64
        )
    }

    Ok(EventList::try_new(detector, elist)?)
}
