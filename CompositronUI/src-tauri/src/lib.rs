use std::sync::Mutex;

use thiserror;
use compositron::importers::waveforms::lecroy_trc_importer::{import_lecroy_trc, TRCImportError};
use compositron::waveforms::waveform::{PulsePolarity, Trace, WaveformArray};

#[derive(Debug, thiserror::Error)]
enum Error {
  #[error(transparent)]
  ImportError(#[from] TRCImportError),
}


impl serde::Serialize for Error {
  fn serialize<S: serde::ser::Serializer>(
      &self, serializer: S
  ) -> Result<S::Ok, S::Error> {
    serializer.serialize_str(self.to_string().as_ref())
  }
}


struct WaveformManager {
    pub filepath: Option<String>,
    pub data: Option<WaveformArray>,
}


impl WaveformManager {
    pub fn default() -> Self {
        WaveformManager {
            filepath: None,
            data: None,
        }
    }
}


#[tauri::command]
fn import_waveforms(
    filepath: &str,
    state: tauri::State<'_, Mutex<WaveformManager>>,
) -> Result<usize, Error> {
    let file = import_lecroy_trc(
        filepath,
        "A".into(),
        PulsePolarity::Negative,
    )?;

    let n_traces = file.traces.len();

    let mut state = state.lock().unwrap();

    state.filepath = Some(filepath.into());
    state.data = Some(file);

    Ok(n_traces)
}


#[tauri::command]
fn get_trace_data(
    index: usize,
    state: tauri::State<'_, Mutex<WaveformManager>>,
) -> Result<Vec<f64>, Error> {
    let state = state.lock().unwrap();

    let trace = state.data.as_ref().unwrap().traces.get(index).unwrap();

    match trace {
        Trace::I8(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
        Trace::I16(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
        Trace::I32(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
        Trace::I64(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
        Trace::F32(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
        Trace::F64(arr) => {
            Ok(arr.iter().map(|&x| x as f64).collect::<Vec<_>>())
        },
    }
}


#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(Mutex::new(WaveformManager::default()))
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            import_waveforms,
            get_trace_data,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
