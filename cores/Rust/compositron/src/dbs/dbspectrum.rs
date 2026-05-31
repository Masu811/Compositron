use crate::core::utils::Spectrum;

pub struct DBSpectrum {
    pub spectrum: Spectrum,
    pub detname: String,
    pub ecal: (f64, f64),
    pub s: f64,
    pub ds: f64,
    pub w: f64,
    pub dw: f64,
    pub v2p: f64,
    pub dv2p: f64,
    pub counts: u64,
    pub dcounts: f64,
    pub peak_counts: f64,
    pub dpeak_counts: f64,
}

impl DBSpectrum {
    pub fn new(spectrum: Spectrum, detname: String, ecal: (f64, f64)) -> Self {
        DBSpectrum {
            spectrum: spectrum,
            detname: detname,
            ecal: ecal,
            s: f64::NAN,
            ds: f64::NAN,
            w: f64::NAN,
            dw: f64::NAN,
            v2p: f64::NAN,
            dv2p: f64::NAN,
            counts: 0,
            dcounts: f64::NAN,
            peak_counts: f64::NAN,
            dpeak_counts: f64::NAN,
        }
    }
}
