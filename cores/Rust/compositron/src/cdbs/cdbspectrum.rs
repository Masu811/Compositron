use crate::core::utils::Spectrum2D;

pub struct CDBSpectrum {
    pub spectrum: Spectrum2D,
    pub detpair: String,
    pub ecal: ((f64, f64), (f64, f64)),
    pub s: f64,
    pub ds: f64,
    pub w: f64,
    pub dw: f64,
    pub counts: u64,
    pub dcounts: f64,
    pub roi_counts: f64,
    pub droi_counts: f64,
}

impl CDBSpectrum {
    pub fn new(
        spectrum: Spectrum2D, detpair: String, ecal: ((f64, f64), (f64, f64))
    ) -> Self {
        CDBSpectrum {
            spectrum: spectrum,
            detpair: detpair,
            ecal: ecal,
            s: f64::NAN,
            ds: f64::NAN,
            w: f64::NAN,
            dw: f64::NAN,
            counts: 0,
            dcounts: f64::NAN,
            roi_counts: f64::NAN,
            droi_counts: f64::NAN,
        }
    }
}
