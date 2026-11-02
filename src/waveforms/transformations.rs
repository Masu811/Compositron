use crate::waveforms::waveform::Waveform;


pub trait TransformWaveform {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}


pub struct InvertY;


impl TransformWaveform for InvertY {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}


pub struct ZeroBaseline;


impl TransformWaveform for ZeroBaseline {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}


pub struct SlidingWindowAverage {
    pub window_size: usize,
}


impl TransformWaveform for SlidingWindowAverage {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}


pub struct DownsampleStride {
    pub stride: usize,
}


impl TransformWaveform for DownsampleStride {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}


pub struct DownsampleAverageRightNeighbors {
    pub n: usize,
}


impl TransformWaveform for DownsampleAverageRightNeighbors {
    fn transform(&self, waveform: &Waveform) -> Waveform {
        todo!();
    }
}
