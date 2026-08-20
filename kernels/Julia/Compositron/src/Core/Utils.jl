using ..Constants: KEV_PER_M0C

# LinearCalibration

mutable struct LinearCalibration
    offset::Float64
    scale::Float64
end

function from_index(ecal::LinearCalibration, index)::Float64
    ecal.scale * index + ecal.offset
end

function to_index_rounded(ecal::LinearCalibration, value)::Int
    Int(round((value - ecal.offset) / ecal.scale))
end

function to_index_trunc(ecal::LinearCalibration, value)::Int
    Int(trunc((value - ecal.offset) / ecal.scale))
end

function to_index_float(ecal::LinearCalibration, value)::Float64
    (value - ecal.offset) / ecal.scale
end

# EcalCorrectionOrder

@enum EcalCorrectionOrder zeroth_order first_order no_corr

# Unit

abstract type Unit end

struct KeV <: Unit
    value::Float64
end

struct M0C <: Unit
    value::Float64
end

struct Eres <: Unit
    value::Float64
end

to_kev(value::KeV, eres::Union{Nothing, Float64})::Float64 = value.value

to_kev(value::M0C, eres::Union{Nothing, Float64})::Float64 = value.value * KEV_PER_M0C

to_kev(value::Eres, eres::Union{Nothing, Float64})::Union{Nothing, Float64} = (
    isnothing(eres) ? nothing : value.value * eres
)

# Detectors

# DBS
mutable struct EnergyDetector
    name::String
    ecal::LinearCalibration
    corrected_ecal::Union{Nothing, LinearCalibration}
    eres::Union{Nothing, Float64}
end

# CDBS
mutable struct EnergyDetectorPair
    name::String
    first_det::EnergyDetector
    second_det::EnergyDetector
    eres::Union{Nothing, Float64}
end

# PALS
mutable struct TimingDetectorPair
    name::String
    tcal::LinearCalibration
    corrected_tcal::Union{Nothing, LinearCalibration}
    tres::Union{Nothing, Float64}
end

# Min-Type Spectra

function min_dtype(arr::AbstractArray{<:Integer})::Type{<:Unsigned}
    max_val = maximum(arr)
    if max_val <= typemax(UInt8)
        return UInt8
    elseif max_val <= typemax(UInt16)
        return UInt16
    elseif max_val <= typemax(UInt32)
        return UInt32
    else
        return UInt64
    end
end

function min_spectrum(arr::AbstractVector{<:Integer})::Vector{<:Unsigned}
    dtype = min_dtype(arr)
    Vector{dtype}(max.(arr, 0))
end

function min_spectrum2d(arr::AbstractMatrix{<:Integer})::Matrix{<:Unsigned}
    dtype = min_dtype(arr)
    Matrix{dtype}(max.(arr, 0))
end
