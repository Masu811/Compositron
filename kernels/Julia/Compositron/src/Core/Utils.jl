import Base: copy

using ..Constants: KEV_PER_M0C

export LinearCalibration, from_index, to_index_rounded, to_index_float,
    to_index_trunc, EcalCorrectionOrder, EcalCorrOrder_Zeroth,
    EcalCorrOrder_First, EcalCorrOrder_None, Unit, KeV, M0C, Eres, to_kev,
    EnergyDetector, EnergyDetectorPair, TimingDetectorPair,
    min_dtype, min_spectrum, min_spectrum2d


# LinearCalibration

mutable struct LinearCalibration
    offset::Float64
    scale::Float64
end


function copy(ecal::LinearCalibration)::LinearCalibration
    LinearCalibration(ecal.offset, ecal.scale)
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

@enum EcalCorrectionOrder begin
    EcalCorrOrder_Zeroth
    EcalCorrOrder_First
    EcalCorrOrder_None
end


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


function to_kev(value::KeV, eres::Union{Nothing, Float64})::Float64
    value.value
end


function to_kev(value::M0C, eres::Union{Nothing, Float64})::Float64
    value.value * KEV_PER_M0C
end


function to_kev(value::Eres, eres::Union{Nothing, Float64})::Float64
    if isnothing(eres)
        throw(ErrorException(
            "Could not convert units of eres to keV due to missing eres"
        ))
    end

    value.value * eres
end


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
