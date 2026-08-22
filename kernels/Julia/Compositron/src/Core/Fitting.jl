struct SimpleFitParam
    val::Float64
    err::Float64
end

mutable struct FitParam
    val::Float64
    err::Float64
    min::Float64
    max::Float64
    vary::Bool
    default_initialized::Bool
end

function FitParam(val::Number)
    FitParam(Float64(val), NaN, -Inf, Inf, true, false)
end

function FitParam(val::Number, min::Number)
    FitParam(Float64(val), NaN, Float64(min), Inf, true, false)
end

function FitParam(val::Number, vary::Bool)
    FitParam(Float64(val), NaN, -Inf, Inf, vary, false)
end

function FitParam(val::Number, vary::Bool, min::Number)
    FitParam(Float64(val), NaN, Float64(min), Inf, vary, false)
end

function FitParam(val::Number, min::Number, max::Number)
    FitParam(Float64(val), NaN, Float64(min), Float64(max), true, false)
end

function FitParam(val::Number, vary::Bool, min::Number, max::Number)
    FitParam(Float64(val), NaN, Float64(min), Float64(max), vary, false)
end
