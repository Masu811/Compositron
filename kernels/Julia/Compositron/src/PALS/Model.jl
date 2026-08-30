using ..Utils

export LifetimeComponent, ResolutionComponent, LifetimeModel,
    FitParamInitializerList


mutable struct LifetimeComponent
    lifetime::FitParam
    intensity::FitParam
end


mutable struct ResolutionComponent
    fwhm::FitParam
    intensity::FitParam
    t0::FitParam
end


mutable struct LifetimeModel
    lifetime_components::Vector{LifetimeComponent}
    resolution_components::Vector{ResolutionComponent}
    background_component::FitParam
    scale_component::FitParam
    shift_component::FitParam
end


const FitParamInitializerList = Union{
    Tuple{<:Number},
    Tuple{<:Number, <:Number},
    Tuple{<:Number, Bool},
    Tuple{<:Number, <:Number, <:Number},
    Tuple{<:Number, Bool, <:Number},
    Tuple{<:Number, Bool, <:Number, <:Number}
}


function LifetimeModel(params::Pair{String, <:FitParamInitializerList}...)
    lifetime_components = Vector{LifetimeComponent}()
    n_l = 1

    model = Dict(params...)

    while true
        lt_key = "lifetime_$n_l"

        !haskey(model, lt_key) && break

        lifetime = FitParam(model[lt_key]...)

        i_key = "intensity_$n_l"

        intensity = if haskey(model, i_key)
            FitParam(model[i_key]...)
        else
            FitParam(0.5, NaN, 0, 1, true, true)
        end

        push!(lifetime_components, LifetimeComponent(lifetime, intensity))

        n_l += 1
    end

    resolution_components = Vector{ResolutionComponent}()
    n_r = 1

    while true
        fwhm_key = "res_fwhm_$n_r"

        !haskey(model, fwhm_key) && break

        fwhm = FitParam(model[fwhm_key]...)

        i_key = "res_intensity_$n_r"

        intensity = if haskey(model, i_key)
            FitParam(model[i_key]...)
        else
            FitParam(0.5, NaN, 0, 1, true, true)
        end

        t0_key = "res_to_$n_r"

        t0 = if haskey(model, t0_key)
            FitParam(model[t0_key]...)
        else
            FitParam(0., NaN, -Inf, Inf, true, true)
        end

        push!(resolution_components, ResolutionComponent(fwhm, intensity, t0))

        n_r += 1
    end

    if n_l == 1
        throw(ErrorException("Lifetime Model is missing lifetime components"))
    end
    if n_r == 1
        throw(ErrorException("Lifetime Model is missing resolution components"))
    end

    bg_key = "background"

    background_component = if haskey(model, bg_key)
        FitParam(model[bg_key]...)
    else
        FitParam(0., NaN, 0, Inf, false, true)
    end

    scale_key = "N"

    scale_component = if haskey(model, scale_key)
        FitParam(model[scale_key]...)
    else
        FitParam(1e7, NaN, 0, Inf, true, true)
    end

    shift_key = "t0"

    shift_component = if haskey(model, shift_key)
        FitParam(model[shift_key]...)
    else
        FitParam(0., NaN, -Inf, Inf, true, true)
    end

    LifetimeModel(
        lifetime_components,
        resolution_components,
        background_component,
        scale_component,
        shift_component
    )
end
