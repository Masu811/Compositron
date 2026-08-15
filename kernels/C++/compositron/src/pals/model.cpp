#include "compositron/pals/model.hpp"
#include "compositron/core/fitting.hpp"
#include <format>
#include <stdexcept>


namespace compositron::pals::model {


// class LifetimeModel

core::fitting::FitParam get_or_default(
    std::map<std::string, core::fitting::FitParam> model,
    std::string key,
    double val
) {
    if (model.contains(key)) {
        return model.at(key);
    } else {
        core::fitting::FitParam param(val);
        param.default_initialized = true;
        return param;
    }
}

core::fitting::FitParam get_or_default(
    std::map<std::string, core::fitting::FitParam> model,
    std::string key,
    double val,
    double min
) {
    if (model.contains(key)) {
        return model.at(key);
    } else {
        core::fitting::FitParam param(val, min);
        param.default_initialized = true;
        return param;
    }
}

core::fitting::FitParam get_or_default(
    std::map<std::string, core::fitting::FitParam> model,
    std::string key,
    double val,
    bool vary
) {
    if (model.contains(key)) {
        return model.at(key);
    } else {
        core::fitting::FitParam param(val, vary);
        param.default_initialized = true;
        return param;
    }
}

core::fitting::FitParam get_or_default(
    std::map<std::string, core::fitting::FitParam> model,
    std::string key,
    double val,
    double min,
    double max
) {
    if (model.contains(key)) {
        return model.at(key);
    } else {
        core::fitting::FitParam param(val, min, max);
        param.default_initialized = true;
        return param;
    }
}

LifetimeModel::LifetimeModel(
    std::map<std::string, core::fitting::FitParam> model
) {
    for (size_t n_l = 1;; ++n_l) {
        std::string lt_key = std::format("lifetime_{}", n_l);
        std::string i_key = std::format("intensity_{}", n_l);

        if (!model.contains(lt_key)) break;

        lifetime_components.push_back(
            LifetimeComponent(
                model.at(lt_key),
                get_or_default(model, i_key, 0.5, 0., 1.)
            )
        );
    }

    for (size_t n_r = 1;; ++n_r) {
        std::string fwhm_key = std::format("res_fwhm_{}", n_r);
        std::string i_key = std::format("res_intensity_{}", n_r);
        std::string t0_key = std::format("res_t0_{}", n_r);

        if (!model.contains(fwhm_key)) break;

        resolution_components.push_back(
            ResolutionComponent(
                model.at(fwhm_key),
                get_or_default(model, i_key, 0.5, 0., 1.),
                get_or_default(model, t0_key, 0.)
            )
        );
    }

    background_component = get_or_default(model, "background", 0., false);
    scale_component = get_or_default(model, "N", 1e7, 0.);
    shift_component = get_or_default(model, "t0", 0.);
}


} // namespace compositron::pals::model
