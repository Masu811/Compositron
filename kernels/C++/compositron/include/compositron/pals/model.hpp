#pragma once

#include <vector>

#include "compositron/core/fitting.hpp"


namespace compositron::pals::model {


struct LifetimeComponent {
    core::fitting::FitParam lifetime;
    core::fitting::FitParam intensity;
};


struct ResolutionComponent {
    core::fitting::FitParam fwhm;
    core::fitting::FitParam intensity;
    core::fitting::FitParam t0;
};


class LifetimeModel {
public:
    std::vector<LifetimeComponent> lifetime_components;
    std::vector<ResolutionComponent> resolution_components;
    core::fitting::FitParam background_component;
    core::fitting::FitParam scale_component;
    core::fitting::FitParam shift_component;

    LifetimeModel(std::map<std::string, core::fitting::FitParam> model);
};


} // namespace compositron::pals::model
