#pragma once

#include "abaqusJobBackend.hpp"
#include "calculixJobBackend.hpp"

#include <memory>

namespace lamopt {

enum class AnalysisBackendKind {
    Calculix,
    AbaqusCompatibility
};

struct DefaultAnalysisBackendSetup {
    AnalysisBackendKind backendKind = AnalysisBackendKind::Calculix;
    std::optional<CalculixJobSetup> calculix;
    std::optional<AbaqusJobConfig> abaqusCompatibility;
};

[[nodiscard]] inline std::unique_ptr<AnalysisBackend>
MakeDefaultAnalysisBackend(const CalculixJobSetup& setup) {
    return std::make_unique<CalculixJobBackend>(MakeDefaultCalculixJobConfig(setup));
}

[[nodiscard]] inline std::unique_ptr<AnalysisBackend>
MakeConfiguredAnalysisBackend(const DefaultAnalysisBackendSetup& setup) {
    switch (setup.backendKind) {
        case AnalysisBackendKind::Calculix:
            if (!setup.calculix.has_value()) {
                throw std::runtime_error("CalculiX backend setup is required for the default analysis backend.");
            }
            return MakeDefaultAnalysisBackend(*setup.calculix);
        case AnalysisBackendKind::AbaqusCompatibility:
            if (!setup.abaqusCompatibility.has_value()) {
                throw std::runtime_error("Abaqus compatibility setup is required for the configured analysis backend.");
            }
            return std::make_unique<AbaqusJobBackend>(*setup.abaqusCompatibility);
    }

    throw std::runtime_error("Unknown analysis backend kind.");
}

}  // namespace lamopt
