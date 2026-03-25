#pragma once

#include "analysis.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lamopt {

enum class LaminateTensorBlockKind {
    A,
    D,
    B
};

struct CanonicalLaminateVariableSlot {
    Eigen::Index globalIndex = 0;
    Eigen::Index sectionLocalIndex = 0;
    LaminateTensorBlockKind blockKind = LaminateTensorBlockKind::A;
    Eigen::Index laminationParameterIndex = 0;
};

struct CanonicalLaminateSectionLayout {
    std::size_t sectionIndex = 0;
    Eigen::Index variableOffset = 0;
    Eigen::Index sectionSize = 0;
    Eigen::Index laminationParameterCount = 0;
    Eigen::Index thicknessIndex = 0;
    bool isBalanced = false;
    bool isSymmetric = true;
    int sublaminateCount = 0;
    std::vector<CanonicalLaminateVariableSlot> laminationParameterSlots;
};

struct LaminateTemplateParameterisationOptions {
    std::string sectionTokenPrefix = "SEC";
    bool includeSectionMetadata = true;
};

[[nodiscard]] inline Eigen::Index CanonicalLaminationParameterCount(const LaminateSectionState& sectionState) {
    return sectionState.isBalanced ? 2 : 4;
}

[[nodiscard]] inline Eigen::Index CanonicalLaminateTensorBlockCount(const LaminateSectionState& sectionState) {
    return sectionState.isSymmetric ? 2 : 3;
}

[[nodiscard]] inline Eigen::Index CanonicalLaminateSectionVariableCount(const LaminateSectionState& sectionState) {
    return CanonicalLaminationParameterCount(sectionState) * CanonicalLaminateTensorBlockCount(sectionState) + 1;
}

[[nodiscard]] inline std::vector<LaminateTensorBlockKind>
CanonicalLaminateTensorBlockOrder(const LaminateSectionState& sectionState) {
    std::vector<LaminateTensorBlockKind> blocks = {
        LaminateTensorBlockKind::A,
        LaminateTensorBlockKind::D
    };
    if (!sectionState.isSymmetric) {
        blocks.push_back(LaminateTensorBlockKind::B);
    }
    return blocks;
}

[[nodiscard]] inline const char* LaminateTensorBlockTokenName(const LaminateTensorBlockKind blockKind) {
    switch (blockKind) {
        case LaminateTensorBlockKind::A:
            return "A";
        case LaminateTensorBlockKind::D:
            return "D";
        case LaminateTensorBlockKind::B:
            return "B";
    }
    return "UNKNOWN";
}

inline bool BuildCanonicalLaminateSectionLayout(const AnalysisRequest& request,
                                                const std::size_t sectionIndex,
                                                CanonicalLaminateSectionLayout& layout,
                                                std::string& message) {
    if (sectionIndex >= request.laminateSections.size()) {
        message = "Laminate section index is out of bounds for canonical layout construction.";
        return false;
    }

    const LaminateSectionState& sectionState = request.laminateSections[sectionIndex];
    if (sectionState.sublaminateCount <= 0) {
        message = "Laminate section state must define a positive sublaminate count.";
        return false;
    }

    const Eigen::Index laminationParameterCount = CanonicalLaminationParameterCount(sectionState);
    const Eigen::Index sectionSize = CanonicalLaminateSectionVariableCount(sectionState);
    if (sectionState.laminationParameters.size() != 0
        && sectionState.laminationParameters.size() != sectionSize - 1) {
        message = "Laminate section lamination-parameter size does not match the canonical section layout.";
        return false;
    }
    if (sectionState.variableOffset < 0
        || sectionState.variableOffset + sectionSize > request.designVariables.size()) {
        message = "Laminate section slice exceeds the design vector in the canonical section layout.";
        return false;
    }

    layout = {};
    layout.sectionIndex = sectionIndex;
    layout.variableOffset = sectionState.variableOffset;
    layout.sectionSize = sectionSize;
    layout.laminationParameterCount = laminationParameterCount;
    layout.thicknessIndex = sectionState.variableOffset + sectionSize - 1;
    layout.isBalanced = sectionState.isBalanced;
    layout.isSymmetric = sectionState.isSymmetric;
    layout.sublaminateCount = sectionState.sublaminateCount;

    const std::vector<LaminateTensorBlockKind> blocks = CanonicalLaminateTensorBlockOrder(sectionState);
    layout.laminationParameterSlots.reserve(static_cast<std::size_t>(sectionSize - 1));

    Eigen::Index localIndex = 0;
    for (const LaminateTensorBlockKind blockKind : blocks) {
        for (Eigen::Index iParameter = 0; iParameter < laminationParameterCount; ++iParameter) {
            layout.laminationParameterSlots.push_back({
                sectionState.variableOffset + localIndex,
                localIndex,
                blockKind,
                iParameter
            });
            ++localIndex;
        }
    }

    return true;
}

inline bool BuildCanonicalLaminateSectionLayouts(const AnalysisRequest& request,
                                                 std::vector<CanonicalLaminateSectionLayout>& layouts,
                                                 std::string& message) {
    layouts.clear();
    layouts.reserve(request.laminateSections.size());

    std::vector<int> variableOwners(static_cast<std::size_t>(request.designVariables.size()), -1);
    for (std::size_t sectionIndex = 0; sectionIndex < request.laminateSections.size(); ++sectionIndex) {
        CanonicalLaminateSectionLayout layout;
        if (!BuildCanonicalLaminateSectionLayout(request, sectionIndex, layout, message)) {
            return false;
        }

        for (Eigen::Index iVar = 0; iVar < layout.sectionSize; ++iVar) {
            const Eigen::Index globalIndex = layout.variableOffset + iVar;
            if (variableOwners[static_cast<std::size_t>(globalIndex)] != -1) {
                message = "Laminate section blocks must not overlap in the canonical section layout.";
                return false;
            }
            variableOwners[static_cast<std::size_t>(globalIndex)] = static_cast<int>(sectionIndex);
        }

        layouts.push_back(std::move(layout));
    }

    return true;
}

[[nodiscard]] inline std::string FormatLaminateTemplateDouble(const double value) {
    std::ostringstream stream;
    stream.setf(std::ios::scientific);
    stream.precision(16);
    stream << value;
    return stream.str();
}

[[nodiscard]] inline std::string MakeLaminateSectionToken(const LaminateTemplateParameterisationOptions& options,
                                                          const std::size_t sectionIndex,
                                                          const std::string& suffix) {
    return "{{" + options.sectionTokenPrefix + std::to_string(sectionIndex) + "_" + suffix + "}}";
}

inline bool BuildLaminateTemplateParameters(const AnalysisRequest& request,
                                            std::vector<ResolvedTemplateParameter>& parameters,
                                            std::string& message,
                                            const LaminateTemplateParameterisationOptions& options = {}) {
    std::vector<CanonicalLaminateSectionLayout> layouts;
    if (!BuildCanonicalLaminateSectionLayouts(request, layouts, message)) {
        return false;
    }

    parameters.clear();
    for (const CanonicalLaminateSectionLayout& layout : layouts) {
        if (options.includeSectionMetadata) {
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "OFFSET"),
                                  std::to_string(layout.variableOffset)});
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "SIZE"),
                                  std::to_string(layout.sectionSize)});
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "LP_COUNT"),
                                  std::to_string(layout.laminationParameterCount)});
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "SUBLAMINATE_COUNT"),
                                  std::to_string(layout.sublaminateCount)});
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "IS_BALANCED"),
                                  layout.isBalanced ? "1" : "0"});
            parameters.push_back({MakeLaminateSectionToken(options, layout.sectionIndex, "IS_SYMMETRIC"),
                                  layout.isSymmetric ? "1" : "0"});
        }

        parameters.push_back({
            MakeLaminateSectionToken(options, layout.sectionIndex, "THICKNESS"),
            FormatLaminateTemplateDouble(request.designVariables(layout.thicknessIndex))
        });

        for (const CanonicalLaminateVariableSlot& slot : layout.laminationParameterSlots) {
            parameters.push_back({
                MakeLaminateSectionToken(options,
                                         layout.sectionIndex,
                                         std::string(LaminateTensorBlockTokenName(slot.blockKind))
                                             + "_LP" + std::to_string(slot.laminationParameterIndex)),
                FormatLaminateTemplateDouble(request.designVariables(slot.globalIndex))
            });
        }
    }

    return true;
}

[[nodiscard]] inline std::vector<ResolvedTemplateParameter>
MakeLaminateTemplateParameters(const AnalysisRequest& request,
                               const LaminateTemplateParameterisationOptions& options = {}) {
    std::vector<ResolvedTemplateParameter> parameters;
    std::string message;
    if (!BuildLaminateTemplateParameters(request, parameters, message, options)) {
        throw std::runtime_error(message);
    }
    return parameters;
}

}  // namespace lamopt
