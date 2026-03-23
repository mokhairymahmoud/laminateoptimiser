#pragma once

#include "analysis.hpp"
#include "../Laminate/laminateSection.hpp"
#include "../Section/section.hpp"

#include <memory>
#include <string>
#include <vector>

namespace lamopt {

struct LaminateSectionBinding {
    std::unique_ptr<optsection::section<double>> section;
    Eigen::Index offset = 0;
    Eigen::Index size = 0;
};

struct LaminateSectionBuildResult {
    bool success = false;
    LaminateSectionBinding binding;
    std::string message;
};

namespace detail {

template<bool IsBalanced, bool IsSymmetric, int NSUBLAM>
LaminateSectionBuildResult MakeLaminateSectionBinding(const LaminateSectionState& sectionState) {
    using Section =
        lampar::laminateSection<lampar::SingleMaterial, IsBalanced, IsSymmetric, NSUBLAM, double>;

    LaminateSectionBuildResult result;
    if (sectionState.laminationParameters.size() != 0
        && sectionState.laminationParameters.size() != Section::Size - 1) {
        result.message = "Laminate section lamination-parameter size does not match the selected section type.";
        return result;
    }

    auto section = std::make_unique<Section>();
    section->setBoundThickness(sectionState.thicknessLowerBound, sectionState.thicknessUpperBound);

    result.success = true;
    result.binding.section = std::move(section);
    result.binding.offset = sectionState.variableOffset;
    result.binding.size = Section::Size;
    return result;
}

template<bool IsBalanced, bool IsSymmetric>
LaminateSectionBuildResult DispatchLaminateSectionBinding(const LaminateSectionState& sectionState) {
    switch (sectionState.sublaminateCount) {
        case 1:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 1>(sectionState);
        case 2:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 2>(sectionState);
        case 4:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 4>(sectionState);
        case 5:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 5>(sectionState);
        case 6:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 6>(sectionState);
        case 8:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 8>(sectionState);
        case 10:
            return MakeLaminateSectionBinding<IsBalanced, IsSymmetric, 10>(sectionState);
        default: {
            LaminateSectionBuildResult result;
            result.message = "Unsupported laminate section sublaminate count: "
                           + std::to_string(sectionState.sublaminateCount);
            return result;
        }
    }
}

}  // namespace detail

inline LaminateSectionBuildResult BuildLaminateSectionBinding(const LaminateSectionState& sectionState) {
    if (sectionState.sublaminateCount <= 0) {
        LaminateSectionBuildResult result;
        result.message = "Laminate section state must define a positive sublaminate count.";
        return result;
    }

    if (sectionState.isBalanced) {
        if (sectionState.isSymmetric) {
            return detail::DispatchLaminateSectionBinding<true, true>(sectionState);
        }
        return detail::DispatchLaminateSectionBinding<true, false>(sectionState);
    }

    if (sectionState.isSymmetric) {
        return detail::DispatchLaminateSectionBinding<false, true>(sectionState);
    }
    return detail::DispatchLaminateSectionBinding<false, false>(sectionState);
}

inline bool BuildLaminateSectionBindings(const std::vector<LaminateSectionState>& sectionStates,
                                         const Eigen::Index variableCount,
                                         std::vector<LaminateSectionBinding>& bindings,
                                         std::string& message) {
    bindings.clear();
    bindings.reserve(sectionStates.size());

    std::vector<int> variableOwners(static_cast<size_t>(variableCount), -1);

    for (size_t iSection = 0; iSection < sectionStates.size(); ++iSection) {
        LaminateSectionBuildResult sectionBuild = BuildLaminateSectionBinding(sectionStates[iSection]);
        if (!sectionBuild.success) {
            message = sectionBuild.message;
            return false;
        }

        const LaminateSectionBinding& binding = sectionBuild.binding;
        if (binding.offset < 0 || binding.offset + binding.size > variableCount) {
            message = "Laminate section slice exceeds the reference design size.";
            return false;
        }

        for (int iVar = 0; iVar < binding.size; ++iVar) {
            const int variableIndex = static_cast<int>(binding.offset) + iVar;
            if (variableOwners[static_cast<size_t>(variableIndex)] != -1) {
                message = "Laminate section blocks must not overlap in the design vector.";
                return false;
            }
            variableOwners[static_cast<size_t>(variableIndex)] = static_cast<int>(iSection);
        }

        bindings.push_back(std::move(sectionBuild.binding));
    }

    return true;
}

inline void ExtractLaminateSectionPointersAndOffsets(
    std::vector<LaminateSectionBinding>& bindings,
    std::vector<optsection::section<double>*>& sectionPointers,
    std::vector<int>& sectionOffsets) {
    sectionPointers.clear();
    sectionOffsets.clear();
    sectionPointers.reserve(bindings.size());
    sectionOffsets.reserve(bindings.size());

    for (LaminateSectionBinding& binding : bindings) {
        sectionPointers.push_back(binding.section.get());
        sectionOffsets.push_back(static_cast<int>(binding.offset));
    }
}

}  // namespace lamopt
