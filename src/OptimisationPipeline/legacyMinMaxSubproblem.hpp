#pragma once

#include "subproblem.hpp"
#include "../BoundSDP/boundSDP.hpp"
#include "../GlobalOptimiser/approxFunction.hpp"
#include "../GlobalOptimiser/scminmaxProb.hpp"
#include "../Laminate/laminateSection.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace lamopt {

class DynamicLinearApproximationFunction {
public:
    using Vector_r = Eigen::VectorXd;
    using Vector_v = Eigen::VectorXd;
    using Matrix_t = Eigen::MatrixXd;
    using Hessian_t = Eigen::MatrixXd;
    using Hessian_r = Eigen::MatrixXd;
    using Matrix_r_i = Eigen::MatrixXi;
    using Vector_r_i = Eigen::VectorXi;
    using reduced_Hessian_r = Eigen::MatrixXd;
    using reduced_Vector_r = Eigen::VectorXd;

    int NVAR = 0;
    int NRESP = 0;

    DynamicLinearApproximationFunction() = default;

    DynamicLinearApproximationFunction(const Eigen::Index responseCount, const Eigen::Index variableCount) {
        resize(responseCount, variableCount);
    }

    void resize(const Eigen::Index responseCount, const Eigen::Index variableCount) {
        NRESP = static_cast<int>(responseCount);
        NVAR = static_cast<int>(variableCount);
        freeTerms = Vector_r::Zero(responseCount);
        linearGradients = Matrix_t::Zero(variableCount, responseCount);
        objectiveMask = Vector_r::Zero(responseCount);
        objectiveCount = 0;
    }

    int getBooleanVector(Vector_r& booleanVector) {
        booleanVector = objectiveMask;
        return objectiveCount;
    }

    void ConfigureLinearModel(const Vector_v& referenceDesign,
                              const Vector_r& referenceResponses,
                              const Matrix_t& gradients,
                              const Vector_r& mask,
                              const int objectives) {
        resize(referenceResponses.size(), referenceDesign.size());
        linearGradients = gradients;
        objectiveMask = mask;
        objectiveCount = objectives;
        freeTerms = referenceResponses - gradients.transpose() * referenceDesign;
    }

    void Eval(Vector_v primalVar,
              Vector_r /*dualVar*/,
              Vector_r& responses,
              Matrix_t& gradients,
              Hessian_t& /*hessian*/) {
        responses = freeTerms + linearGradients.transpose() * primalVar;
        gradients = linearGradients;
    }

    void Eval(Vector_v primalVar, Vector_r& responses) {
        responses = freeTerms + linearGradients.transpose() * primalVar;
    }

private:
    Vector_r freeTerms;
    Matrix_t linearGradients;
    Vector_r objectiveMask;
    int objectiveCount = 0;
};

struct LegacyMinMax2Var4RespOptions {
    double defaultLowerBound = 1.0e-6;
    double strictLowerBoundFloor = 1.0e-9;
    double objectiveScaleFloor = 1.0e-12;
    bool verbose = false;
};

class LegacyMinMax2Var4RespSubproblemSolver : public SubproblemSolver {
public:
    explicit LegacyMinMax2Var4RespSubproblemSolver(LegacyMinMax2Var4RespOptions options = {})
        : m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        using ApproximationFunction = globopt::approxFunction<4, 2, double>;
        using SideConstraint = lampar::boundSDP<lampar::Lower, double>;

        SubproblemResult result;

        if (problem.referenceDesign.size() != 2) {
            result.message = "Legacy min-max adapter currently supports exactly 2 design variables.";
            return result;
        }
        if (problem.objectiveValues.size() != 1 || problem.constraintValues.size() != 3) {
            result.message = "Legacy min-max adapter currently supports 1 objective and 3 constraints.";
            return result;
        }
        if (!problem.objectiveGradients.has_value() || !problem.constraintGradients.has_value()) {
            result.message = "Legacy min-max adapter requires objective and constraint gradients.";
            return result;
        }
        if (problem.objectiveGradients->rows() != 2 || problem.objectiveGradients->cols() != 1) {
            result.message = "Objective gradient dimensions do not match the legacy adapter requirements.";
            return result;
        }
        if (problem.constraintGradients->rows() != 2 || problem.constraintGradients->cols() != 3) {
            result.message = "Constraint gradient dimensions do not match the legacy adapter requirements.";
            return result;
        }
        if (!problem.laminateSections.empty()) {
            result.message =
                "Laminate section side constraints are not yet wired into the legacy min-max adapter.";
            return result;
        }

        ApproximationFunction approximationFunction;
        typename ApproximationFunction::Vector_v design = problem.referenceDesign;
        typename ApproximationFunction::Vector_r referenceResponses;
        typename ApproximationFunction::Matrix_t gradients = ApproximationFunction::Matrix_t::Zero();
        typename ApproximationFunction::Vector_r objectiveMask;
        objectiveMask << 1.0, 0.0, 0.0, 0.0;

        const double objectiveScale =
            std::max(std::abs(problem.objectiveValues(0)), m_options.objectiveScaleFloor);
        referenceResponses << problem.objectiveValues(0) / objectiveScale,
                              problem.constraintValues(0),
                              problem.constraintValues(1),
                              problem.constraintValues(2);

        gradients.col(0) = problem.objectiveGradients->col(0) / objectiveScale;
        gradients.col(1) = problem.constraintGradients->col(0);
        gradients.col(2) = problem.constraintGradients->col(1);
        gradients.col(3) = problem.constraintGradients->col(2);

        approximationFunction.ConfigureConLinModel(design, referenceResponses, gradients, objectiveMask, 1);

        SideConstraint sideConstraints[2];
        for (int iVar = 0; iVar < 2; ++iVar) {
            double lowerBound = m_options.defaultLowerBound;
            if (problem.lowerBounds.has_value()) {
                if (problem.lowerBounds->size() != 2) {
                    result.message = "Lower bound vector size does not match the legacy adapter requirements.";
                    return result;
                }
                lowerBound = (*problem.lowerBounds)(iVar);
            }
            lowerBound = std::max(lowerBound, m_options.strictLowerBoundFloor);
            sideConstraints[iVar].setBound(lowerBound);
        }

        globopt::responseInterface<ApproximationFunction, SideConstraint> solver;
        solver.scMinMaxProb(&approximationFunction, sideConstraints);
        solver.setVerbose(m_options.verbose);
        solver.Solver(design);

        Eigen::VectorXd candidate = design;
        if (problem.lowerBounds.has_value()) {
            Eigen::VectorXd strictLowerBounds = *problem.lowerBounds;
            strictLowerBounds = strictLowerBounds.unaryExpr([this](double value) {
                return std::max(value, m_options.strictLowerBoundFloor);
            });
            candidate = candidate.cwiseMax(strictLowerBounds);
        } else {
            candidate = candidate.unaryExpr([this](double value) {
                return std::max(value, m_options.strictLowerBoundFloor);
            });
        }
        if (problem.upperBounds.has_value()) {
            if (problem.upperBounds->size() != 2) {
                result.message = "Upper bound vector size does not match the legacy adapter requirements.";
                return result;
            }
            candidate = candidate.cwiseMin(*problem.upperBounds);
        }

        typename ApproximationFunction::Vector_v clampedDesign = candidate;
        typename ApproximationFunction::Vector_r predictedResponses;
        approximationFunction.Eval(clampedDesign, predictedResponses);

        result.success = true;
        result.iterations = solver.getLastIterationCount();
        result.candidateDesign = candidate;
        result.predictedObjectives = Eigen::VectorXd::Constant(1, predictedResponses(0) * objectiveScale);
        result.predictedConstraints = Eigen::VectorXd(3);
        result.predictedConstraints << predictedResponses(1), predictedResponses(2), predictedResponses(3);
        result.message = "Candidate produced by the legacy scMinMaxProb adapter.";
        return result;
    }

private:
    LegacyMinMax2Var4RespOptions m_options;
};

struct LegacyLaminateSection1RespOptions {
    double defaultLowerBound = -10.0;
    double strictBoundMargin = 1.0e-9;
    double objectiveScaleFloor = 1.0e-12;
    bool verbose = false;
};

class LegacyLaminateSection1RespSubproblemSolver : public SubproblemSolver {
public:
    explicit LegacyLaminateSection1RespSubproblemSolver(LegacyLaminateSection1RespOptions options = {})
        : m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        using ApproximationFunction = DynamicLinearApproximationFunction;
        using SideConstraint = lampar::boundSDP<lampar::Lower, double>;

        SubproblemResult result;
        const Eigen::Index variableCount = problem.referenceDesign.size();
        const Eigen::Index constraintCount = problem.constraintValues.size();

        if (problem.objectiveValues.size() != 1) {
            result.message = "Legacy laminate adapter currently supports exactly 1 objective response.";
            return result;
        }
        if (!problem.objectiveGradients.has_value()) {
            result.message = "Legacy laminate adapter requires objective gradients.";
            return result;
        }
        if (problem.objectiveGradients->rows() != variableCount || problem.objectiveGradients->cols() != 1) {
            result.message = "Objective gradient dimensions do not match the laminate section adapter requirements.";
            return result;
        }
        if (constraintCount != 0 && !problem.constraintGradients.has_value()) {
            result.message = "Constraint gradients are required when response constraints are present.";
            return result;
        }
        if (problem.constraintGradients.has_value()) {
            if (problem.constraintGradients->rows() != variableCount
                || problem.constraintGradients->cols() != constraintCount) {
                result.message =
                    "Constraint gradient dimensions do not match the laminate section adapter requirements.";
                return result;
            }
        }
        if (problem.laminateSections.empty()) {
            result.message = "Legacy laminate adapter requires at least one laminate section block.";
            return result;
        }
        if (problem.lowerBounds.has_value() && problem.lowerBounds->size() != variableCount) {
            result.message = "Lower bound vector size does not match the laminate section adapter requirements.";
            return result;
        }
        if (problem.upperBounds.has_value() && problem.upperBounds->size() != variableCount) {
            result.message = "Upper bound vector size does not match the laminate section adapter requirements.";
            return result;
        }

        std::vector<SectionBinding> sectionBindings;
        sectionBindings.reserve(problem.laminateSections.size());
        std::vector<int> variableOwners(static_cast<size_t>(variableCount), -1);

        for (size_t iSection = 0; iSection < problem.laminateSections.size(); ++iSection) {
            SectionBuildResult sectionBuild = buildSectionBinding(problem.laminateSections[iSection]);
            if (!sectionBuild.success) {
                result.message = sectionBuild.message;
                return result;
            }

            const SectionBinding& binding = sectionBuild.binding;
            if (binding.offset < 0 || binding.offset + binding.size > variableCount) {
                result.message = "Laminate section slice exceeds the reference design size.";
                return result;
            }

            for (int iVar = 0; iVar < binding.size; ++iVar) {
                const int variableIndex = binding.offset + iVar;
                if (variableOwners[static_cast<size_t>(variableIndex)] != -1) {
                    result.message = "Laminate section blocks must not overlap in the design vector.";
                    return result;
                }
                variableOwners[static_cast<size_t>(variableIndex)] = static_cast<int>(iSection);
            }

            sectionBindings.push_back(std::move(sectionBuild.binding));
        }

        const Eigen::Index responseCount = std::max<Eigen::Index>(2, 1 + constraintCount);
        ApproximationFunction approximationFunction(responseCount, variableCount);
        ApproximationFunction::Vector_v design = problem.referenceDesign;
        ApproximationFunction::Vector_r referenceResponses =
            ApproximationFunction::Vector_r::Zero(responseCount);
        ApproximationFunction::Matrix_t gradients =
            ApproximationFunction::Matrix_t::Zero(variableCount, responseCount);
        ApproximationFunction::Vector_r objectiveMask =
            ApproximationFunction::Vector_r::Zero(responseCount);
        objectiveMask(0) = 1.0;

        const double objectiveScale =
            std::max(std::abs(problem.objectiveValues(0)), m_options.objectiveScaleFloor);
        referenceResponses(0) = problem.objectiveValues(0) / objectiveScale;
        gradients.col(0) = problem.objectiveGradients->col(0) / objectiveScale;
        if (constraintCount != 0) {
            referenceResponses.segment(1, constraintCount) = problem.constraintValues;
            gradients.block(0, 1, variableCount, constraintCount) = *problem.constraintGradients;
        } else {
            referenceResponses(1) = -1.0;
        }
        approximationFunction.ConfigureLinearModel(design, referenceResponses, gradients, objectiveMask, 1);

        std::vector<SideConstraint> sideConstraints(static_cast<size_t>(variableCount));
        for (Eigen::Index iVar = 0; iVar < variableCount; ++iVar) {
            double lowerBound = m_options.defaultLowerBound;
            if (problem.lowerBounds.has_value()) {
                lowerBound = (*problem.lowerBounds)(iVar);
            }
            lowerBound = std::min(lowerBound, design(iVar) - m_options.strictBoundMargin);
            sideConstraints[static_cast<size_t>(iVar)].setBound(lowerBound);
        }

        std::vector<optsection::section<double>*> sectionPointers;
        std::vector<int> sectionOffsets;
        sectionPointers.reserve(sectionBindings.size());
        sectionOffsets.reserve(sectionBindings.size());
        for (SectionBinding& binding : sectionBindings) {
            sectionPointers.push_back(binding.section.get());
            sectionOffsets.push_back(static_cast<int>(binding.offset));
        }

        globopt::responseInterface<ApproximationFunction, SideConstraint> solver;
        solver.scMinMaxProb(&approximationFunction, sideConstraints.data());
        solver.setSections(sectionPointers, sectionOffsets);
        solver.setVerbose(m_options.verbose);
        solver.Solver(design);

        Eigen::VectorXd candidate = design;
        if (problem.lowerBounds.has_value()) {
            candidate = candidate.cwiseMax(*problem.lowerBounds);
        }
        if (problem.upperBounds.has_value()) {
            candidate = candidate.cwiseMin(*problem.upperBounds);
        }

        ApproximationFunction::Vector_v clampedDesign = candidate;
        ApproximationFunction::Vector_r predictedResponses;
        approximationFunction.Eval(clampedDesign, predictedResponses);

        result.success = true;
        result.iterations = solver.getLastIterationCount();
        result.candidateDesign = candidate;
        result.predictedObjectives = Eigen::VectorXd::Constant(1, predictedResponses(0) * objectiveScale);
        result.predictedConstraints = constraintCount == 0
            ? Eigen::VectorXd()
            : predictedResponses.segment(1, constraintCount);
        result.message = "Candidate produced by the legacy laminate section adapter.";
        return result;
    }

private:
    struct SectionBinding {
        std::unique_ptr<optsection::section<double>> section;
        Eigen::Index offset = 0;
        Eigen::Index size = 0;
    };

    struct SectionBuildResult {
        bool success = false;
        SectionBinding binding;
        std::string message;
    };

    template<bool IsBalanced, bool IsSymmetric>
    SectionBuildResult dispatchBySubLaminates(const LaminateSectionState& sectionState) const {
        switch (sectionState.sublaminateCount) {
            case 1:
                return makeSectionBinding<IsBalanced, IsSymmetric, 1>(sectionState);
            case 2:
                return makeSectionBinding<IsBalanced, IsSymmetric, 2>(sectionState);
            case 4:
                return makeSectionBinding<IsBalanced, IsSymmetric, 4>(sectionState);
            case 5:
                return makeSectionBinding<IsBalanced, IsSymmetric, 5>(sectionState);
            case 6:
                return makeSectionBinding<IsBalanced, IsSymmetric, 6>(sectionState);
            case 8:
                return makeSectionBinding<IsBalanced, IsSymmetric, 8>(sectionState);
            case 10:
                return makeSectionBinding<IsBalanced, IsSymmetric, 10>(sectionState);
            default: {
                SectionBuildResult result;
                result.message = "Unsupported laminate section sublaminate count: "
                               + std::to_string(sectionState.sublaminateCount);
                return result;
            }
        }
    }

    template<bool IsBalanced, bool IsSymmetric, int NSUBLAM>
    SectionBuildResult makeSectionBinding(const LaminateSectionState& sectionState) const {
        using Section =
            lampar::laminateSection<lampar::SingleMaterial, IsBalanced, IsSymmetric, NSUBLAM, double>;

        SectionBuildResult result;
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

    SectionBuildResult buildSectionBinding(const LaminateSectionState& sectionState) const {
        if (sectionState.isBalanced) {
            if (sectionState.isSymmetric) {
                return dispatchBySubLaminates<true, true>(sectionState);
            }
            return dispatchBySubLaminates<true, false>(sectionState);
        }

        if (sectionState.isSymmetric) {
            return dispatchBySubLaminates<false, true>(sectionState);
        }
        return dispatchBySubLaminates<false, false>(sectionState);
    }

    LegacyLaminateSection1RespOptions m_options;
};

}  // namespace lamopt
