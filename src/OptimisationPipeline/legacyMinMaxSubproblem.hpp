#pragma once

#include "subproblem.hpp"
#include "../BoundSDP/boundSDP.hpp"
#include "../GlobalOptimiser/approxFunction.hpp"
#include "../GlobalOptimiser/scminmaxProb.hpp"
#include "../Laminate/laminateSection.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace lamopt {

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
        SubproblemResult result;

        if (problem.objectiveValues.size() != 1) {
            result.message = "Legacy laminate adapter currently supports exactly 1 objective response.";
            return result;
        }
        if (problem.constraintValues.size() != 0) {
            result.message = "Legacy laminate adapter does not yet support response constraints.";
            return result;
        }
        if (!problem.objectiveGradients.has_value()) {
            result.message = "Legacy laminate adapter requires objective gradients.";
            return result;
        }
        if (problem.laminateSections.size() != 1) {
            result.message = "Legacy laminate adapter currently supports exactly one laminate section block.";
            return result;
        }

        const LaminateSectionState& sectionState = problem.laminateSections.front();
        if (sectionState.variableOffset != 0) {
            result.message = "Legacy laminate adapter currently requires the laminate section to start at offset 0.";
            return result;
        }

        if (sectionState.isBalanced) {
            if (sectionState.isSymmetric) {
                return dispatchBySubLaminates<true, true>(problem, sectionState);
            }
            return dispatchBySubLaminates<true, false>(problem, sectionState);
        }

        if (sectionState.isSymmetric) {
            return dispatchBySubLaminates<false, true>(problem, sectionState);
        }
        return dispatchBySubLaminates<false, false>(problem, sectionState);
    }

private:
    template<bool IsBalanced, bool IsSymmetric>
    SubproblemResult dispatchBySubLaminates(const ApproximationProblem& problem,
                                            const LaminateSectionState& sectionState) const {
        switch (sectionState.sublaminateCount) {
            case 1:
                return solveSection<IsBalanced, IsSymmetric, 1>(problem, sectionState);
            case 2:
                return solveSection<IsBalanced, IsSymmetric, 2>(problem, sectionState);
            case 4:
                return solveSection<IsBalanced, IsSymmetric, 4>(problem, sectionState);
            case 5:
                return solveSection<IsBalanced, IsSymmetric, 5>(problem, sectionState);
            case 6:
                return solveSection<IsBalanced, IsSymmetric, 6>(problem, sectionState);
            case 8:
                return solveSection<IsBalanced, IsSymmetric, 8>(problem, sectionState);
            case 10:
                return solveSection<IsBalanced, IsSymmetric, 10>(problem, sectionState);
            default: {
                SubproblemResult result;
                result.message = "Unsupported laminate section sublaminate count: "
                               + std::to_string(sectionState.sublaminateCount);
                return result;
            }
        }
    }

    template<bool IsBalanced, bool IsSymmetric, int NSUBLAM>
    SubproblemResult solveSection(const ApproximationProblem& problem,
                                  const LaminateSectionState& sectionState) const {
        using ApproximationFunction =
            globopt::approxFunction<2,
                                    lampar::laminateSection<lampar::SingleMaterial,
                                                            IsBalanced,
                                                            IsSymmetric,
                                                            NSUBLAM,
                                                            double>::Size,
                                    double>;
        using SideConstraint = lampar::boundSDP<lampar::Lower, double>;
        using Section =
            lampar::laminateSection<lampar::SingleMaterial, IsBalanced, IsSymmetric, NSUBLAM, double>;

        SubproblemResult result;
        constexpr int kSectionSize = Section::Size;

        if (problem.referenceDesign.size() != kSectionSize) {
            result.message = "Reference design size does not match the laminate section adapter requirements.";
            return result;
        }
        if (problem.objectiveGradients->rows() != kSectionSize || problem.objectiveGradients->cols() != 1) {
            result.message = "Objective gradient dimensions do not match the laminate section adapter requirements.";
            return result;
        }

        ApproximationFunction approximationFunction;
        typename ApproximationFunction::Vector_v design = problem.referenceDesign;
        typename ApproximationFunction::Vector_r referenceResponses;
        typename ApproximationFunction::Matrix_t gradients = ApproximationFunction::Matrix_t::Zero();
        typename ApproximationFunction::Vector_r objectiveMask;
        objectiveMask << 1.0, 0.0;

        const double objectiveScale =
            std::max(std::abs(problem.objectiveValues(0)), m_options.objectiveScaleFloor);
        referenceResponses << problem.objectiveValues(0) / objectiveScale, -1.0;
        gradients.col(0) = problem.objectiveGradients->col(0) / objectiveScale;
        approximationFunction.ConfigureLinearModel(design, referenceResponses, gradients, objectiveMask, 1);

        SideConstraint sideConstraints[kSectionSize];
        for (int iVar = 0; iVar < kSectionSize; ++iVar) {
            double lowerBound = m_options.defaultLowerBound;
            if (problem.lowerBounds.has_value()) {
                if (problem.lowerBounds->size() != kSectionSize) {
                    result.message =
                        "Lower bound vector size does not match the laminate section adapter requirements.";
                    return result;
                }
                lowerBound = (*problem.lowerBounds)(iVar);
            }
            lowerBound = std::min(lowerBound, design(iVar) - m_options.strictBoundMargin);
            sideConstraints[iVar].setBound(lowerBound);
        }

        Section laminateSection;
        laminateSection.setBoundThickness(sectionState.thicknessLowerBound, sectionState.thicknessUpperBound);

        globopt::responseInterface<ApproximationFunction, SideConstraint> solver;
        solver.scMinMaxProb(&approximationFunction, sideConstraints);
        solver.setSections({static_cast<optsection::section<double>*>(&laminateSection)}, {0});
        solver.setVerbose(m_options.verbose);
        solver.Solver(design);

        Eigen::VectorXd candidate = design;
        if (problem.lowerBounds.has_value()) {
            candidate = candidate.cwiseMax(*problem.lowerBounds);
        }
        if (problem.upperBounds.has_value()) {
            if (problem.upperBounds->size() != kSectionSize) {
                result.message =
                    "Upper bound vector size does not match the laminate section adapter requirements.";
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
        result.predictedConstraints = Eigen::VectorXd();
        result.message = "Candidate produced by the legacy laminate section adapter.";
        return result;
    }

    LegacyLaminateSection1RespOptions m_options;
};

}  // namespace lamopt
