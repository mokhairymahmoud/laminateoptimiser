#pragma once

#include "laminateSectionBinding.hpp"
#include "subproblem.hpp"
#include "../BoundSDP/boundSDP.hpp"
#include "../GlobalOptimiser/approxFunction.hpp"
#include "../GlobalOptimiser/scminmaxProb.hpp"

#include <algorithm>
#include <cmath>
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
        quadraticCurvature = Matrix_t::Zero(variableCount, responseCount);
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
        quadraticCurvature = Matrix_t::Zero(referenceDesign.size(), referenceResponses.size());
        objectiveMask = mask;
        objectiveCount = objectives;
        freeTerms = referenceResponses - gradients.transpose() * referenceDesign;
    }

    void ConfigureQuadraticModel(const Vector_v& referenceDesign,
                                 const Vector_r& referenceResponses,
                                 const Matrix_t& gradients,
                                 const Matrix_t& curvature,
                                 const Vector_r& mask,
                                 const int objectives) {
        resize(referenceResponses.size(), referenceDesign.size());
        linearGradients = gradients;
        quadraticCurvature = curvature;
        objectiveMask = mask;
        objectiveCount = objectives;

        for (Eigen::Index iResp = 0; iResp < referenceResponses.size(); ++iResp) {
            double freeTerm = referenceResponses(iResp);
            for (Eigen::Index iVar = 0; iVar < referenceDesign.size(); ++iVar) {
                const double x0 = referenceDesign(iVar);
                freeTerm -= gradients(iVar, iResp) * x0;
                freeTerm += 0.5 * curvature(iVar, iResp) * x0 * x0;
                linearGradients(iVar, iResp) -= curvature(iVar, iResp) * x0;
            }
            freeTerms(iResp) = freeTerm;
        }
    }

    void Eval(Vector_v primalVar,
              Vector_r dualVar,
              Vector_r& responses,
              Matrix_t& gradients,
              Hessian_t& hessian) {
        responses = freeTerms + linearGradients.transpose() * primalVar;
        gradients = linearGradients;
        if (!quadraticCurvature.isZero(0.0)) {
            for (Eigen::Index iResp = 0; iResp < NRESP; ++iResp) {
                for (Eigen::Index iVar = 0; iVar < NVAR; ++iVar) {
                    const double curvature = quadraticCurvature(iVar, iResp);
                    responses(iResp) += 0.5 * curvature * primalVar(iVar) * primalVar(iVar);
                    gradients(iVar, iResp) += curvature * primalVar(iVar);
                    hessian(iVar, iVar) += curvature * dualVar(iResp);
                }
            }
        }
    }

    void Eval(Vector_v primalVar, Vector_r& responses) {
        responses = freeTerms + linearGradients.transpose() * primalVar;
        if (!quadraticCurvature.isZero(0.0)) {
            for (Eigen::Index iResp = 0; iResp < NRESP; ++iResp) {
                for (Eigen::Index iVar = 0; iVar < NVAR; ++iVar) {
                    responses(iResp) += 0.5 * quadraticCurvature(iVar, iResp) * primalVar(iVar) * primalVar(iVar);
                }
            }
        }
    }

private:
    Vector_r freeTerms;
    Matrix_t linearGradients;
    Matrix_t quadraticCurvature;
    Vector_r objectiveMask;
    int objectiveCount = 0;
};

struct CoreMinMax2Var4RespOptions {
    double defaultLowerBound = 1.0e-6;
    double strictLowerBoundFloor = 1.0e-9;
    double objectiveScaleFloor = 1.0e-12;
    bool verbose = false;
};

class CoreMinMax2Var4RespSubproblemSolver : public SubproblemSolver {
public:
    explicit CoreMinMax2Var4RespSubproblemSolver(CoreMinMax2Var4RespOptions options = {})
        : m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        using ApproximationFunction = globopt::approxFunction<4, 2, double>;
        using SideConstraint = lampar::boundSDP<lampar::Lower, double>;

        SubproblemResult result;

        if (problem.referenceDesign.size() != 2) {
            result.message = "Core min-max adapter currently supports exactly 2 design variables.";
            return result;
        }
        if (problem.objectiveValues.size() != 1 || problem.constraintValues.size() != 3) {
            result.message = "Core min-max adapter currently supports 1 objective and 3 constraints.";
            return result;
        }
        if (!problem.objectiveGradients.has_value() || !problem.constraintGradients.has_value()) {
            result.message = "Core min-max adapter requires objective and constraint gradients.";
            return result;
        }
        if (problem.objectiveGradients->rows() != 2 || problem.objectiveGradients->cols() != 1) {
            result.message = "Objective gradient dimensions do not match the core adapter requirements.";
            return result;
        }
        if (problem.constraintGradients->rows() != 2 || problem.constraintGradients->cols() != 3) {
            result.message = "Constraint gradient dimensions do not match the core adapter requirements.";
            return result;
        }
        if (!problem.laminateSections.empty()) {
            result.message =
                "Laminate section side constraints are not yet wired into the core min-max benchmark adapter.";
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
                    result.message = "Lower bound vector size does not match the core adapter requirements.";
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
                result.message = "Upper bound vector size does not match the core adapter requirements.";
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
        result.message = "Candidate produced by the core scMinMaxProb adapter.";
        return result;
    }

private:
    CoreMinMax2Var4RespOptions m_options;
};

struct CoreLaminateSection1RespOptions {
    double defaultLowerBound = -10.0;
    double strictBoundMargin = 1.0e-9;
    double objectiveScaleFloor = 1.0e-12;
    bool verbose = false;
};

class CoreLaminateSection1RespSubproblemSolver : public SubproblemSolver {
public:
    explicit CoreLaminateSection1RespSubproblemSolver(CoreLaminateSection1RespOptions options = {})
        : m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        using ApproximationFunction = DynamicLinearApproximationFunction;
        using SideConstraint = lampar::boundSDP<lampar::Lower, double>;

        SubproblemResult result;
        const Eigen::Index variableCount = problem.referenceDesign.size();
        const Eigen::Index objectiveCount = problem.objectiveValues.size();
        const Eigen::Index constraintCount = problem.constraintValues.size();

        if (objectiveCount == 0) {
            result.message = "Core laminate adapter requires at least one objective response.";
            return result;
        }
        if (!problem.objectiveGradients.has_value()) {
            result.message = "Core laminate adapter requires objective gradients.";
            return result;
        }
        if (problem.objectiveGradients->rows() != variableCount || problem.objectiveGradients->cols() != objectiveCount) {
            result.message = "Objective gradient dimensions do not match the laminate section adapter requirements.";
            return result;
        }
        if (problem.objectiveCurvature.has_value()) {
            if (problem.objectiveCurvature->rows() != variableCount
                || problem.objectiveCurvature->cols() != objectiveCount) {
                result.message = "Objective curvature dimensions do not match the laminate section adapter requirements.";
                return result;
            }
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
            result.message = "Core laminate adapter requires at least one laminate section block.";
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

        std::vector<LaminateSectionBinding> sectionBindings;
        if (!BuildLaminateSectionBindings(problem.laminateSections,
                                          variableCount,
                                          sectionBindings,
                                          result.message)) {
            return result;
        }

        const Eigen::Index responseCount = std::max<Eigen::Index>(2, objectiveCount + constraintCount);
        ApproximationFunction approximationFunction(responseCount, variableCount);
        ApproximationFunction::Vector_v design = problem.referenceDesign;
        ApproximationFunction::Vector_r referenceResponses =
            ApproximationFunction::Vector_r::Zero(responseCount);
        ApproximationFunction::Matrix_t gradients =
            ApproximationFunction::Matrix_t::Zero(variableCount, responseCount);
        ApproximationFunction::Matrix_t curvature =
            ApproximationFunction::Matrix_t::Zero(variableCount, responseCount);
        ApproximationFunction::Vector_r objectiveMask =
            ApproximationFunction::Vector_r::Zero(responseCount);
        objectiveMask.head(objectiveCount).setOnes();

        Eigen::VectorXd objectiveScales = Eigen::VectorXd::Constant(objectiveCount, m_options.objectiveScaleFloor);
        for (Eigen::Index iObjective = 0; iObjective < objectiveCount; ++iObjective) {
            objectiveScales(iObjective) =
                std::max(std::abs(problem.objectiveValues(iObjective)), m_options.objectiveScaleFloor);
            referenceResponses(iObjective) = problem.objectiveValues(iObjective) / objectiveScales(iObjective);
            gradients.col(iObjective) = problem.objectiveGradients->col(iObjective) / objectiveScales(iObjective);
            if (problem.objectiveCurvature.has_value()) {
                curvature.col(iObjective) = problem.objectiveCurvature->col(iObjective) / objectiveScales(iObjective);
            }
        }
        if (constraintCount != 0) {
            referenceResponses.segment(objectiveCount, constraintCount) = problem.constraintValues;
            gradients.block(0, objectiveCount, variableCount, constraintCount) = *problem.constraintGradients;
        } else if (responseCount > objectiveCount) {
            referenceResponses(objectiveCount) = -1.0;
        }
        if (problem.objectiveCurvature.has_value()) {
            approximationFunction.ConfigureQuadraticModel(
                design,
                referenceResponses,
                gradients,
                curvature,
                objectiveMask,
                objectiveCount
            );
        } else {
            approximationFunction.ConfigureLinearModel(design, referenceResponses, gradients, objectiveMask, objectiveCount);
        }

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
        ExtractLaminateSectionPointersAndOffsets(sectionBindings, sectionPointers, sectionOffsets);

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
        result.predictedObjectives = predictedResponses.head(objectiveCount);
        for (Eigen::Index iObjective = 0; iObjective < objectiveCount; ++iObjective) {
            result.predictedObjectives(iObjective) *= objectiveScales(iObjective);
        }
        result.predictedConstraints = constraintCount == 0
            ? Eigen::VectorXd()
            : predictedResponses.segment(objectiveCount, constraintCount);
        result.message = "Candidate produced by the core laminate section adapter.";
        return result;
    }

private:
    CoreLaminateSection1RespOptions m_options;
};

}  // namespace lamopt
