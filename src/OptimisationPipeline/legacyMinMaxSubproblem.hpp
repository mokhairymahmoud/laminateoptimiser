#pragma once

#include "subproblem.hpp"
#include "../BoundSDP/boundSDP.hpp"
#include "../GlobalOptimiser/approxFunction.hpp"
#include "../GlobalOptimiser/scminmaxProb.hpp"

#include <algorithm>
#include <cmath>
#include <string>

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

}  // namespace lamopt
