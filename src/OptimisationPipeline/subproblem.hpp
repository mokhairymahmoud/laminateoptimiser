#pragma once

#include "approximation.hpp"

#include <cmath>
#include <string>

namespace lamopt {

struct SubproblemResult {
    bool success = false;
    Eigen::VectorXd candidateDesign;
    Eigen::VectorXd predictedObjectives;
    Eigen::VectorXd predictedConstraints;
    int iterations = 0;
    std::string message;
};

class SubproblemSolver {
public:
    virtual ~SubproblemSolver() = default;
    virtual SubproblemResult solve(const ApproximationProblem& problem) = 0;
};

struct GradientPenaltySubproblemOptions {
    double initialStep = 0.25;
    double constraintPenalty = 5.0;
    double minimumDirectionNorm = 1.0e-12;
};

class GradientPenaltySubproblemSolver : public SubproblemSolver {
public:
    explicit GradientPenaltySubproblemSolver(GradientPenaltySubproblemOptions options = {})
        : m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        SubproblemResult result;
        result.iterations = 1;

        if (problem.objectiveValues.size() == 0) {
            result.message = "At least one objective value is required.";
            return result;
        }

        if (!problem.objectiveGradients.has_value()) {
            result.message = "Objective gradients are required for the gradient-penalty solver.";
            return result;
        }

        const Eigen::MatrixXd& objectiveGradients = *problem.objectiveGradients;
        if (objectiveGradients.rows() != problem.referenceDesign.size()) {
            result.message = "Objective gradient dimensions do not match the design vector.";
            return result;
        }

        problem.objectiveValues.maxCoeff(&m_activeObjectiveIndex);

        Eigen::VectorXd direction = -objectiveGradients.col(m_activeObjectiveIndex);
        if (problem.constraintGradients.has_value()) {
            const Eigen::MatrixXd& constraintGradients = *problem.constraintGradients;
            if (constraintGradients.rows() != problem.referenceDesign.size()) {
                result.message = "Constraint gradient dimensions do not match the design vector.";
                return result;
            }
            for (Eigen::Index i = 0; i < problem.constraintValues.size(); ++i) {
                if (problem.constraintValues(i) > 0.0) {
                    direction -= m_options.constraintPenalty * (1.0 + problem.constraintValues(i))
                               * constraintGradients.col(i);
                }
            }
        } else if (problem.constraintValues.size() != 0 && problem.constraintValues.maxCoeff() > 0.0) {
            result.message = "Constraint gradients are required when violated constraints are present.";
            return result;
        }

        const double directionNorm = direction.norm();
        if (directionNorm < m_options.minimumDirectionNorm) {
            result.success = true;
            result.candidateDesign = problem.referenceDesign;
            result.predictedObjectives = problem.objectiveValues;
            result.predictedConstraints = problem.constraintValues;
            result.message = "Zero search direction.";
            return result;
        }

        const double dampingScale = problem.responseDampingFactors.size() == 0
            ? 1.0
            : 1.0 + problem.responseDampingFactors.mean();
        const double step = m_options.initialStep / dampingScale;
        Eigen::VectorXd delta = step * direction / directionNorm;
        result.candidateDesign = problem.referenceDesign + delta;

        if (problem.lowerBounds.has_value()) {
            result.candidateDesign = result.candidateDesign.cwiseMax(*problem.lowerBounds);
        }
        if (problem.upperBounds.has_value()) {
            result.candidateDesign = result.candidateDesign.cwiseMin(*problem.upperBounds);
        }

        const Eigen::VectorXd realisedDelta = result.candidateDesign - problem.referenceDesign;
        result.predictedObjectives = problem.objectiveValues;
        if (problem.objectiveValues.size() != 0) {
            result.predictedObjectives += objectiveGradients.transpose() * realisedDelta;
        }

        result.predictedConstraints = problem.constraintValues;
        if (problem.constraintGradients.has_value() && problem.constraintValues.size() != 0) {
            result.predictedConstraints += problem.constraintGradients->transpose() * realisedDelta;
        }

        result.success = true;
        result.message = "Candidate produced by the gradient-penalty solver.";
        return result;
    }

private:
    GradientPenaltySubproblemOptions m_options;
    Eigen::Index m_activeObjectiveIndex = 0;
};

}  // namespace lamopt
