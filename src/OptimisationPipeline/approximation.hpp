#pragma once

#include "analysis.hpp"

#include <optional>
#include <vector>

namespace lamopt {

struct ApproximationProblem {
    Eigen::VectorXd referenceDesign;
    Eigen::VectorXd objectiveValues;
    Eigen::VectorXd constraintValues;
    std::optional<Eigen::MatrixXd> objectiveGradients;
    std::optional<Eigen::MatrixXd> constraintGradients;
    std::optional<Eigen::MatrixXd> objectiveCurvature;
    Eigen::VectorXd responseDampingFactors;
    Eigen::VectorXd designDampingVector;
    std::optional<Eigen::VectorXd> lowerBounds;
    std::optional<Eigen::VectorXd> upperBounds;
    std::vector<LaminateSectionState> laminateSections;
};

class ApproximationBuilder {
public:
    virtual ~ApproximationBuilder() = default;

    virtual ApproximationProblem build(const AnalysisRequest& request,
                                       const AnalysisResult& referenceResult,
                                       const Eigen::VectorXd& designDampingVector,
                                       const Eigen::VectorXd& responseDampingFactors) = 0;
};

class LinearApproximationBuilder : public ApproximationBuilder {
public:
    ApproximationProblem build(const AnalysisRequest& request,
                               const AnalysisResult& referenceResult,
                               const Eigen::VectorXd& designDampingVector,
                               const Eigen::VectorXd& responseDampingFactors) override {
        ApproximationProblem problem;
        problem.referenceDesign = request.designVariables;
        problem.objectiveValues = referenceResult.objectives;
        problem.constraintValues = referenceResult.constraints;
        problem.objectiveGradients = referenceResult.objectiveGradients;
        problem.constraintGradients = referenceResult.constraintGradients;
        problem.objectiveCurvature = referenceResult.objectiveCurvature;
        problem.designDampingVector = designDampingVector;
        problem.responseDampingFactors = responseDampingFactors;
        problem.lowerBounds = request.lowerBounds;
        problem.upperBounds = request.upperBounds;
        problem.laminateSections = request.laminateSections;
        return problem;
    }
};

}  // namespace lamopt
