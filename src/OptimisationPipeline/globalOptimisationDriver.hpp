#pragma once

#include "analysis.hpp"
#include "approximation.hpp"
#include "laminationParameterDerivatives.hpp"
#include "subproblem.hpp"
#include "../GlobalOptimiser/dampingUtils.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace lamopt {

enum class GradientFallbackMode {
    Disabled,
    FiniteDifference
};

struct DriverOptions {
    int maxOuterIterations = 25;
    int maxSubIterations = 8;
    double objectiveTolerance = 1.0e-6;
    double constraintTolerance = 1.0e-6;
    double stagnationTolerance = 1.0e-8;
    double finiteDifferenceStep = 1.0e-6;
    double dampingFloor = 1.0e-8;
    bool requestSensitivities = true;
    GradientFallbackMode gradientFallbackMode = GradientFallbackMode::Disabled;
};

struct IterationRecord {
    int outerIteration = 0;
    int subIteration = 0;
    bool accepted = false;
    double referenceObjective = 0.0;
    double candidateObjective = 0.0;
    double maxConstraint = 0.0;
    Eigen::VectorXd design;
    Eigen::VectorXd dampingFactorsBefore;
    Eigen::VectorXd dampingFactorsAfter;
    std::string message;
};

struct GlobalOptimisationResult {
    bool converged = false;
    Eigen::VectorXd design;
    AnalysisResult analysis;
    std::vector<IterationRecord> history;
    std::string message;
};

class GlobalOptimisationDriver {
public:
    GlobalOptimisationDriver(AnalysisBackend& analysisBackend,
                             ApproximationBuilder& approximationBuilder,
                             SubproblemSolver& subproblemSolver,
                             DriverOptions options = {},
                             LaminationParameterDerivativeProvider* laminationDerivativeProvider = nullptr)
        : m_analysisBackend(analysisBackend)
        , m_approximationBuilder(approximationBuilder)
        , m_subproblemSolver(subproblemSolver)
        , m_options(options)
        , m_laminationDerivativeProvider(laminationDerivativeProvider) {}

    GlobalOptimisationResult optimise(const AnalysisRequest& initialRequest) {
        GlobalOptimisationResult result;

        AnalysisRequest referenceRequest = initialRequest;
        referenceRequest.requestSensitivities = m_options.requestSensitivities;

        AnalysisResult referenceResult = evaluateWithFallback(referenceRequest);
        if (!referenceResult.isSuccessful()) {
            result.design = referenceRequest.designVariables;
            result.analysis = referenceResult;
            result.message = "Initial analysis failed: " + referenceResult.diagnostics.message;
            return result;
        }

        const Eigen::Index totalResponses = referenceResult.objectives.size() + referenceResult.constraints.size();
        Eigen::VectorXd dampingFactors = Eigen::VectorXd::Ones(totalResponses);
        Eigen::VectorXd designDampingVector = Eigen::VectorXd::Ones(referenceRequest.designVariables.size());

        for (int outerIteration = 0; outerIteration < m_options.maxOuterIterations; ++outerIteration) {
            bool acceptedCandidate = false;
            const double referenceObjective = MaxObjectiveValue(referenceResult);

            for (int subIteration = 0; subIteration < m_options.maxSubIterations; ++subIteration) {
                ApproximationProblem approximation = m_approximationBuilder.build(
                    referenceRequest,
                    referenceResult,
                    designDampingVector,
                    dampingFactors
                );

                SubproblemResult subproblemResult = m_subproblemSolver.solve(approximation);
                IterationRecord record;
                record.outerIteration = outerIteration;
                record.subIteration = subIteration;
                record.referenceObjective = referenceObjective;
                record.dampingFactorsBefore = dampingFactors;
                record.message = subproblemResult.message;

                if (!subproblemResult.success) {
                    record.dampingFactorsAfter = dampingFactors;
                    result.history.push_back(record);
                    result.design = referenceRequest.designVariables;
                    result.analysis = referenceResult;
                    result.message = "Subproblem solve failed: " + subproblemResult.message;
                    return result;
                }

                AnalysisRequest candidateRequest = referenceRequest;
                candidateRequest.designVariables = subproblemResult.candidateDesign;
                candidateRequest.requestSensitivities = m_options.requestSensitivities;
                candidateRequest.workDirectory = makeSubdirectory(referenceRequest.workDirectory, outerIteration, subIteration);

                AnalysisResult candidateResult = evaluateWithFallback(candidateRequest);
                record.design = subproblemResult.candidateDesign;
                record.candidateObjective = MaxObjectiveValue(candidateResult);
                record.maxConstraint = MaxConstraintValue(candidateResult);

                if (!candidateResult.isSuccessful()) {
                    record.dampingFactorsAfter = dampingFactors;
                    record.message += " Analysis failed: " + candidateResult.diagnostics.message;
                    result.history.push_back(record);
                    result.design = referenceRequest.designVariables;
                    result.analysis = candidateResult;
                    result.message = "Candidate analysis failed.";
                    return result;
                }

                if (isAccepted(referenceResult, candidateResult)) {
                    record.accepted = true;
                    record.dampingFactorsAfter = dampingFactors;
                    result.history.push_back(record);

                    const Eigen::VectorXd oldReferenceDesign = referenceRequest.designVariables;
                    referenceRequest = candidateRequest;
                    referenceResult = candidateResult;
                    designDampingVector = (referenceRequest.designVariables - oldReferenceDesign).cwiseAbs();
                    designDampingVector = designDampingVector.unaryExpr([this](double value) {
                        return std::max(value, m_options.dampingFloor);
                    });

                    const double objectiveImprovement = std::abs(record.referenceObjective - record.candidateObjective);
                    const double designShift = (referenceRequest.designVariables - oldReferenceDesign).norm();
                    if (objectiveImprovement <= m_options.stagnationTolerance
                        || designShift <= m_options.stagnationTolerance) {
                        result.converged = true;
                        result.design = referenceRequest.designVariables;
                        result.analysis = referenceResult;
                        result.message = "Converged after accepted improvement step.";
                        return result;
                    }

                    acceptedCandidate = true;
                    break;
                }

                updateDampingFactors(candidateResult,
                                     subproblemResult,
                                     candidateRequest.designVariables,
                                     referenceRequest.designVariables,
                                     dampingFactors);
                record.dampingFactorsAfter = dampingFactors;
                result.history.push_back(record);
            }

            if (!acceptedCandidate) {
                result.design = referenceRequest.designVariables;
                result.analysis = referenceResult;
                result.message = "Reached the sub-iteration limit without an accepted candidate.";
                return result;
            }
        }

        result.design = referenceRequest.designVariables;
        result.analysis = referenceResult;
        result.message = "Reached the outer-iteration limit without convergence.";
        return result;
    }

private:
    AnalysisBackend& m_analysisBackend;
    ApproximationBuilder& m_approximationBuilder;
    SubproblemSolver& m_subproblemSolver;
    DriverOptions m_options;
    LaminationParameterDerivativeProvider* m_laminationDerivativeProvider = nullptr;

    AnalysisResult evaluateWithFallback(const AnalysisRequest& request) {
        AnalysisResult baseResult = m_analysisBackend.evaluate(request);
        if (!baseResult.isSuccessful()) {
            return baseResult;
        }

        if (!request.requestSensitivities || baseResult.hasAllGradients()) {
            return baseResult;
        }

        attachOptimiserSideLaminateGradients(request, baseResult);
        if (baseResult.hasAllGradients()) {
            return baseResult;
        }

        if (m_options.gradientFallbackMode != GradientFallbackMode::FiniteDifference) {
            baseResult.status = AnalysisStatus::MissingGradients;
            appendDiagnosticMessage(baseResult.diagnostics.message,
                                    "Analysis backend returned incomplete gradients and finite-difference fallback is disabled.");
            return baseResult;
        }

        const Eigen::Index nVar = request.designVariables.size();
        const Eigen::Index nObj = baseResult.objectives.size();
        const Eigen::Index nCon = baseResult.constraints.size();
        const bool needObjectiveGradients = !baseResult.hasObjectiveGradients();
        const bool needConstraintGradients = !baseResult.hasConstraintGradients();
        Eigen::MatrixXd objectiveGradients =
            Eigen::MatrixXd::Zero(nVar, nObj);
        Eigen::MatrixXd constraintGradients =
            Eigen::MatrixXd::Zero(nVar, nCon);
        if (!needObjectiveGradients && baseResult.objectiveGradients.has_value()) {
            objectiveGradients = *baseResult.objectiveGradients;
        }
        if (!needConstraintGradients && baseResult.constraintGradients.has_value()) {
            constraintGradients = *baseResult.constraintGradients;
        }

        for (Eigen::Index iVar = 0; iVar < nVar; ++iVar) {
            AnalysisRequest perturbedRequest = request;
            perturbedRequest.requestSensitivities = false;
            perturbedRequest.designVariables(iVar) += m_options.finiteDifferenceStep;
            perturbedRequest.workDirectory = makeFiniteDifferenceDirectory(request.workDirectory, iVar);

            AnalysisResult perturbedResult = m_analysisBackend.evaluate(perturbedRequest);
            if (!perturbedResult.isSuccessful()) {
                perturbedResult.status = AnalysisStatus::MissingGradients;
                perturbedResult.diagnostics.message =
                    "Finite-difference gradient evaluation failed at variable " + std::to_string(iVar) + ".";
                return perturbedResult;
            }

            if (needObjectiveGradients && nObj != 0) {
                objectiveGradients.row(iVar) =
                    ((perturbedResult.objectives - baseResult.objectives) / m_options.finiteDifferenceStep).transpose();
            }
            if (needConstraintGradients && nCon != 0) {
                constraintGradients.row(iVar) =
                    ((perturbedResult.constraints - baseResult.constraints) / m_options.finiteDifferenceStep).transpose();
            }
        }

        if (needObjectiveGradients) {
            baseResult.objectiveGradients = objectiveGradients;
        }
        if (needConstraintGradients) {
            baseResult.constraintGradients = constraintGradients;
        }
        baseResult.status = AnalysisStatus::Success;
        appendDiagnosticMessage(baseResult.diagnostics.message, "Finite-difference gradients attached.");
        return baseResult;
    }

    [[nodiscard]] bool isAccepted(const AnalysisResult& referenceResult,
                                  const AnalysisResult& candidateResult) const {
        return MaxObjectiveValue(candidateResult) <= MaxObjectiveValue(referenceResult) + m_options.objectiveTolerance
            && MaxConstraintValue(candidateResult) <= m_options.constraintTolerance;
    }

    void updateDampingFactors(const AnalysisResult& referenceResult,
                              const SubproblemResult& subproblemResult,
                              const Eigen::VectorXd& candidateDesign,
                              const Eigen::VectorXd& referenceDesign,
                              Eigen::VectorXd& dampingFactors) const {
        const Eigen::VectorXd actual = concatenate(referenceResult.objectives, referenceResult.constraints);
        const Eigen::VectorXd predicted = concatenate(subproblemResult.predictedObjectives,
                                                      subproblemResult.predictedConstraints);

        const double movement = std::max((candidateDesign - referenceDesign).cwiseAbs().maxCoeff(),
                                         m_options.dampingFloor);

        if (actual.size() == predicted.size() && actual.size() == dampingFactors.size()) {
            lampar::GlobalOptimiser::UpdateDamping(actual, predicted, movement, dampingFactors);
        } else {
            dampingFactors *= 2.0;
        }
    }

    [[nodiscard]] static Eigen::VectorXd concatenate(const Eigen::VectorXd& lhs,
                                                     const Eigen::VectorXd& rhs) {
        Eigen::VectorXd combined(lhs.size() + rhs.size());
        if (lhs.size() != 0) {
            combined.head(lhs.size()) = lhs;
        }
        if (rhs.size() != 0) {
            combined.tail(rhs.size()) = rhs;
        }
        return combined;
    }

    [[nodiscard]] static std::filesystem::path makeSubdirectory(const std::filesystem::path& root,
                                                                int outerIteration,
                                                                int subIteration) {
        if (root.empty()) {
            return {};
        }
        return root / ("outer_" + std::to_string(outerIteration) + "_sub_" + std::to_string(subIteration));
    }

    [[nodiscard]] static std::filesystem::path makeFiniteDifferenceDirectory(const std::filesystem::path& root,
                                                                             Eigen::Index variableIndex) {
        if (root.empty()) {
            return {};
        }
        return root / ("fd_" + std::to_string(variableIndex));
    }

    void attachOptimiserSideLaminateGradients(const AnalysisRequest& request,
                                              AnalysisResult& result) const {
        if (m_laminationDerivativeProvider == nullptr || result.hasAllGradients()) {
            return;
        }

        LaminationParameterDerivativeResult derivativeResult =
            m_laminationDerivativeProvider->compute(request, result);
        if (!derivativeResult.success) {
            if (!derivativeResult.message.empty()) {
                appendDiagnosticMessage(result.diagnostics.message, derivativeResult.message);
            }
            return;
        }

        if (!result.hasObjectiveGradients() && derivativeResult.objectiveGradients.has_value()) {
            result.objectiveGradients = derivativeResult.objectiveGradients;
        }
        if (!result.hasConstraintGradients() && derivativeResult.constraintGradients.has_value()) {
            result.constraintGradients = derivativeResult.constraintGradients;
        }
        if (!derivativeResult.message.empty()) {
            appendDiagnosticMessage(result.diagnostics.message, derivativeResult.message);
        }
    }

    static void appendDiagnosticMessage(std::string& message, const std::string& suffix) {
        if (suffix.empty()) {
            return;
        }
        if (!message.empty() && message.back() != ' ') {
            message += ' ';
        }
        message += suffix;
    }
};

inline void WriteCheckpoint(const std::filesystem::path& checkpointPath,
                            const GlobalOptimisationResult& result) {
    std::ofstream checkpoint(checkpointPath);
    checkpoint << "converged=" << (result.converged ? "true" : "false") << '\n';
    checkpoint << "message=" << result.message << '\n';
    checkpoint << "status=" << AnalysisStatusToString(result.analysis.status) << '\n';
    checkpoint << "design=";
    for (Eigen::Index i = 0; i < result.design.size(); ++i) {
        if (i != 0) {
            checkpoint << ',';
        }
        checkpoint << result.design(i);
    }
    checkpoint << '\n';
    checkpoint << "history_count=" << result.history.size() << '\n';
}

}  // namespace lamopt
