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
#include <stdexcept>
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

struct CheckpointState {
    bool converged = false;
    AnalysisStatus status = AnalysisStatus::BackendFailure;
    Eigen::VectorXd design;
    std::size_t historyCount = 0;
    std::string message;
};

inline CheckpointState ReadCheckpoint(const std::filesystem::path& checkpointPath);

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

    GlobalOptimisationResult optimiseFromCheckpoint(const std::filesystem::path& checkpointPath,
                                                    const AnalysisRequest& requestTemplate) {
        const CheckpointState checkpoint = ReadCheckpoint(checkpointPath);
        AnalysisRequest resumedRequest = requestTemplate;
        resumedRequest.designVariables = checkpoint.design;

        GlobalOptimisationResult result = optimise(resumedRequest);
        if (checkpoint.converged) {
            appendDiagnosticMessage(result.message, "Checkpoint was already marked converged.");
        } else {
            appendDiagnosticMessage(result.message,
                                    "Resumed from checkpoint with " + std::to_string(checkpoint.historyCount)
                                        + " recorded iterations.");
        }
        return result;
    }

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

        if (baseResult.sensitivityPolicy.has_value() && !baseResult.hasConsistentSensitivityPolicy()) {
            appendDiagnosticMessage(baseResult.diagnostics.message,
                                    "Response sensitivity policy does not match the objective/constraint counts.");
        }

        attachOptimiserSideLaminateGradients(request, baseResult);
        if (baseResult.hasAllGradients()) {
            return baseResult;
        }

        if (m_options.gradientFallbackMode != GradientFallbackMode::FiniteDifference) {
            baseResult.status = AnalysisStatus::MissingGradients;
            if (baseResult.sensitivityPolicy.has_value()) {
                appendDiagnosticMessage(
                    baseResult.diagnostics.message,
                    "Configured backend sensitivity policy: "
                        + DescribeResponseSensitivityPolicy(*baseResult.sensitivityPolicy) + "."
                );
            }
            appendDiagnosticMessage(baseResult.diagnostics.message,
                                    "Analysis backend returned incomplete gradients and finite-difference fallback is disabled.");
            return baseResult;
        }

        const Eigen::Index nVar = request.designVariables.size();
        const Eigen::Index nObj = baseResult.objectives.size();
        const Eigen::Index nCon = baseResult.constraints.size();
        const bool needObjectiveGradients = !baseResult.hasObjectiveGradients();
        const bool needConstraintGradients = !baseResult.hasConstraintGradients();
        if (needObjectiveGradients) {
            ensureGradientStorage(baseResult.objectiveGradients,
                                  baseResult.objectiveGradientMask,
                                  nVar,
                                  nObj);
        }
        if (needConstraintGradients) {
            ensureGradientStorage(baseResult.constraintGradients,
                                  baseResult.constraintGradientMask,
                                  nVar,
                                  nCon);
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
                const Eigen::RowVectorXd objectiveRow =
                    ((perturbedResult.objectives - baseResult.objectives) / m_options.finiteDifferenceStep).transpose();
                fillMissingGradientRow(*baseResult.objectiveGradients,
                                       *baseResult.objectiveGradientMask,
                                       iVar,
                                       objectiveRow);
            }
            if (needConstraintGradients && nCon != 0) {
                const Eigen::RowVectorXd constraintRow =
                    ((perturbedResult.constraints - baseResult.constraints) / m_options.finiteDifferenceStep).transpose();
                fillMissingGradientRow(*baseResult.constraintGradients,
                                       *baseResult.constraintGradientMask,
                                       iVar,
                                       constraintRow);
            }
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

        const Eigen::Index nVar = request.designVariables.size();
        if (!mergeGradientContribution(result.objectiveGradients,
                                       result.objectiveGradientMask,
                                       derivativeResult.objectiveGradients,
                                       nVar,
                                       result.objectives.size(),
                                       derivativeResult.message)) {
            appendDiagnosticMessage(result.diagnostics.message, derivativeResult.message);
            return;
        }
        if (!mergeGradientContribution(result.constraintGradients,
                                       result.constraintGradientMask,
                                       derivativeResult.constraintGradients,
                                       nVar,
                                       result.constraints.size(),
                                       derivativeResult.message)) {
            appendDiagnosticMessage(result.diagnostics.message, derivativeResult.message);
            return;
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

    static void ensureGradientStorage(std::optional<Eigen::MatrixXd>& gradients,
                                      std::optional<Eigen::MatrixXi>& gradientMask,
                                      Eigen::Index nVar,
                                      Eigen::Index nResp) {
        const bool hadGradients = gradients.has_value();
        if (!gradients.has_value()) {
            gradients = Eigen::MatrixXd::Zero(nVar, nResp);
        }
        if (!gradientMask.has_value()) {
            if (hadGradients) {
                gradientMask = Eigen::MatrixXi::Ones(nVar, nResp);
            } else {
                gradientMask = Eigen::MatrixXi::Zero(nVar, nResp);
            }
        }
    }

    static void fillMissingGradientRow(Eigen::MatrixXd& gradients,
                                       Eigen::MatrixXi& gradientMask,
                                       Eigen::Index rowIndex,
                                       const Eigen::RowVectorXd& rowValues) {
        for (Eigen::Index iResp = 0; iResp < rowValues.size(); ++iResp) {
            if (gradientMask(rowIndex, iResp) == 0) {
                gradients(rowIndex, iResp) = rowValues(iResp);
                gradientMask(rowIndex, iResp) = 1;
            }
        }
    }

    static bool mergeGradientContribution(std::optional<Eigen::MatrixXd>& gradients,
                                          std::optional<Eigen::MatrixXi>& gradientMask,
                                          const std::optional<GradientContributionMatrix>& contribution,
                                          Eigen::Index nVar,
                                          Eigen::Index nResp,
                                          std::string& message) {
        if (!contribution.has_value()) {
            return true;
        }
        if (contribution->values.rows() != nVar || contribution->values.cols() != nResp
            || contribution->mask.rows() != nVar || contribution->mask.cols() != nResp) {
            message = "Optimiser-side laminate derivative contribution dimensions do not match the analysis gradients.";
            return false;
        }

        ensureGradientStorage(gradients, gradientMask, nVar, nResp);
        for (Eigen::Index iVar = 0; iVar < nVar; ++iVar) {
            for (Eigen::Index iResp = 0; iResp < nResp; ++iResp) {
                if (contribution->mask(iVar, iResp) != 0 && (*gradientMask)(iVar, iResp) == 0) {
                    (*gradients)(iVar, iResp) = contribution->values(iVar, iResp);
                    (*gradientMask)(iVar, iResp) = 1;
                }
            }
        }

        return true;
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

inline CheckpointState ReadCheckpoint(const std::filesystem::path& checkpointPath) {
    std::ifstream checkpoint(checkpointPath);
    if (!checkpoint) {
        throw std::runtime_error("Unable to open checkpoint file: " + checkpointPath.string());
    }

    CheckpointState state;
    bool sawConverged = false;
    bool sawStatus = false;
    bool sawDesign = false;
    bool sawHistoryCount = false;

    std::string line;
    while (std::getline(checkpoint, line)) {
        const std::size_t separator = line.find('=');
        if (separator == std::string::npos) {
            continue;
        }

        const std::string key = line.substr(0, separator);
        const std::string value = line.substr(separator + 1);

        if (key == "converged") {
            if (value == "true") {
                state.converged = true;
            } else if (value == "false") {
                state.converged = false;
            } else {
                throw std::runtime_error("Invalid checkpoint converged value: " + value);
            }
            sawConverged = true;
            continue;
        }

        if (key == "message") {
            state.message = value;
            continue;
        }

        if (key == "status") {
            const std::optional<AnalysisStatus> status = ParseAnalysisStatus(value);
            if (!status.has_value()) {
                throw std::runtime_error("Invalid checkpoint analysis status: " + value);
            }
            state.status = *status;
            sawStatus = true;
            continue;
        }

        if (key == "design") {
            if (value.empty()) {
                state.design = Eigen::VectorXd();
            } else {
                std::vector<double> entries;
                std::stringstream values(value);
                std::string token;
                while (std::getline(values, token, ',')) {
                    entries.push_back(std::stod(token));
                }
                state.design = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(entries.size()));
                for (std::size_t index = 0; index < entries.size(); ++index) {
                    state.design(static_cast<Eigen::Index>(index)) = entries[index];
                }
            }
            sawDesign = true;
            continue;
        }

        if (key == "history_count") {
            state.historyCount = static_cast<std::size_t>(std::stoull(value));
            sawHistoryCount = true;
        }
    }

    if (!sawConverged || !sawStatus || !sawDesign || !sawHistoryCount) {
        throw std::runtime_error("Checkpoint file is missing one or more required fields: "
                                 + checkpointPath.string());
    }

    return state;
}

}  // namespace lamopt
