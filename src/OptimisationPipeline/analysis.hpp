#pragma once

#include "responseSchema.hpp"

#include <Eigen/Dense>

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace lamopt {

enum class AnalysisStatus {
    Success,
    BackendFailure,
    Timeout,
    InvalidOutput,
    MissingGradients
};

enum class ResponseSensitivitySource {
    Unspecified,
    BackendNative,
    OptimiserSideLaminate,
    FiniteDifference
};

struct ResponseSensitivityPolicy {
    std::vector<ResponseSensitivitySource> objectiveSources;
    std::vector<ResponseSensitivitySource> constraintSources;
};

inline std::string ResponseSensitivitySourceToString(ResponseSensitivitySource source) {
    switch (source) {
        case ResponseSensitivitySource::Unspecified:
            return "unspecified";
        case ResponseSensitivitySource::BackendNative:
            return "backend_native";
        case ResponseSensitivitySource::OptimiserSideLaminate:
            return "optimiser_side_laminate";
        case ResponseSensitivitySource::FiniteDifference:
            return "finite_difference";
    }
    return "unknown";
}

inline std::string DescribeResponseSensitivityPolicy(const ResponseSensitivityPolicy& policy) {
    std::string description;

    for (std::size_t index = 0; index < policy.objectiveSources.size(); ++index) {
        if (!description.empty()) {
            description += "; ";
        }
        description += "objective[" + std::to_string(index) + "]=";
        description += ResponseSensitivitySourceToString(policy.objectiveSources[index]);
    }
    for (std::size_t index = 0; index < policy.constraintSources.size(); ++index) {
        if (!description.empty()) {
            description += "; ";
        }
        description += "constraint[" + std::to_string(index) + "]=";
        description += ResponseSensitivitySourceToString(policy.constraintSources[index]);
    }

    if (description.empty()) {
        return "no_response_sensitivity_policy";
    }
    return description;
}

struct AnalysisDiagnostics {
    std::string backendName;
    std::string message;
    std::string command;
    std::filesystem::path runDirectory;
    std::filesystem::path standardOutputPath;
    std::filesystem::path standardErrorPath;
    int exitCode = 0;
    int attempts = 0;
};

struct LaminateSectionState {
    Eigen::Index variableOffset = 0;
    double thickness = 0.0;
    double thicknessLowerBound = 0.0;
    double thicknessUpperBound = 1.0;
    Eigen::VectorXd laminationParameters;
    bool isBalanced = false;
    bool isSymmetric = true;
    int sublaminateCount = 0;
};

struct AnalysisRequest {
    Eigen::VectorXd designVariables;
    std::filesystem::path workDirectory;
    bool requestSensitivities = false;
    std::optional<std::string> runLabel;
    std::optional<Eigen::VectorXd> lowerBounds;
    std::optional<Eigen::VectorXd> upperBounds;
    std::vector<LaminateSectionState> laminateSections;
};

struct AnalysisResult {
    AnalysisStatus status = AnalysisStatus::BackendFailure;
    Eigen::VectorXd objectives;
    Eigen::VectorXd constraints;
    NamedScalarResponseMap extractedScalarValues;
    std::optional<ResponseSensitivityPolicy> sensitivityPolicy;
    std::optional<Eigen::MatrixXd> objectiveGradients;
    std::optional<Eigen::MatrixXi> objectiveGradientMask;
    std::optional<Eigen::MatrixXd> constraintGradients;
    std::optional<Eigen::MatrixXi> constraintGradientMask;
    std::optional<Eigen::MatrixXd> objectiveCurvature;
    AnalysisDiagnostics diagnostics;

    [[nodiscard]] bool isSuccessful() const {
        return status == AnalysisStatus::Success;
    }

    [[nodiscard]] bool hasObjectiveGradients() const {
        return hasCompleteGradientMatrix(objectives, objectiveGradients, objectiveGradientMask);
    }

    [[nodiscard]] bool hasConstraintGradients() const {
        return hasCompleteGradientMatrix(constraints, constraintGradients, constraintGradientMask);
    }

    [[nodiscard]] bool hasAllGradients() const {
        return hasObjectiveGradients() && hasConstraintGradients();
    }

    [[nodiscard]] bool hasConsistentSensitivityPolicy() const {
        return !sensitivityPolicy.has_value()
            || (sensitivityPolicy->objectiveSources.size() == static_cast<std::size_t>(objectives.size())
                && sensitivityPolicy->constraintSources.size() == static_cast<std::size_t>(constraints.size()));
    }

private:
    [[nodiscard]] static bool hasCompleteGradientMatrix(const Eigen::VectorXd& responses,
                                                        const std::optional<Eigen::MatrixXd>& gradients,
                                                        const std::optional<Eigen::MatrixXi>& gradientMask) {
        if (responses.size() == 0) {
            return true;
        }
        if (!gradients.has_value() || gradients->cols() != responses.size()) {
            return false;
        }
        if (!gradientMask.has_value()) {
            return true;
        }
        return gradientMask->rows() == gradients->rows()
            && gradientMask->cols() == gradients->cols()
            && gradientMask->minCoeff() != 0;
    }
};

class AnalysisBackend {
public:
    virtual ~AnalysisBackend() = default;
    virtual AnalysisResult evaluate(const AnalysisRequest& request) = 0;
};

inline double MaxObjectiveValue(const AnalysisResult& result) {
    return result.objectives.size() == 0 ? 0.0 : result.objectives.maxCoeff();
}

inline double MaxConstraintValue(const AnalysisResult& result) {
    return result.constraints.size() == 0 ? 0.0 : result.constraints.maxCoeff();
}

inline std::string AnalysisStatusToString(AnalysisStatus status) {
    switch (status) {
        case AnalysisStatus::Success:
            return "success";
        case AnalysisStatus::BackendFailure:
            return "backend_failure";
        case AnalysisStatus::Timeout:
            return "timeout";
        case AnalysisStatus::InvalidOutput:
            return "invalid_output";
        case AnalysisStatus::MissingGradients:
            return "missing_gradients";
    }
    return "unknown";
}

inline std::optional<AnalysisStatus> ParseAnalysisStatus(const std::string& status) {
    if (status == "success") {
        return AnalysisStatus::Success;
    }
    if (status == "backend_failure") {
        return AnalysisStatus::BackendFailure;
    }
    if (status == "timeout") {
        return AnalysisStatus::Timeout;
    }
    if (status == "invalid_output") {
        return AnalysisStatus::InvalidOutput;
    }
    if (status == "missing_gradients") {
        return AnalysisStatus::MissingGradients;
    }
    return std::nullopt;
}

}  // namespace lamopt
