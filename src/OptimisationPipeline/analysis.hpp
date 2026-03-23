#pragma once

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

struct AnalysisDiagnostics {
    std::string backendName;
    std::string message;
    std::string command;
    std::filesystem::path runDirectory;
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
    std::optional<Eigen::MatrixXd> objectiveGradients;
    std::optional<Eigen::MatrixXd> constraintGradients;
    std::optional<Eigen::MatrixXd> objectiveCurvature;
    AnalysisDiagnostics diagnostics;

    [[nodiscard]] bool isSuccessful() const {
        return status == AnalysisStatus::Success;
    }

    [[nodiscard]] bool hasObjectiveGradients() const {
        return objectives.size() == 0 || objectiveGradients.has_value();
    }

    [[nodiscard]] bool hasConstraintGradients() const {
        return constraints.size() == 0 || constraintGradients.has_value();
    }

    [[nodiscard]] bool hasAllGradients() const {
        return hasObjectiveGradients() && hasConstraintGradients();
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

}  // namespace lamopt
