#pragma once

#include "analysis.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace lamopt {

struct GradientContributionMatrix {
    Eigen::MatrixXd values;
    Eigen::MatrixXi mask;
};

struct LaminationParameterDerivativeResult {
    bool success = false;
    std::optional<GradientContributionMatrix> objectiveGradients;
    std::optional<GradientContributionMatrix> constraintGradients;
    std::string message;
};

class LaminationParameterDerivativeProvider {
public:
    virtual ~LaminationParameterDerivativeProvider() = default;

    [[nodiscard]] virtual bool supports(const AnalysisRequest& request,
                                        const AnalysisResult& result) const = 0;

    [[nodiscard]] virtual LaminationParameterDerivativeResult compute(const AnalysisRequest& request,
                                                                      const AnalysisResult& result) const = 0;
};

inline Eigen::Index LaminateSectionVariableCount(const LaminateSectionState& sectionState) {
    const Eigen::Index laminationParameterCount = sectionState.isBalanced ? 2 : 4;
    const Eigen::Index variableTypes = sectionState.isSymmetric ? 2 : 3;
    return laminationParameterCount * variableTypes + 1;
}

struct DenseDesignResponseTerm {
    Eigen::VectorXd linearCoefficients;
    Eigen::VectorXd quadraticCoefficients;
};

struct SeparableLaminateSectionResponseTerm {
    Eigen::VectorXd linearCoefficients;
    Eigen::VectorXd quadraticCoefficients;
};

struct SeparableLaminateResponse {
    double constant = 0.0;
    DenseDesignResponseTerm designTerms;
    std::vector<SeparableLaminateSectionResponseTerm> sectionTerms;
};

struct NamedSeparableLaminateResponse {
    std::string id;
    SeparableLaminateResponse response;
};

class SeparableLaminateResponseModel {
public:
    std::vector<SeparableLaminateResponse> objectives;
    std::vector<SeparableLaminateResponse> constraints;

    [[nodiscard]] bool validate(const AnalysisRequest& request,
                                const AnalysisResult& result,
                                std::string& message) const {
        if (request.laminateSections.empty()) {
            message = "At least one laminate section is required for optimiser-side laminate derivatives.";
            return false;
        }

        std::vector<Eigen::Index> sectionSizes;
        if (!validateSectionSlices(request, sectionSizes, message)) {
            return false;
        }

        if (result.objectives.size() != static_cast<Eigen::Index>(objectives.size())) {
            message = "Objective count does not match the separable laminate response model.";
            return false;
        }
        if (result.constraints.size() != static_cast<Eigen::Index>(constraints.size())) {
            message = "Constraint count does not match the separable laminate response model.";
            return false;
        }

        if (!validateResponses(objectives, request.designVariables.size(), sectionSizes, "objective", message)) {
            return false;
        }
        if (!validateResponses(constraints, request.designVariables.size(), sectionSizes, "constraint", message)) {
            return false;
        }

        return true;
    }

    [[nodiscard]] Eigen::VectorXd evaluateObjectives(const AnalysisRequest& request) const {
        return evaluateResponses(objectives, request);
    }

    [[nodiscard]] Eigen::VectorXd evaluateConstraints(const AnalysisRequest& request) const {
        return evaluateResponses(constraints, request);
    }

    [[nodiscard]] Eigen::MatrixXd objectiveDesignGradients(const AnalysisRequest& request) const {
        return evaluateDesignResponseGradients(objectives, request);
    }

    [[nodiscard]] Eigen::MatrixXd constraintDesignGradients(const AnalysisRequest& request) const {
        return evaluateDesignResponseGradients(constraints, request);
    }

    [[nodiscard]] Eigen::MatrixXd objectiveSectionGradients(const AnalysisRequest& request) const {
        return evaluateSectionResponseGradients(objectives, request);
    }

    [[nodiscard]] Eigen::MatrixXd constraintSectionGradients(const AnalysisRequest& request) const {
        return evaluateSectionResponseGradients(constraints, request);
    }

    [[nodiscard]] Eigen::MatrixXd objectiveGradients(const AnalysisRequest& request) const {
        return objectiveDesignGradients(request) + objectiveSectionGradients(request);
    }

    [[nodiscard]] Eigen::MatrixXd constraintGradients(const AnalysisRequest& request) const {
        return constraintDesignGradients(request) + constraintSectionGradients(request);
    }

public:
    [[nodiscard]] static bool validateSectionSlices(const AnalysisRequest& request,
                                                    std::vector<Eigen::Index>& sectionSizes,
                                                    std::string& message) {
        sectionSizes.clear();
        sectionSizes.reserve(request.laminateSections.size());

        std::vector<int> variableOwners(static_cast<size_t>(request.designVariables.size()), -1);

        for (size_t iSection = 0; iSection < request.laminateSections.size(); ++iSection) {
            const LaminateSectionState& sectionState = request.laminateSections[iSection];
            if (sectionState.sublaminateCount <= 0) {
                message = "Laminate section state must define a positive sublaminate count.";
                return false;
            }

            const Eigen::Index sectionSize = LaminateSectionVariableCount(sectionState);
            if (sectionState.laminationParameters.size() != 0
                && sectionState.laminationParameters.size() != sectionSize - 1) {
                message =
                    "Laminate section lamination-parameter size does not match the selected section type.";
                return false;
            }
            if (sectionState.variableOffset < 0
                || sectionState.variableOffset + sectionSize > request.designVariables.size()) {
                message = "Laminate section slice exceeds the reference design size.";
                return false;
            }

            for (Eigen::Index iVar = 0; iVar < sectionSize; ++iVar) {
                const Eigen::Index variableIndex = sectionState.variableOffset + iVar;
                if (variableOwners[static_cast<size_t>(variableIndex)] != -1) {
                    message = "Laminate section blocks must not overlap in the design vector.";
                    return false;
                }
                variableOwners[static_cast<size_t>(variableIndex)] = static_cast<int>(iSection);
            }

            sectionSizes.push_back(sectionSize);
        }

        return true;
    }

    [[nodiscard]] static bool validateResponses(const std::vector<SeparableLaminateResponse>& responses,
                                                const Eigen::Index designVariableCount,
                                                const std::vector<Eigen::Index>& sectionSizes,
                                                const char* responseKind,
                                                std::string& message) {
        for (size_t iResponse = 0; iResponse < responses.size(); ++iResponse) {
            const SeparableLaminateResponse& response = responses[iResponse];
            if (!validateDenseTerm(response.designTerms, designVariableCount, responseKind, message)) {
                return false;
            }
            if (!response.sectionTerms.empty() && response.sectionTerms.size() != sectionSizes.size()) {
                message = std::string("Each ") + responseKind
                        + " model must define zero section terms or one term per laminate section.";
                return false;
            }

            for (size_t iSection = 0; iSection < response.sectionTerms.size(); ++iSection) {
                const SeparableLaminateSectionResponseTerm& sectionTerm = response.sectionTerms[iSection];
                if (sectionTerm.linearCoefficients.size() != sectionSizes[iSection]
                    || sectionTerm.quadraticCoefficients.size() != sectionSizes[iSection]) {
                    message = std::string("A ") + responseKind
                            + " model term does not match the laminate section variable count.";
                    return false;
                }
            }
        }

        return true;
    }

    [[nodiscard]] static bool validateDenseTerm(const DenseDesignResponseTerm& designTerms,
                                                const Eigen::Index designVariableCount,
                                                const char* responseKind,
                                                std::string& message) {
        const bool linearValid =
            designTerms.linearCoefficients.size() == 0
            || designTerms.linearCoefficients.size() == designVariableCount;
        const bool quadraticValid =
            designTerms.quadraticCoefficients.size() == 0
            || designTerms.quadraticCoefficients.size() == designVariableCount;

        if (linearValid && quadraticValid) {
            return true;
        }

        message = std::string("The dense design term in a ") + responseKind
                + " model does not match the design-variable count.";
        return false;
    }

    [[nodiscard]] static Eigen::VectorXd evaluateResponses(const std::vector<SeparableLaminateResponse>& responses,
                                                           const AnalysisRequest& request) {
        Eigen::VectorXd values = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(responses.size()));

        for (size_t iResponse = 0; iResponse < responses.size(); ++iResponse) {
            values(static_cast<Eigen::Index>(iResponse)) =
                evaluateResponse(responses[iResponse], request);
        }

        return values;
    }

    [[nodiscard]] static double evaluateResponse(const SeparableLaminateResponse& response,
                                                 const AnalysisRequest& request) {
        double value = response.constant;
        value += evaluateDenseTerm(response.designTerms, request.designVariables);

        for (size_t iSection = 0; iSection < response.sectionTerms.size(); ++iSection) {
            const LaminateSectionState& sectionState = request.laminateSections[iSection];
            const SeparableLaminateSectionResponseTerm& sectionTerm = response.sectionTerms[iSection];
            const Eigen::Index sectionSize = sectionTerm.linearCoefficients.size();
            const Eigen::VectorXd sectionVariables =
                request.designVariables.segment(sectionState.variableOffset, sectionSize);

            value += sectionTerm.linearCoefficients.dot(sectionVariables);
            value += 0.5
                  * sectionTerm.quadraticCoefficients.dot(sectionVariables.cwiseProduct(sectionVariables));
        }

        return value;
    }

    [[nodiscard]] static double evaluateDenseTerm(const DenseDesignResponseTerm& designTerms,
                                                  const Eigen::VectorXd& designVariables) {
        double value = 0.0;

        if (designTerms.linearCoefficients.size() != 0) {
            value += designTerms.linearCoefficients.dot(designVariables);
        }
        if (designTerms.quadraticCoefficients.size() != 0) {
            value += 0.5
                  * designTerms.quadraticCoefficients.dot(designVariables.cwiseProduct(designVariables));
        }

        return value;
    }

    [[nodiscard]] static Eigen::MatrixXd evaluateDesignResponseGradients(
        const std::vector<SeparableLaminateResponse>& responses,
        const AnalysisRequest& request) {
        Eigen::MatrixXd gradients =
            Eigen::MatrixXd::Zero(request.designVariables.size(),
                                  static_cast<Eigen::Index>(responses.size()));

        for (size_t iResponse = 0; iResponse < responses.size(); ++iResponse) {
            gradients.col(static_cast<Eigen::Index>(iResponse)) =
                evaluateDenseTermGradient(responses[iResponse].designTerms, request.designVariables);
        }

        return gradients;
    }

    [[nodiscard]] static Eigen::VectorXd evaluateDenseTermGradient(const DenseDesignResponseTerm& designTerms,
                                                                   const Eigen::VectorXd& designVariables) {
        Eigen::VectorXd gradient = Eigen::VectorXd::Zero(designVariables.size());

        if (designTerms.linearCoefficients.size() != 0) {
            gradient += designTerms.linearCoefficients;
        }
        if (designTerms.quadraticCoefficients.size() != 0) {
            gradient += designTerms.quadraticCoefficients.cwiseProduct(designVariables);
        }

        return gradient;
    }

    [[nodiscard]] static Eigen::MatrixXd evaluateSectionResponseGradients(
        const std::vector<SeparableLaminateResponse>& responses,
        const AnalysisRequest& request) {
        Eigen::MatrixXd gradients =
            Eigen::MatrixXd::Zero(request.designVariables.size(),
                                  static_cast<Eigen::Index>(responses.size()));

        for (size_t iResponse = 0; iResponse < responses.size(); ++iResponse) {
            for (size_t iSection = 0; iSection < responses[iResponse].sectionTerms.size(); ++iSection) {
                const LaminateSectionState& sectionState = request.laminateSections[iSection];
                const SeparableLaminateSectionResponseTerm& sectionTerm =
                    responses[iResponse].sectionTerms[iSection];
                const Eigen::Index sectionSize = sectionTerm.linearCoefficients.size();
                const Eigen::VectorXd sectionVariables =
                    request.designVariables.segment(sectionState.variableOffset, sectionSize);

                gradients.block(sectionState.variableOffset,
                                static_cast<Eigen::Index>(iResponse),
                                sectionSize,
                                1) += sectionTerm.linearCoefficients
                                    + sectionTerm.quadraticCoefficients.cwiseProduct(sectionVariables);
            }
        }

        return gradients;
    }
};

class SeparableLaminateQuantityModel {
public:
    std::vector<NamedSeparableLaminateResponse> quantities;

    [[nodiscard]] bool validate(const AnalysisRequest& request,
                                const AnalysisResult& result,
                                std::string& message) const {
        if (request.laminateSections.empty()) {
            message = "At least one laminate section is required for optimiser-side laminate derivatives.";
            return false;
        }

        std::vector<Eigen::Index> sectionSizes;
        if (!SeparableLaminateResponseModel::validateSectionSlices(request, sectionSizes, message)) {
            return false;
        }

        std::unordered_map<std::string, std::size_t> ids;
        for (const NamedSeparableLaminateResponse& quantity : quantities) {
            if (quantity.id.empty()) {
                message = "Laminate quantity ids must not be empty.";
                return false;
            }
            if (!ids.emplace(quantity.id, ids.size()).second) {
                message = "Laminate quantity ids must be unique.";
                return false;
            }
            const std::vector<SeparableLaminateResponse> response = {quantity.response};
            if (!SeparableLaminateResponseModel::validateResponses(response,
                                                                  request.designVariables.size(),
                                                                  sectionSizes,
                                                                  "quantity",
                                                                  message)) {
                return false;
            }
            if (result.extractedScalarValues.find(quantity.id) == result.extractedScalarValues.end()) {
                message = "Extracted FE quantity not found for laminate derivative assembly: " + quantity.id;
                return false;
            }
        }

        return true;
    }

    [[nodiscard]] std::optional<std::size_t> findIndex(const std::string& id) const {
        for (std::size_t index = 0; index < quantities.size(); ++index) {
            if (quantities[index].id == id) {
                return index;
            }
        }
        return std::nullopt;
    }

    [[nodiscard]] Eigen::VectorXd evaluateQuantities(const AnalysisRequest& request) const {
        Eigen::VectorXd values = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(quantities.size()));
        for (std::size_t index = 0; index < quantities.size(); ++index) {
            values(static_cast<Eigen::Index>(index)) =
                SeparableLaminateResponseModel::evaluateResponse(quantities[index].response, request);
        }
        return values;
    }

    [[nodiscard]] Eigen::MatrixXd quantityGradients(const AnalysisRequest& request) const {
        Eigen::MatrixXd gradients =
            Eigen::MatrixXd::Zero(request.designVariables.size(), static_cast<Eigen::Index>(quantities.size()));
        for (std::size_t index = 0; index < quantities.size(); ++index) {
            gradients.col(static_cast<Eigen::Index>(index)) =
                responseGradient(quantities[index].response, request);
        }
        return gradients;
    }

    [[nodiscard]] Eigen::MatrixXd quantitySectionGradients(const AnalysisRequest& request) const {
        Eigen::MatrixXd gradients =
            Eigen::MatrixXd::Zero(request.designVariables.size(), static_cast<Eigen::Index>(quantities.size()));
        for (std::size_t index = 0; index < quantities.size(); ++index) {
            gradients.col(static_cast<Eigen::Index>(index)) =
                responseSectionGradient(quantities[index].response, request);
        }
        return gradients;
    }

private:
    [[nodiscard]] static Eigen::VectorXd responseGradient(const SeparableLaminateResponse& response,
                                                          const AnalysisRequest& request) {
        return responseDesignGradient(response, request) + responseSectionGradient(response, request);
    }

    [[nodiscard]] static Eigen::VectorXd responseDesignGradient(const SeparableLaminateResponse& response,
                                                                const AnalysisRequest& request) {
        return SeparableLaminateResponseModel::evaluateDenseTermGradient(response.designTerms,
                                                                         request.designVariables);
    }

    [[nodiscard]] static Eigen::VectorXd responseSectionGradient(const SeparableLaminateResponse& response,
                                                                 const AnalysisRequest& request) {
        Eigen::VectorXd gradient = Eigen::VectorXd::Zero(request.designVariables.size());
        for (std::size_t iSection = 0; iSection < response.sectionTerms.size(); ++iSection) {
            const LaminateSectionState& sectionState = request.laminateSections[iSection];
            const SeparableLaminateSectionResponseTerm& sectionTerm = response.sectionTerms[iSection];
            const Eigen::Index sectionSize = sectionTerm.linearCoefficients.size();
            const Eigen::VectorXd sectionVariables =
                request.designVariables.segment(sectionState.variableOffset, sectionSize);

            gradient.segment(sectionState.variableOffset, sectionSize) +=
                sectionTerm.linearCoefficients
                + sectionTerm.quadraticCoefficients.cwiseProduct(sectionVariables);
        }
        return gradient;
    }
};

enum class LaminationParameterCoverageMode {
    CompleteResponse,
    LaminateSectionRowsOnly
};

struct SeparableLaminationParameterDerivativeOptions {
    LaminationParameterCoverageMode coverageMode = LaminationParameterCoverageMode::CompleteResponse;
};

class SeparableLaminationParameterDerivativeProvider final : public LaminationParameterDerivativeProvider {
public:
    explicit SeparableLaminationParameterDerivativeProvider(
        SeparableLaminateResponseModel model,
        SeparableLaminationParameterDerivativeOptions options = {})
        : m_model(std::move(model))
        , m_options(options) {}

    [[nodiscard]] bool supports(const AnalysisRequest& request,
                                const AnalysisResult& result) const override {
        std::string message;
        return m_model.validate(request, result, message);
    }

    [[nodiscard]] LaminationParameterDerivativeResult compute(const AnalysisRequest& request,
                                                              const AnalysisResult& result) const override {
        LaminationParameterDerivativeResult derivativeResult;
        if (!m_model.validate(request, result, derivativeResult.message)) {
            return derivativeResult;
        }

        if (!result.hasObjectiveGradients()) {
            derivativeResult.objectiveGradients =
                makeContribution(request, m_model.objectives.size(), true);
        }
        if (!result.hasConstraintGradients()) {
            derivativeResult.constraintGradients =
                makeContribution(request, m_model.constraints.size(), false);
        }

        derivativeResult.success = true;
        derivativeResult.message = "Optimiser-side lamination-parameter gradients attached.";
        return derivativeResult;
    }

    [[nodiscard]] const SeparableLaminateResponseModel& model() const {
        return m_model;
    }

private:
    [[nodiscard]] GradientContributionMatrix makeContribution(const AnalysisRequest& request,
                                                              const size_t responseCount,
                                                              const bool objectiveResponse) const {
        GradientContributionMatrix contribution;

        if (m_options.coverageMode == LaminationParameterCoverageMode::CompleteResponse) {
            contribution.values = objectiveResponse
                ? m_model.objectiveGradients(request)
                : m_model.constraintGradients(request);
            contribution.mask = Eigen::MatrixXi::Ones(contribution.values.rows(), contribution.values.cols());
            return contribution;
        }

        contribution.values = objectiveResponse
            ? m_model.objectiveSectionGradients(request)
            : m_model.constraintSectionGradients(request);
        contribution.mask =
            Eigen::MatrixXi::Zero(request.designVariables.size(), static_cast<Eigen::Index>(responseCount));

        for (const LaminateSectionState& sectionState : request.laminateSections) {
            const Eigen::Index sectionSize = LaminateSectionVariableCount(sectionState);
            contribution.mask.block(sectionState.variableOffset,
                                    0,
                                    sectionSize,
                                    static_cast<Eigen::Index>(responseCount)).setOnes();
        }

        return contribution;
    }

    SeparableLaminateResponseModel m_model;
    SeparableLaminationParameterDerivativeOptions m_options;
};

class AssembledLaminationParameterDerivativeProvider final : public LaminationParameterDerivativeProvider {
public:
    AssembledLaminationParameterDerivativeProvider(SeparableLaminateQuantityModel quantityModel,
                                                   std::vector<DerivedScalarResponseRule> objectiveRules,
                                                   std::vector<DerivedScalarResponseRule> constraintRules,
                                                   SeparableLaminationParameterDerivativeOptions options = {})
        : m_quantityModel(std::move(quantityModel))
        , m_objectiveRules(std::move(objectiveRules))
        , m_constraintRules(std::move(constraintRules))
        , m_options(options) {}

    [[nodiscard]] bool supports(const AnalysisRequest& request,
                                const AnalysisResult& result) const override {
        std::string message;
        return validate(request, result, message);
    }

    [[nodiscard]] LaminationParameterDerivativeResult compute(const AnalysisRequest& request,
                                                              const AnalysisResult& result) const override {
        LaminationParameterDerivativeResult derivativeResult;
        if (!validate(request, result, derivativeResult.message)) {
            return derivativeResult;
        }

        const Eigen::MatrixXd quantityGradients = m_options.coverageMode == LaminationParameterCoverageMode::CompleteResponse
            ? m_quantityModel.quantityGradients(request)
            : m_quantityModel.quantitySectionGradients(request);

        if (!result.hasObjectiveGradients()) {
            derivativeResult.objectiveGradients =
                assembleContribution(quantityGradients,
                                     request,
                                     result.extractedScalarValues,
                                     m_objectiveRules);
        }
        if (!result.hasConstraintGradients()) {
            derivativeResult.constraintGradients =
                assembleContribution(quantityGradients,
                                     request,
                                     result.extractedScalarValues,
                                     m_constraintRules);
        }

        derivativeResult.success = true;
        derivativeResult.message = "Optimiser-side lamination-parameter gradients attached from extracted FE quantities.";
        return derivativeResult;
    }

private:
    [[nodiscard]] bool validate(const AnalysisRequest& request,
                                const AnalysisResult& result,
                                std::string& message) const {
        if (!m_quantityModel.validate(request, result, message)) {
            return false;
        }
        if (result.objectives.size() != static_cast<Eigen::Index>(m_objectiveRules.size())) {
            message = "Objective count does not match the assembled laminate derivative rules.";
            return false;
        }
        if (result.constraints.size() != static_cast<Eigen::Index>(m_constraintRules.size())) {
            message = "Constraint count does not match the assembled laminate derivative rules.";
            return false;
        }

        for (const DerivedScalarResponseRule& rule : m_objectiveRules) {
            if (!m_quantityModel.findIndex(rule.sourceId).has_value()) {
                message = "Assembled laminate objective rule references an unknown quantity id: " + rule.sourceId;
                return false;
            }
        }
        for (const DerivedScalarResponseRule& rule : m_constraintRules) {
            if (!m_quantityModel.findIndex(rule.sourceId).has_value()) {
                message = "Assembled laminate constraint rule references an unknown quantity id: " + rule.sourceId;
                return false;
            }
        }

        return true;
    }

    [[nodiscard]] GradientContributionMatrix
    assembleContribution(const Eigen::MatrixXd& quantityGradients,
                         const AnalysisRequest& request,
                         const NamedScalarResponseMap& extractedScalars,
                         const std::vector<DerivedScalarResponseRule>& rules) const {
        GradientContributionMatrix contribution;
        contribution.values =
            Eigen::MatrixXd::Zero(request.designVariables.size(), static_cast<Eigen::Index>(rules.size()));
        contribution.mask =
            Eigen::MatrixXi::Zero(request.designVariables.size(), static_cast<Eigen::Index>(rules.size()));

        for (std::size_t iRule = 0; iRule < rules.size(); ++iRule) {
            const auto quantityIndex = m_quantityModel.findIndex(rules[iRule].sourceId);
            const double chainFactor = EvaluateDerivedScalarResponseDerivative(rules[iRule], extractedScalars);
            contribution.values.col(static_cast<Eigen::Index>(iRule)) =
                chainFactor * quantityGradients.col(static_cast<Eigen::Index>(*quantityIndex));
        }

        if (m_options.coverageMode == LaminationParameterCoverageMode::CompleteResponse) {
            contribution.mask.setOnes();
        } else {
            for (const LaminateSectionState& sectionState : request.laminateSections) {
                const Eigen::Index sectionSize = LaminateSectionVariableCount(sectionState);
                contribution.mask.block(sectionState.variableOffset,
                                        0,
                                        sectionSize,
                                        contribution.mask.cols()).setOnes();
            }
        }

        return contribution;
    }

    SeparableLaminateQuantityModel m_quantityModel;
    std::vector<DerivedScalarResponseRule> m_objectiveRules;
    std::vector<DerivedScalarResponseRule> m_constraintRules;
    SeparableLaminationParameterDerivativeOptions m_options;
};

}  // namespace lamopt
