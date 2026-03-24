#include "../src/OptimisationPipeline/calculixComposipyBenchmark.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"
#include "../src/OptimisationPipeline/laminationParameterDerivatives.hpp"
#include "../src/OptimisationPipeline/responseSchema.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <utility>

namespace {

std::filesystem::path TempPath(const std::string& name) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("laminateoptimiser_" + name);
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
}

lamopt::LaminateSectionState MakeBalancedSymmetricSection(const Eigen::Index variableOffset) {
    lamopt::LaminateSectionState laminateSection;
    laminateSection.variableOffset = variableOffset;
    laminateSection.isBalanced = true;
    laminateSection.isSymmetric = true;
    laminateSection.sublaminateCount = 6;
    laminateSection.thicknessLowerBound = 0.5;
    laminateSection.thicknessUpperBound = 2.0;
    laminateSection.thickness = 1.1;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    return laminateSection;
}

lamopt::SeparableLaminateResponseModel MakeBenchmarkModel() {
    const Eigen::Index designSize = 7;
    const lamopt::LaminateSectionState section = MakeBalancedSymmetricSection(1);
    const Eigen::Index sectionSize = lamopt::LaminateSectionVariableCount(section);

    lamopt::DenseDesignResponseTerm objectiveDesign;
    objectiveDesign.linearCoefficients = Eigen::VectorXd::Zero(designSize);
    objectiveDesign.quadraticCoefficients = Eigen::VectorXd::Zero(designSize);
    objectiveDesign.linearCoefficients(0) = -0.6;
    objectiveDesign.linearCoefficients(6) = 0.2;
    objectiveDesign.quadraticCoefficients(0) = 0.1;
    objectiveDesign.quadraticCoefficients(6) = -0.05;

    lamopt::SeparableLaminateSectionResponseTerm objectiveTerm;
    objectiveTerm.linearCoefficients = Eigen::VectorXd::Zero(sectionSize);
    objectiveTerm.quadraticCoefficients = Eigen::VectorXd::Zero(sectionSize);
    objectiveTerm.linearCoefficients << 0.2, -0.4, 0.1, 0.3, -1.0;
    objectiveTerm.quadraticCoefficients << 0.5, 0.0, -0.25, 0.0, 0.8;

    lamopt::DenseDesignResponseTerm constraintDesign;
    constraintDesign.linearCoefficients = Eigen::VectorXd::Zero(designSize);
    constraintDesign.quadraticCoefficients = Eigen::VectorXd::Zero(designSize);
    constraintDesign.linearCoefficients(0) = 0.05;
    constraintDesign.linearCoefficients(6) = -0.1;
    constraintDesign.quadraticCoefficients(0) = 0.02;
    constraintDesign.quadraticCoefficients(6) = 0.03;

    lamopt::SeparableLaminateSectionResponseTerm constraintTerm;
    constraintTerm.linearCoefficients = Eigen::VectorXd::Zero(sectionSize);
    constraintTerm.quadraticCoefficients = Eigen::VectorXd::Zero(sectionSize);
    constraintTerm.linearCoefficients << 0.1, 0.0, 0.2, -0.1, 0.15;
    constraintTerm.quadraticCoefficients << 0.0, 0.4, 0.0, 0.1, 0.2;

    lamopt::SeparableLaminateResponseModel model;
    model.objectives.push_back({0.7, objectiveDesign, {objectiveTerm}});
    model.constraints.push_back({-1.2, constraintDesign, {constraintTerm}});
    return model;
}

lamopt::SeparableLaminateQuantityModel MakeAssembledQuantityModel() {
    const Eigen::Index designSize = 7;
    const lamopt::LaminateSectionState section = MakeBalancedSymmetricSection(1);
    const Eigen::Index sectionSize = lamopt::LaminateSectionVariableCount(section);

    auto makeResponse = [designSize, sectionSize](double constant,
                                                  std::initializer_list<double> designLinear,
                                                  std::initializer_list<double> designQuadratic,
                                                  std::initializer_list<double> sectionLinear,
                                                  std::initializer_list<double> sectionQuadratic) {
        lamopt::DenseDesignResponseTerm designTerms;
        designTerms.linearCoefficients = Eigen::VectorXd::Zero(designSize);
        designTerms.quadraticCoefficients = Eigen::VectorXd::Zero(designSize);

        Eigen::Index index = 0;
        for (double value : designLinear) {
            designTerms.linearCoefficients(index++) = value;
        }
        index = 0;
        for (double value : designQuadratic) {
            designTerms.quadraticCoefficients(index++) = value;
        }

        lamopt::SeparableLaminateSectionResponseTerm sectionTerms;
        sectionTerms.linearCoefficients = Eigen::VectorXd::Zero(sectionSize);
        sectionTerms.quadraticCoefficients = Eigen::VectorXd::Zero(sectionSize);
        index = 0;
        for (double value : sectionLinear) {
            sectionTerms.linearCoefficients(index++) = value;
        }
        index = 0;
        for (double value : sectionQuadratic) {
            sectionTerms.quadraticCoefficients(index++) = value;
        }

        return lamopt::SeparableLaminateResponse{constant, designTerms, {sectionTerms}};
    };

    lamopt::SeparableLaminateQuantityModel model;
    model.quantities.push_back({"mass", makeResponse(9.0,
                                                     {0.3, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1},
                                                     {0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03},
                                                     {0.1, -0.05, 0.08, 0.0, 0.4},
                                                     {0.0, 0.02, 0.0, 0.0, 0.1})});
    model.quantities.push_back({"buckling_lambda_1", makeResponse(2.5,
                                                                  {-0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15},
                                                                  {0.04, 0.0, 0.0, 0.0, 0.0, 0.0, -0.01},
                                                                  {0.2, 0.05, -0.1, 0.04, 0.25},
                                                                  {0.03, 0.0, 0.01, 0.0, 0.05})});
    model.quantities.push_back({"tip_u3", makeResponse(0.6,
                                                       {0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2},
                                                       {0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                                       {0.05, -0.03, 0.04, 0.02, 0.3},
                                                       {0.0, 0.01, 0.0, 0.0, 0.04})});
    return model;
}

std::vector<lamopt::DerivedScalarResponseRule> MakeObjectiveAssemblyRules() {
    return {
        {"mass", lamopt::DerivedScalarResponseTransform::Identity, 1.0, 0.0, "objective_mass"}
    };
}

std::vector<lamopt::DerivedScalarResponseRule> MakeConstraintAssemblyRules() {
    return {
        {"buckling_lambda_1", lamopt::DerivedScalarResponseTransform::InverseAffine, 1.0, -1.0, "buckling_margin"},
        {"tip_u3", lamopt::DerivedScalarResponseTransform::AbsoluteAffine, 1.0 / 3.0, -1.0, "tip_displacement_margin"}
    };
}

class ZeroMovementSolver final : public lamopt::SubproblemSolver {
public:
    lamopt::SubproblemResult solve(const lamopt::ApproximationProblem& problem) override {
        lamopt::SubproblemResult result;
        result.success = true;
        result.iterations = 1;
        result.candidateDesign = problem.referenceDesign;
        result.predictedObjectives = problem.objectiveValues;
        result.predictedConstraints = problem.constraintValues;
        result.message = "zero-movement candidate";
        return result;
    }
};

class SeparableLaminateBackend final : public lamopt::AnalysisBackend {
public:
    enum class GradientMode {
        None,
        DesignOnlyPartial,
        ObjectiveOnlyBiased,
        All
    };

    SeparableLaminateBackend(lamopt::SeparableLaminateResponseModel model,
                             GradientMode gradientMode)
        : m_model(std::move(model))
        , m_gradientMode(gradientMode) {}

    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        ++evaluationCount;
        if (!request.workDirectory.empty()) {
            std::filesystem::create_directories(request.workDirectory);
        }

        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;
        result.objectives = m_model.evaluateObjectives(request);
        result.constraints = m_model.evaluateConstraints(request);
        result.diagnostics.message = "separable laminate benchmark";

        if (m_gradientMode == GradientMode::All || m_gradientMode == GradientMode::ObjectiveOnlyBiased) {
            result.objectiveGradients = m_model.objectiveGradients(request);
        }
        if (m_gradientMode == GradientMode::ObjectiveOnlyBiased) {
            (*result.objectiveGradients)(request.laminateSections.front().variableOffset, 0) += 0.125;
        }
        if (m_gradientMode == GradientMode::DesignOnlyPartial) {
            result.objectiveGradients = m_model.objectiveDesignGradients(request);
            result.objectiveGradientMask = makeDesignOnlyMask(request, result.objectives.size());
            result.constraintGradients = m_model.constraintDesignGradients(request);
            result.constraintGradientMask = makeDesignOnlyMask(request, result.constraints.size());
        }
        if (m_gradientMode == GradientMode::All) {
            result.constraintGradients = m_model.constraintGradients(request);
        }

        return result;
    }

    int evaluationCount = 0;

private:
    static Eigen::MatrixXi makeDesignOnlyMask(const lamopt::AnalysisRequest& request,
                                              const Eigen::Index responseCount) {
        Eigen::MatrixXi mask = Eigen::MatrixXi::Ones(request.designVariables.size(), responseCount);
        const lamopt::LaminateSectionState& sectionState = request.laminateSections.front();
        const Eigen::Index sectionSize = lamopt::LaminateSectionVariableCount(sectionState);
        mask.block(sectionState.variableOffset, 0, sectionSize, responseCount).setZero();
        return mask;
    }

    lamopt::SeparableLaminateResponseModel m_model;
    GradientMode m_gradientMode = GradientMode::None;
};

class AssembledLaminateBackend final : public lamopt::AnalysisBackend {
public:
    enum class GradientMode {
        None,
        DesignOnlyPartial
    };

    AssembledLaminateBackend(lamopt::SeparableLaminateQuantityModel quantityModel,
                             std::vector<lamopt::DerivedScalarResponseRule> objectiveRules,
                             std::vector<lamopt::DerivedScalarResponseRule> constraintRules,
                             GradientMode gradientMode)
        : m_quantityModel(std::move(quantityModel))
        , m_objectiveRules(std::move(objectiveRules))
        , m_constraintRules(std::move(constraintRules))
        , m_gradientMode(gradientMode) {}

    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;

        const Eigen::VectorXd rawValues = m_quantityModel.evaluateQuantities(request);
        for (std::size_t index = 0; index < m_quantityModel.quantities.size(); ++index) {
            result.extractedScalarValues.emplace(m_quantityModel.quantities[index].id,
                                                 rawValues(static_cast<Eigen::Index>(index)));
        }

        result.objectives = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(m_objectiveRules.size()));
        for (std::size_t index = 0; index < m_objectiveRules.size(); ++index) {
            result.objectives(static_cast<Eigen::Index>(index)) =
                lamopt::EvaluateDerivedScalarResponseRule(m_objectiveRules[index], result.extractedScalarValues);
        }

        result.constraints = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(m_constraintRules.size()));
        for (std::size_t index = 0; index < m_constraintRules.size(); ++index) {
            result.constraints(static_cast<Eigen::Index>(index)) =
                lamopt::EvaluateDerivedScalarResponseRule(m_constraintRules[index], result.extractedScalarValues);
        }

        if (m_gradientMode == GradientMode::DesignOnlyPartial) {
            result.objectiveGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), result.objectives.size());
            result.constraintGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), result.constraints.size());
            result.objectiveGradientMask = makeDesignOnlyMask(request, result.objectives.size());
            result.constraintGradientMask = makeDesignOnlyMask(request, result.constraints.size());
        }

        result.diagnostics.message = "assembled laminate benchmark";
        return result;
    }

private:
    static Eigen::MatrixXi makeDesignOnlyMask(const lamopt::AnalysisRequest& request,
                                              const Eigen::Index responseCount) {
        Eigen::MatrixXi mask = Eigen::MatrixXi::Ones(request.designVariables.size(), responseCount);
        const lamopt::LaminateSectionState& sectionState = request.laminateSections.front();
        const Eigen::Index sectionSize = lamopt::LaminateSectionVariableCount(sectionState);
        mask.block(sectionState.variableOffset, 0, sectionSize, responseCount).setZero();
        return mask;
    }

    lamopt::SeparableLaminateQuantityModel m_quantityModel;
    std::vector<lamopt::DerivedScalarResponseRule> m_objectiveRules;
    std::vector<lamopt::DerivedScalarResponseRule> m_constraintRules;
    GradientMode m_gradientMode = GradientMode::None;
};

class ComposipyPowerLawBackend final : public lamopt::AnalysisBackend {
public:
    explicit ComposipyPowerLawBackend(lamopt::CalculixComposipyBenchmarkSpec spec = {})
        : m_spec(std::move(spec)) {}

    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;

        const double thickness = request.designVariables(0);
        const double massPerArea = static_cast<double>(m_spec.plyAngles.size()) * thickness;
        const double buckling = 2.12358 * std::pow(thickness / m_spec.referencePlyThickness, 3.0);
        const double tip = -0.714186 * std::pow(m_spec.referencePlyThickness / thickness, 3.0);

        result.extractedScalarValues["mass_per_area"] = massPerArea;
        result.extractedScalarValues["buckling_lambda_1"] = buckling;
        result.extractedScalarValues["tip_u3"] = tip;

        const auto objectiveRules = lamopt::MakeCalculixComposipyBenchmarkObjectiveRules(m_spec);
        const auto constraintRules = lamopt::MakeCalculixComposipyBenchmarkConstraintRules(m_spec);
        result.objectives = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(objectiveRules.size()));
        result.constraints = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(constraintRules.size()));

        for (std::size_t index = 0; index < objectiveRules.size(); ++index) {
            result.objectives(static_cast<Eigen::Index>(index)) =
                lamopt::EvaluateDerivedScalarResponseRule(objectiveRules[index], result.extractedScalarValues);
        }
        for (std::size_t index = 0; index < constraintRules.size(); ++index) {
            result.constraints(static_cast<Eigen::Index>(index)) =
                lamopt::EvaluateDerivedScalarResponseRule(constraintRules[index], result.extractedScalarValues);
        }

        result.sensitivityPolicy = lamopt::MakeCalculixComposipyBenchmarkSensitivityPolicy(m_spec);
        return result;
    }

private:
    lamopt::CalculixComposipyBenchmarkSpec m_spec;
};

lamopt::AnalysisRequest MakeBenchmarkRequest(const std::string& name) {
    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Zero(7);
    request.designVariables << 2.0, 0.1, -0.2, 0.05, 0.3, 1.1, -3.0;
    request.requestSensitivities = true;
    request.workDirectory = TempPath(name);
    request.laminateSections.push_back(MakeBalancedSymmetricSection(1));
    return request;
}

lamopt::DriverOptions MakeDriverOptions() {
    lamopt::DriverOptions options;
    options.maxOuterIterations = 1;
    options.maxSubIterations = 1;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-12;
    options.finiteDifferenceStep = 1.0e-7;
    return options;
}

}  // namespace

TEST(LaminationParameterDerivativesTest, OptimiserSideProviderMatchesFiniteDifferenceBenchmark) {
    const lamopt::SeparableLaminateResponseModel model = MakeBenchmarkModel();
    lamopt::SeparableLaminationParameterDerivativeProvider derivativeProvider(model);
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions providerOptions = MakeDriverOptions();
    providerOptions.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    SeparableLaminateBackend providerBackend(model, SeparableLaminateBackend::GradientMode::None);
    lamopt::GlobalOptimisationDriver providerDriver(providerBackend,
                                                    approximationBuilder,
                                                    subproblemSolver,
                                                    providerOptions,
                                                    &derivativeProvider);

    const lamopt::AnalysisRequest providerRequest = MakeBenchmarkRequest("laminate_derivative_provider");
    const lamopt::GlobalOptimisationResult providerResult = providerDriver.optimise(providerRequest);

    ASSERT_TRUE(providerResult.converged);
    ASSERT_TRUE(providerResult.analysis.hasAllGradients());
    ASSERT_TRUE(providerResult.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(providerResult.analysis.constraintGradients.has_value());
    EXPECT_TRUE(providerResult.analysis.objectiveGradients->isApprox(model.objectiveGradients(providerRequest), 1.0e-12));
    EXPECT_TRUE(providerResult.analysis.constraintGradients->isApprox(model.constraintGradients(providerRequest), 1.0e-12));

    lamopt::DriverOptions finiteDifferenceOptions = MakeDriverOptions();
    finiteDifferenceOptions.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;

    SeparableLaminateBackend finiteDifferenceBackend(model, SeparableLaminateBackend::GradientMode::None);
    lamopt::GlobalOptimisationDriver finiteDifferenceDriver(finiteDifferenceBackend,
                                                            approximationBuilder,
                                                            subproblemSolver,
                                                            finiteDifferenceOptions);

    const lamopt::AnalysisRequest finiteDifferenceRequest =
        MakeBenchmarkRequest("laminate_derivative_fd");
    const lamopt::GlobalOptimisationResult finiteDifferenceResult =
        finiteDifferenceDriver.optimise(finiteDifferenceRequest);

    ASSERT_TRUE(finiteDifferenceResult.converged);
    ASSERT_TRUE(finiteDifferenceResult.analysis.hasAllGradients());
    ASSERT_TRUE(finiteDifferenceResult.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(finiteDifferenceResult.analysis.constraintGradients.has_value());

    EXPECT_TRUE(providerResult.analysis.objectiveGradients->isApprox(
        *finiteDifferenceResult.analysis.objectiveGradients,
        1.0e-5));
    EXPECT_TRUE(providerResult.analysis.constraintGradients->isApprox(
        *finiteDifferenceResult.analysis.constraintGradients,
        1.0e-5));
    EXPECT_LT(providerBackend.evaluationCount, finiteDifferenceBackend.evaluationCount);
}

TEST(LaminationParameterDerivativesTest, ProviderPreservesBackendGradientsAndFillsMissingOnes) {
    const lamopt::SeparableLaminateResponseModel model = MakeBenchmarkModel();
    lamopt::SeparableLaminationParameterDerivativeProvider derivativeProvider(model);
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions options = MakeDriverOptions();
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    SeparableLaminateBackend backend(model, SeparableLaminateBackend::GradientMode::ObjectiveOnlyBiased);
    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);

    const lamopt::AnalysisRequest request = MakeBenchmarkRequest("laminate_derivative_partial");
    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(result.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(result.analysis.constraintGradients.has_value());

    Eigen::MatrixXd biasedObjectiveGradients = model.objectiveGradients(request);
    biasedObjectiveGradients(request.laminateSections.front().variableOffset, 0) += 0.125;

    EXPECT_TRUE(result.analysis.objectiveGradients->isApprox(biasedObjectiveGradients, 1.0e-12));
    EXPECT_TRUE(result.analysis.constraintGradients->isApprox(model.constraintGradients(request), 1.0e-12));
    EXPECT_NE(result.analysis.diagnostics.message.find("Optimiser-side lamination-parameter gradients attached."),
              std::string::npos);
}

TEST(LaminationParameterDerivativesTest, ProviderCombinesLaminateRowsWithPartialBackendGradients) {
    const lamopt::SeparableLaminateResponseModel model = MakeBenchmarkModel();
    lamopt::SeparableLaminationParameterDerivativeOptions derivativeOptions;
    derivativeOptions.coverageMode = lamopt::LaminationParameterCoverageMode::LaminateSectionRowsOnly;
    lamopt::SeparableLaminationParameterDerivativeProvider derivativeProvider(model, derivativeOptions);
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions options = MakeDriverOptions();
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    SeparableLaminateBackend backend(model, SeparableLaminateBackend::GradientMode::DesignOnlyPartial);
    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);

    const lamopt::AnalysisRequest request = MakeBenchmarkRequest("laminate_derivative_assembled");
    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(result.analysis.hasAllGradients());
    ASSERT_TRUE(result.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(result.analysis.constraintGradients.has_value());
    ASSERT_TRUE(result.analysis.objectiveGradientMask.has_value());
    ASSERT_TRUE(result.analysis.constraintGradientMask.has_value());

    EXPECT_TRUE(result.analysis.objectiveGradients->isApprox(model.objectiveGradients(request), 1.0e-12));
    EXPECT_TRUE(result.analysis.constraintGradients->isApprox(model.constraintGradients(request), 1.0e-12));
    EXPECT_EQ(result.analysis.objectiveGradientMask->minCoeff(), 1);
    EXPECT_EQ(result.analysis.constraintGradientMask->minCoeff(), 1);
    EXPECT_NE(result.analysis.diagnostics.message.find("Optimiser-side lamination-parameter gradients attached."),
              std::string::npos);
}

TEST(LaminationParameterDerivativesTest, AssembledProviderBindsExtractedQuantitiesThroughChainRule) {
    const lamopt::SeparableLaminateQuantityModel quantityModel = MakeAssembledQuantityModel();
    const auto objectiveRules = MakeObjectiveAssemblyRules();
    const auto constraintRules = MakeConstraintAssemblyRules();

    lamopt::AssembledLaminationParameterDerivativeProvider derivativeProvider(
        quantityModel,
        objectiveRules,
        constraintRules
    );
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions options = MakeDriverOptions();
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    AssembledLaminateBackend backend(quantityModel,
                                     objectiveRules,
                                     constraintRules,
                                     AssembledLaminateBackend::GradientMode::None);
    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);

    const lamopt::AnalysisRequest request = MakeBenchmarkRequest("laminate_derivative_schema_binding");
    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(result.analysis.hasAllGradients());
    ASSERT_TRUE(result.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(result.analysis.constraintGradients.has_value());

    const Eigen::MatrixXd quantityGradients = quantityModel.quantityGradients(request);
    const double bucklingValue = result.analysis.extractedScalarValues.at("buckling_lambda_1");
    const double tipValue = result.analysis.extractedScalarValues.at("tip_u3");

    Eigen::MatrixXd expectedObjectiveGradients(request.designVariables.size(), 1);
    expectedObjectiveGradients.col(0) = quantityGradients.col(0);

    Eigen::MatrixXd expectedConstraintGradients(request.designVariables.size(), 2);
    expectedConstraintGradients.col(0) = (-1.0 / (bucklingValue * bucklingValue)) * quantityGradients.col(1);
    expectedConstraintGradients.col(1) = (1.0 / 3.0) * quantityGradients.col(2);

    EXPECT_TRUE(result.analysis.objectiveGradients->isApprox(expectedObjectiveGradients, 1.0e-12));
    EXPECT_TRUE(result.analysis.constraintGradients->isApprox(expectedConstraintGradients, 1.0e-12));
    EXPECT_NE(result.analysis.diagnostics.message.find("extracted FE quantities"), std::string::npos);
}

TEST(LaminationParameterDerivativesTest, AssembledProviderCanFillLaminateRowsOnlyForSchemaResponses) {
    const lamopt::SeparableLaminateQuantityModel quantityModel = MakeAssembledQuantityModel();
    const auto objectiveRules = MakeObjectiveAssemblyRules();
    const auto constraintRules = MakeConstraintAssemblyRules();

    lamopt::SeparableLaminationParameterDerivativeOptions derivativeOptions;
    derivativeOptions.coverageMode = lamopt::LaminationParameterCoverageMode::LaminateSectionRowsOnly;

    lamopt::AssembledLaminationParameterDerivativeProvider derivativeProvider(
        quantityModel,
        objectiveRules,
        constraintRules,
        derivativeOptions
    );
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions options = MakeDriverOptions();
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    AssembledLaminateBackend backend(quantityModel,
                                     objectiveRules,
                                     constraintRules,
                                     AssembledLaminateBackend::GradientMode::DesignOnlyPartial);
    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);

    const lamopt::AnalysisRequest request = MakeBenchmarkRequest("laminate_derivative_schema_partial");
    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(result.analysis.hasAllGradients());
    ASSERT_TRUE(result.analysis.objectiveGradientMask.has_value());
    ASSERT_TRUE(result.analysis.constraintGradientMask.has_value());
    EXPECT_EQ(result.analysis.objectiveGradientMask->minCoeff(), 1);
    EXPECT_EQ(result.analysis.constraintGradientMask->minCoeff(), 1);
}

TEST(LaminationParameterDerivativesTest, ProductionComposipyProviderBindsRealExtractedQuantityRules) {
    const lamopt::CalculixComposipyBenchmarkSpec spec;
    auto derivativeProvider = lamopt::MakeCalculixComposipyBenchmarkDerivativeProvider(spec);
    ComposipyPowerLawBackend backend(spec);
    lamopt::LinearApproximationBuilder approximationBuilder;
    ZeroMovementSolver subproblemSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 1;
    options.maxSubIterations = 1;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;
    options.stagnationTolerance = 1.0e-12;

    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, spec.referencePlyThickness);
    request.requestSensitivities = true;
    request.workDirectory = TempPath("laminate_derivative_composipy_provider");

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(result.analysis.hasAllGradients());
    ASSERT_TRUE(result.analysis.objectiveGradients.has_value());
    ASSERT_TRUE(result.analysis.constraintGradients.has_value());

    const double thickness = request.designVariables(0);
    const double massPerArea = result.analysis.extractedScalarValues.at("mass_per_area");
    const double buckling = result.analysis.extractedScalarValues.at("buckling_lambda_1");
    const double tip = result.analysis.extractedScalarValues.at("tip_u3");

    Eigen::MatrixXd expectedObjectiveGradients(1, 1);
    expectedObjectiveGradients(0, 0) = massPerArea / thickness;

    Eigen::MatrixXd expectedConstraintGradients(1, 2);
    expectedConstraintGradients(0, 0) =
        (-spec.requiredBucklingFactor / (buckling * buckling)) * (3.0 * buckling / thickness);
    expectedConstraintGradients(0, 1) =
        (-1.0 / spec.tipDisplacementLimit) * (-3.0 * tip / thickness);

    EXPECT_TRUE(result.analysis.objectiveGradients->isApprox(expectedObjectiveGradients, 1.0e-12));
    EXPECT_TRUE(result.analysis.constraintGradients->isApprox(expectedConstraintGradients, 1.0e-12));
    EXPECT_NE(result.analysis.diagnostics.message.find("power-law rules"), std::string::npos);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
