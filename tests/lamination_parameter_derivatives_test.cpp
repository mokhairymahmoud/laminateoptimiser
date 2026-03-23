#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"
#include "../src/OptimisationPipeline/laminationParameterDerivatives.hpp"

#include <gtest/gtest.h>

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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
