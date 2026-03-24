#include "../src/OptimisationPipeline/defaultLaminateSubproblem.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"
#include "../src/OptimisationPipeline/laminateAwareSubproblem.hpp"
#include "../src/OptimisationPipeline/coreMinMaxSubproblem.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>

namespace {

class QuadraticBackend final : public lamopt::AnalysisBackend {
public:
    explicit QuadraticBackend(bool provideGradients)
        : m_provideGradients(provideGradients) {}

    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        ++evaluationCount;
        if (!request.workDirectory.empty()) {
            std::filesystem::create_directories(request.workDirectory);
        }

        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;
        result.objectives = Eigen::VectorXd::Constant(1, request.designVariables.squaredNorm());
        result.constraints = Eigen::VectorXd::Constant(1, request.designVariables.sum() - 10.0);
        result.diagnostics.message = "quadratic backend";

        if (m_provideGradients) {
            Eigen::MatrixXd objectiveGradient(request.designVariables.size(), 1);
            objectiveGradient.col(0) = 2.0 * request.designVariables;
            result.objectiveGradients = objectiveGradient;

            Eigen::MatrixXd constraintGradient(request.designVariables.size(), 1);
            constraintGradient.col(0) = Eigen::VectorXd::Ones(request.designVariables.size());
            result.constraintGradients = constraintGradient;
        }

        return result;
    }

    int evaluationCount = 0;

private:
    bool m_provideGradients = true;
};

class PolicyAwarePartialGradientBackend final : public lamopt::AnalysisBackend {
public:
    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;
        result.objectives = Eigen::VectorXd::Constant(1, request.designVariables(0));
        result.constraints = Eigen::VectorXd::Constant(1, request.designVariables(0) - 2.0);
        result.objectiveGradients = Eigen::MatrixXd::Constant(1, 1, 1.0);
        result.objectiveGradientMask = Eigen::MatrixXi::Ones(1, 1);
        result.sensitivityPolicy = lamopt::ResponseSensitivityPolicy{
            {lamopt::ResponseSensitivitySource::BackendNative},
            {lamopt::ResponseSensitivitySource::FiniteDifference}
        };
        result.diagnostics.message = "policy-aware backend";
        return result;
    }
};

class ScriptedSolver final : public lamopt::SubproblemSolver {
public:
    lamopt::SubproblemResult solve(const lamopt::ApproximationProblem& problem) override {
        lamopt::SubproblemResult result;
        result.success = true;
        result.iterations = 1;

        if (callCount == 0) {
            result.candidateDesign = problem.referenceDesign + Eigen::VectorXd::Constant(problem.referenceDesign.size(), 1.0);
            result.predictedObjectives = Eigen::VectorXd::Constant(1, 1.0);
            result.predictedConstraints = Eigen::VectorXd::Constant(1, -1.0);
            result.message = "forced rejection";
        } else {
            result.candidateDesign = 0.5 * problem.referenceDesign;
            result.predictedObjectives = Eigen::VectorXd::Constant(1, 0.25);
            result.predictedConstraints = Eigen::VectorXd::Constant(1, -1.0);
            result.message = "forced acceptance";
        }

        ++callCount;
        return result;
    }

    int callCount = 0;
};

class HalvingSolver final : public lamopt::SubproblemSolver {
public:
    lamopt::SubproblemResult solve(const lamopt::ApproximationProblem& problem) override {
        lamopt::SubproblemResult result;
        result.success = true;
        result.iterations = 1;
        result.candidateDesign = 0.5 * problem.referenceDesign;
        result.predictedObjectives = 0.25 * problem.objectiveValues;
        result.predictedConstraints = Eigen::VectorXd::Constant(problem.constraintValues.size(), -1.0);
        result.message = "halving candidate";
        return result;
    }
};

class LaminateScriptedSolver final : public lamopt::SubproblemSolver {
public:
    lamopt::SubproblemResult solve(const lamopt::ApproximationProblem& problem) override {
        lamopt::SubproblemResult result;
        result.success = true;
        result.iterations = 1;
        result.candidateDesign = problem.referenceDesign;
        result.candidateDesign(0) = 5.0;
        result.candidateDesign(1) = 5.0;
        result.candidateDesign(2) = 5.0;
        result.candidateDesign(3) = 5.0;
        result.candidateDesign(4) = 5.0;
        result.predictedObjectives = Eigen::VectorXd::Constant(1, -5.0);
        result.predictedConstraints = Eigen::VectorXd();
        result.message = "laminate scripted candidate";
        return result;
    }
};

class Problem1Backend final : public lamopt::AnalysisBackend {
public:
    lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
        lamopt::AnalysisResult result;
        result.status = lamopt::AnalysisStatus::Success;

        const double x1 = request.designVariables(0);
        const double x2 = request.designVariables(1);
        const double sqrt2 = std::sqrt(2.0);
        const double denominator = 2.0 * x1 * x2 + sqrt2 * x1 * x1;
        const double aux = sqrt2 * x1 * x2 + x1 * x1;
        const double common = sqrt2 * x2 + x1;

        result.objectives = Eigen::VectorXd::Constant(1, 100.0 * (2.0 * sqrt2 * x1 + x2));
        result.constraints = Eigen::VectorXd(3);
        result.constraints << (x2 + sqrt2 * x1) / denominator - 1.0,
                              1.0 / (x1 + sqrt2 * x2) - 1.0,
                              4.0 / 3.0 * x2 / denominator - 1.0;

        if (request.requestSensitivities) {
            Eigen::MatrixXd objectiveGradient(2, 1);
            objectiveGradient << 200.0 * sqrt2, 100.0;
            result.objectiveGradients = objectiveGradient;

            Eigen::MatrixXd constraintGradients(2, 3);
            constraintGradients(0, 0) = -(sqrt2 * x1 * x2 + x1 * x1 + x2 * x2) / std::pow(aux, 2);
            constraintGradients(1, 0) = -1.0 / (sqrt2 * std::pow(common, 2));
            constraintGradients(0, 1) = -1.0 / std::pow(common, 2);
            constraintGradients(1, 1) = -sqrt2 / std::pow(common, 2);
            constraintGradients(0, 2) = -4.0 / 3.0 * x2 * (x2 + sqrt2 * x1) / std::pow(aux, 2);
            constraintGradients(1, 2) = 4.0 / 3.0 * 1.0 / (sqrt2 * std::pow(common, 2));
            result.constraintGradients = constraintGradients;
        }

        return result;
    }
};

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
    laminateSection.thickness = 1.0;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    return laminateSection;
}

}  // namespace

TEST(GlobalDriverTest, DriverConvergesWithDirectGradients) {
    QuadraticBackend backend(true);
    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::GradientPenaltySubproblemSolver subproblemSolver({4.0, 5.0, 1.0e-12});

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.stagnationTolerance = 5.0e-2;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(2.0, -1.0);
    request.workDirectory = TempPath("driver_direct_gradients");
    request.lowerBounds = Eigen::Vector2d::Constant(-5.0);
    request.upperBounds = Eigen::Vector2d::Constant(5.0);

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.design.norm(), 0.3);
    EXPECT_FALSE(result.history.empty());
}

TEST(GlobalDriverTest, DriverUsesFiniteDifferenceFallbackWhenGradientsAreMissing) {
    QuadraticBackend backend(false);
    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::GradientPenaltySubproblemSolver subproblemSolver({4.0, 5.0, 1.0e-12});

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;
    options.stagnationTolerance = 5.0e-2;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, -0.5);
    request.workDirectory = TempPath("driver_fd_gradients");

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(result.analysis.hasAllGradients());
    EXPECT_GT(backend.evaluationCount, 3);
}

TEST(GlobalDriverTest, DriverReportsConfiguredSensitivityPolicyWhenFallbackIsDisabled) {
    PolicyAwarePartialGradientBackend backend;
    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::GradientPenaltySubproblemSolver subproblemSolver({4.0, 5.0, 1.0e-12});

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::Disabled;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, 1.0);
    request.workDirectory = TempPath("driver_policy_disabled");

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_EQ(result.analysis.status, lamopt::AnalysisStatus::MissingGradients);
    EXPECT_NE(result.analysis.diagnostics.message.find("Configured backend sensitivity policy: objective[0]=backend_native; constraint[0]=finite_difference."),
              std::string::npos);
}

TEST(GlobalDriverTest, DriverCanResumeFromCheckpointAndMatchUninterruptedRun) {
    QuadraticBackend partialBackend(true);
    lamopt::LinearApproximationBuilder approximationBuilder;
    HalvingSolver subproblemSolver;

    lamopt::DriverOptions partialOptions;
    partialOptions.maxOuterIterations = 2;
    partialOptions.maxSubIterations = 1;
    partialOptions.requestSensitivities = true;
    partialOptions.stagnationTolerance = 1.0e-12;

    lamopt::GlobalOptimisationDriver partialDriver(
        partialBackend,
        approximationBuilder,
        subproblemSolver,
        partialOptions
    );

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, 1.0);
    request.workDirectory = TempPath("driver_resume_partial");

    const lamopt::GlobalOptimisationResult partialResult = partialDriver.optimise(request);
    ASSERT_FALSE(partialResult.converged);
    ASSERT_EQ(partialResult.history.size(), 2U);

    const std::filesystem::path checkpointPath = request.workDirectory / "driver.chk";
    lamopt::WriteCheckpoint(checkpointPath, partialResult);
    const lamopt::CheckpointState checkpoint = lamopt::ReadCheckpoint(checkpointPath);

    EXPECT_FALSE(checkpoint.converged);
    EXPECT_EQ(checkpoint.status, lamopt::AnalysisStatus::Success);
    EXPECT_EQ(checkpoint.historyCount, 2U);
    EXPECT_TRUE(checkpoint.design.isApprox(partialResult.design, 1.0e-12));

    QuadraticBackend resumedBackend(true);
    HalvingSolver resumedSolver;
    lamopt::DriverOptions resumedOptions = partialOptions;
    lamopt::GlobalOptimisationDriver resumedDriver(
        resumedBackend,
        approximationBuilder,
        resumedSolver,
        resumedOptions
    );

    lamopt::AnalysisRequest resumedTemplate = request;
    resumedTemplate.workDirectory = TempPath("driver_resume_continued");
    const lamopt::GlobalOptimisationResult resumedResult =
        resumedDriver.optimiseFromCheckpoint(checkpointPath, resumedTemplate);

    QuadraticBackend directBackend(true);
    HalvingSolver directSolver;
    lamopt::DriverOptions directOptions = partialOptions;
    directOptions.maxOuterIterations = 4;
    lamopt::GlobalOptimisationDriver directDriver(
        directBackend,
        approximationBuilder,
        directSolver,
        directOptions
    );

    lamopt::AnalysisRequest directRequest = request;
    directRequest.workDirectory = TempPath("driver_resume_direct");
    const lamopt::GlobalOptimisationResult directResult = directDriver.optimise(directRequest);

    EXPECT_TRUE(resumedResult.design.isApprox(directResult.design, 1.0e-12));
    EXPECT_NEAR(resumedResult.analysis.objectives(0), directResult.analysis.objectives(0), 1.0e-12);
    EXPECT_NE(resumedResult.message.find("Resumed from checkpoint with 2 recorded iterations."),
              std::string::npos);
}

TEST(GlobalDriverTest, DriverRejectsNonImprovingCandidateAndIncreasesDamping) {
    QuadraticBackend backend(true);
    lamopt::LinearApproximationBuilder approximationBuilder;
    ScriptedSolver scriptedSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 4;
    options.maxSubIterations = 3;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, scriptedSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, 1.0);
    request.workDirectory = TempPath("driver_rejection");

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_GE(result.history.size(), 2U);
    EXPECT_FALSE(result.history[0].accepted);
    EXPECT_TRUE(result.history[1].accepted);
    EXPECT_GT(result.history[0].dampingFactorsAfter(0), result.history[0].dampingFactorsBefore(0));
}

TEST(GlobalDriverTest, CoreMinMaxAdapterSolvesExactLinearSubproblem) {
    lamopt::CoreMinMax2Var4RespSubproblemSolver subproblemSolver;

    Problem1Backend backend;
    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 1.0);
    request.requestSensitivities = true;
    const lamopt::AnalysisResult analysis = backend.evaluate(request);

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = request.designVariables;
    problem.objectiveValues = analysis.objectives;
    problem.constraintValues = analysis.constraints;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(4);
    problem.designDampingVector = Eigen::Vector2d::Ones();
    problem.lowerBounds = Eigen::Vector2d::Zero();
    problem.objectiveGradients = analysis.objectiveGradients;
    problem.constraintGradients = analysis.constraintGradients;

    const lamopt::SubproblemResult result = subproblemSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_NEAR(result.candidateDesign(0), 0.744299, 1.0e-4);
    EXPECT_NEAR(result.candidateDesign(1), 0.569662, 1.0e-4);
    EXPECT_LT(result.predictedObjectives(0), analysis.objectives(0));
    EXPECT_LE(result.predictedConstraints.maxCoeff(), 1.0e-6);
}

TEST(GlobalDriverTest, DriverCanUseCoreMinMaxAdapter) {
    Problem1Backend backend;
    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::CoreMinMax2Var4RespSubproblemSolver subproblemSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 25;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-4;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 1.0);
    request.workDirectory = TempPath("driver_core_minmax");
    request.lowerBounds = Eigen::Vector2d::Zero();

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.design(0), 0.7886, 5.0e-3);
    EXPECT_NEAR(result.design(1), 0.4084, 5.0e-3);
    EXPECT_LE(result.analysis.constraints.maxCoeff(), 1.0e-6);
    EXPECT_LT(result.analysis.objectives(0), 100.0 * (2.0 * std::sqrt(2.0) + 1.0));
}

TEST(GlobalDriverTest, LaminateAwareWrapperProjectsSectionBlocks) {
    LaminateScriptedSolver scriptedSolver;
    lamopt::LaminateAwareSubproblemSolver laminateAwareSolver(scriptedSolver);

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(5);
    problem.referenceDesign(4) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -1.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.responseDampingFactors = Eigen::VectorXd::Ones(1);
    problem.designDampingVector = Eigen::VectorXd::Ones(5);

    lamopt::LaminateSectionState laminateSection;
    laminateSection.variableOffset = 0;
    laminateSection.isBalanced = true;
    laminateSection.isSymmetric = true;
    laminateSection.sublaminateCount = 6;
    laminateSection.thicknessLowerBound = 0.5;
    laminateSection.thicknessUpperBound = 2.0;
    laminateSection.thickness = 1.0;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    problem.laminateSections.push_back(laminateSection);

    const lamopt::SubproblemResult result = laminateAwareSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_LT(result.candidateDesign(4), 2.0);
    EXPECT_GT(result.candidateDesign(4), 0.5);
    EXPECT_GT(result.candidateDesign(4), 1.5);
    EXPECT_LT(result.candidateDesign(4), 5.0);
    EXPECT_TRUE(result.candidateDesign.allFinite());
    EXPECT_NE(result.message.find("laminate-aware projection complete"), std::string::npos);
}

TEST(GlobalDriverTest, DriverCanRunWithLaminateAwareSubproblemWrapper) {
    class ThicknessBackend final : public lamopt::AnalysisBackend {
    public:
        lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
            lamopt::AnalysisResult result;
            result.status = lamopt::AnalysisStatus::Success;
            result.objectives = Eigen::VectorXd::Constant(1, -request.designVariables(4));
            result.constraints = Eigen::VectorXd();
            result.objectiveGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), 1);
            (*result.objectiveGradients)(4, 0) = -1.0;
            return result;
        }
    } backend;

    LaminateScriptedSolver scriptedSolver;
    lamopt::LaminateAwareSubproblemSolver laminateAwareSolver(scriptedSolver);
    lamopt::LinearApproximationBuilder approximationBuilder;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 4;
    options.maxSubIterations = 2;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, laminateAwareSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Zero(5);
    request.designVariables(4) = 1.0;
    request.workDirectory = TempPath("driver_laminate_wrapper");

    lamopt::LaminateSectionState laminateSection;
    laminateSection.variableOffset = 0;
    laminateSection.isBalanced = true;
    laminateSection.isSymmetric = true;
    laminateSection.sublaminateCount = 6;
    laminateSection.thicknessLowerBound = 0.5;
    laminateSection.thicknessUpperBound = 2.0;
    laminateSection.thickness = 1.0;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    request.laminateSections.push_back(laminateSection);

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.design(4), 2.0);
    EXPECT_GT(result.design(4), 0.5);
    EXPECT_GT(result.design(4), 1.5);
    EXPECT_LT(result.analysis.objectives(0), -1.5);
}

TEST(GlobalDriverTest, CoreLaminateAdapterSolvesSectionConstrainedObjective) {
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(5);
    problem.referenceDesign(4) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -1.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.objectiveGradients = Eigen::MatrixXd::Zero(5, 1);
    (*problem.objectiveGradients)(4, 0) = -1.0;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(1);
    problem.designDampingVector = Eigen::VectorXd::Ones(5);

    lamopt::LaminateSectionState laminateSection;
    laminateSection.variableOffset = 0;
    laminateSection.isBalanced = true;
    laminateSection.isSymmetric = true;
    laminateSection.sublaminateCount = 6;
    laminateSection.thicknessLowerBound = 0.5;
    laminateSection.thicknessUpperBound = 2.0;
    laminateSection.thickness = 1.0;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    problem.laminateSections.push_back(laminateSection);

    const lamopt::SubproblemResult result = subproblemSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_TRUE(result.candidateDesign.allFinite());
    EXPECT_GT(result.candidateDesign(4), 1.2);
    EXPECT_LT(result.candidateDesign(4), 2.0);
    EXPECT_LT(result.predictedObjectives(0), -1.2);
}

TEST(GlobalDriverTest, CoreLaminateAdapterSupportsResponseConstraints) {
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(5);
    problem.referenceDesign(4) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -1.0);
    problem.constraintValues = Eigen::VectorXd::Constant(1, -0.4);
    problem.objectiveGradients = Eigen::MatrixXd::Zero(5, 1);
    (*problem.objectiveGradients)(4, 0) = -1.0;
    problem.constraintGradients = Eigen::MatrixXd::Zero(5, 1);
    (*problem.constraintGradients)(4, 0) = 1.0;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(2);
    problem.designDampingVector = Eigen::VectorXd::Ones(5);
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(0));

    const lamopt::SubproblemResult result = subproblemSolver.solve(problem);

    ASSERT_TRUE(result.success);
    ASSERT_EQ(result.predictedConstraints.size(), 1);
    EXPECT_TRUE(result.candidateDesign.allFinite());
    EXPECT_GT(result.candidateDesign(4), 1.1);
    EXPECT_LT(result.candidateDesign(4), 1.45);
    EXPECT_LE(result.predictedConstraints(0), 1.0e-6);
}

TEST(GlobalDriverTest, CoreLaminateAdapterSupportsSectionOffsets) {
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(7);
    problem.referenceDesign(6) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -1.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.objectiveGradients = Eigen::MatrixXd::Zero(7, 1);
    (*problem.objectiveGradients)(6, 0) = -1.0;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(1);
    problem.designDampingVector = Eigen::VectorXd::Ones(7);
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(2));

    const lamopt::SubproblemResult result = subproblemSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_TRUE(result.candidateDesign.allFinite());
    EXPECT_GT(result.candidateDesign(6), 1.2);
    EXPECT_LT(result.candidateDesign(6), 2.0);
}

TEST(GlobalDriverTest, CoreLaminateAdapterSupportsMultipleSectionBlocks) {
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(10);
    problem.referenceDesign(4) = 1.0;
    problem.referenceDesign(9) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -2.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.objectiveGradients = Eigen::MatrixXd::Zero(10, 1);
    (*problem.objectiveGradients)(4, 0) = -1.0;
    (*problem.objectiveGradients)(9, 0) = -1.0;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(1);
    problem.designDampingVector = Eigen::VectorXd::Ones(10);
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(0));
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(5));

    const lamopt::SubproblemResult result = subproblemSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_TRUE(result.candidateDesign.allFinite());
    EXPECT_GT(result.candidateDesign(4), 1.15);
    EXPECT_GT(result.candidateDesign(9), 1.15);
    EXPECT_LT(result.candidateDesign(4), 2.0);
    EXPECT_LT(result.candidateDesign(9), 2.0);
}

TEST(GlobalDriverTest, DefaultLaminateSolverUsesDirectCoreRouteWhenSupported) {
    LaminateScriptedSolver scriptedFallbackSolver;
    lamopt::CoreLaminateSection1RespSubproblemSolver directLaminateSolver;
    lamopt::DefaultLaminateSubproblemSolver routedSolver(scriptedFallbackSolver, directLaminateSolver);

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(5);
    problem.referenceDesign(4) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(1, -1.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.objectiveGradients = Eigen::MatrixXd::Zero(5, 1);
    (*problem.objectiveGradients)(4, 0) = -1.0;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(1);
    problem.designDampingVector = Eigen::VectorXd::Ones(5);
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(0));

    const lamopt::SubproblemResult result = routedSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_GT(result.candidateDesign(4), 1.2);
    EXPECT_LT(result.candidateDesign(4), 1.5);
    EXPECT_NE(result.message.find("route: direct_core_laminate"), std::string::npos);
}

TEST(GlobalDriverTest, DefaultLaminateSolverFallsBackToProjectionRouteWhenUnsupported) {
    LaminateScriptedSolver scriptedFallbackSolver;
    lamopt::CoreLaminateSection1RespSubproblemSolver directLaminateSolver;
    lamopt::DefaultLaminateSubproblemSolver routedSolver(scriptedFallbackSolver, directLaminateSolver);

    lamopt::ApproximationProblem problem;
    problem.referenceDesign = Eigen::VectorXd::Zero(5);
    problem.referenceDesign(4) = 1.0;
    problem.objectiveValues = Eigen::VectorXd::Constant(2, -1.0);
    problem.constraintValues = Eigen::VectorXd();
    problem.objectiveGradients = Eigen::MatrixXd::Zero(5, 2);
    (*problem.objectiveGradients)(4, 0) = -1.0;
    (*problem.objectiveGradients)(4, 1) = -0.5;
    problem.responseDampingFactors = Eigen::VectorXd::Ones(2);
    problem.designDampingVector = Eigen::VectorXd::Ones(5);
    problem.laminateSections.push_back(MakeBalancedSymmetricSection(0));

    const lamopt::SubproblemResult result = routedSolver.solve(problem);

    ASSERT_TRUE(result.success);
    EXPECT_GT(result.candidateDesign(4), 1.5);
    EXPECT_LT(result.candidateDesign(4), 2.0);
    EXPECT_NE(result.message.find("route: laminate_projection_fallback"), std::string::npos);
}

TEST(GlobalDriverTest, DriverImprovesWithCoreLaminateSectionAdapter) {
    class ThicknessBackend final : public lamopt::AnalysisBackend {
    public:
        lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
            lamopt::AnalysisResult result;
            result.status = lamopt::AnalysisStatus::Success;
            result.objectives = Eigen::VectorXd::Constant(1, -request.designVariables(4));
            result.constraints = Eigen::VectorXd();
            result.objectiveGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), 1);
            (*result.objectiveGradients)(4, 0) = -1.0;
            return result;
        }
    } backend;

    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Zero(5);
    request.designVariables(4) = 1.0;
    request.workDirectory = TempPath("driver_core_laminate_adapter");

    lamopt::LaminateSectionState laminateSection;
    laminateSection.variableOffset = 0;
    laminateSection.isBalanced = true;
    laminateSection.isSymmetric = true;
    laminateSection.sublaminateCount = 6;
    laminateSection.thicknessLowerBound = 0.5;
    laminateSection.thicknessUpperBound = 2.0;
    laminateSection.thickness = 1.0;
    laminateSection.laminationParameters = Eigen::VectorXd::Zero(4);
    request.laminateSections.push_back(laminateSection);

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_FALSE(result.history.empty());
    EXPECT_TRUE(std::any_of(result.history.begin(),
                            result.history.end(),
                            [](const lamopt::IterationRecord& record) { return record.accepted; }));
    EXPECT_GT(result.design(4), 1.5);
    EXPECT_LT(result.design(4), 2.0);
    EXPECT_LT(result.analysis.objectives(0), -1.5);
}

TEST(GlobalDriverTest, DriverImprovesWithOffsetCoreLaminateSectionAdapter) {
    class ThicknessBackend final : public lamopt::AnalysisBackend {
    public:
        lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
            lamopt::AnalysisResult result;
            result.status = lamopt::AnalysisStatus::Success;
            result.objectives = Eigen::VectorXd::Constant(1, -request.designVariables(6));
            result.constraints = Eigen::VectorXd();
            result.objectiveGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), 1);
            (*result.objectiveGradients)(6, 0) = -1.0;
            return result;
        }
    } backend;

    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::CoreLaminateSection1RespSubproblemSolver subproblemSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Zero(7);
    request.designVariables(6) = 1.0;
    request.workDirectory = TempPath("driver_core_laminate_adapter_offset");
    request.laminateSections.push_back(MakeBalancedSymmetricSection(2));

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_FALSE(result.history.empty());
    EXPECT_TRUE(std::any_of(result.history.begin(),
                            result.history.end(),
                            [](const lamopt::IterationRecord& record) { return record.accepted; }));
    EXPECT_GT(result.design(6), 1.5);
    EXPECT_LT(result.design(6), 2.0);
    EXPECT_LT(result.analysis.objectives(0), -1.5);
}

TEST(GlobalDriverTest, DriverUsesDirectLaminateRouteByDefault) {
    class ThicknessBackend final : public lamopt::AnalysisBackend {
    public:
        lamopt::AnalysisResult evaluate(const lamopt::AnalysisRequest& request) override {
            lamopt::AnalysisResult result;
            result.status = lamopt::AnalysisStatus::Success;
            result.objectives = Eigen::VectorXd::Constant(1, -request.designVariables(4));
            result.constraints = Eigen::VectorXd();
            result.objectiveGradients = Eigen::MatrixXd::Zero(request.designVariables.size(), 1);
            (*result.objectiveGradients)(4, 0) = -1.0;
            return result;
        }
    } backend;

    LaminateScriptedSolver scriptedFallbackSolver;
    lamopt::CoreLaminateSection1RespSubproblemSolver directLaminateSolver;
    lamopt::DefaultLaminateSubproblemSolver routedSolver(scriptedFallbackSolver, directLaminateSolver);
    lamopt::LinearApproximationBuilder approximationBuilder;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 12;
    options.maxSubIterations = 4;
    options.requestSensitivities = true;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, routedSolver, options);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Zero(5);
    request.designVariables(4) = 1.0;
    request.workDirectory = TempPath("driver_default_laminate_route");
    request.laminateSections.push_back(MakeBalancedSymmetricSection(0));

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    EXPECT_FALSE(result.history.empty());
    EXPECT_NE(result.history.front().message.find("route: direct_core_laminate"), std::string::npos);
    EXPECT_GT(result.design(4), 1.5);
    EXPECT_LT(result.design(4), 2.0);
    EXPECT_LT(result.analysis.objectives(0), -1.5);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
