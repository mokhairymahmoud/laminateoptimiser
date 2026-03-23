#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"

#include <gtest/gtest.h>

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

std::filesystem::path TempPath(const std::string& name) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("laminateoptimiser_" + name);
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
