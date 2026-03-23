#include "../src/OptimisationPipeline/calculixComposipyBenchmark.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <system_error>
#include <vector>

namespace {

std::filesystem::path TempPath(const std::string& name) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("laminateoptimiser_" + name);
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
}

bool CanLaunchCalculix(const std::filesystem::path& executablePath) {
    if (executablePath.empty()) {
        return false;
    }

    const std::string command =
        lamopt::ShellQuoteForCommand(executablePath) + " -v >/dev/null 2>&1";
    const int exitCode = std::system(command.c_str());
    return exitCode == 0 || exitCode == 201;
}

bool IsCalculixCandidate(const std::filesystem::path& path) {
    const std::string filename = path.filename().string();
    return filename == "ccx"
        || filename.rfind("ccx_", 0) == 0
        || filename.rfind("ccx-", 0) == 0;
}

void CollectCalculixCandidates(const std::filesystem::path& root,
                               std::vector<std::filesystem::path>& candidates) {
    std::error_code errorCode;
    if (!std::filesystem::exists(root, errorCode)) {
        return;
    }

    if (std::filesystem::is_regular_file(root, errorCode)) {
        if (IsCalculixCandidate(root)) {
            candidates.push_back(root);
        }
        return;
    }

    for (std::filesystem::recursive_directory_iterator iterator(
             root,
             std::filesystem::directory_options::skip_permission_denied,
             errorCode);
         iterator != std::filesystem::recursive_directory_iterator();
         iterator.increment(errorCode)) {
        if (errorCode) {
            continue;
        }
        if (!iterator->is_regular_file(errorCode)) {
            continue;
        }
        if (IsCalculixCandidate(iterator->path())) {
            candidates.push_back(iterator->path());
        }
    }
}

std::optional<std::filesystem::path> DiscoverCalculixExecutable() {
    if (const char* configured = std::getenv("LAMOPT_CCX_EXECUTABLE")) {
        const std::filesystem::path path(configured);
        if (CanLaunchCalculix(path)) {
            return path;
        }
    }
    if (const char* configured = std::getenv("CCX")) {
        const std::filesystem::path path(configured);
        if (CanLaunchCalculix(path)) {
            return path;
        }
    }

    const std::filesystem::path resolved = lamopt::ResolveCalculixExecutable();
    if (CanLaunchCalculix(resolved)) {
        return resolved;
    }

    std::vector<std::filesystem::path> candidates;
    CollectCalculixCandidates("/opt/homebrew/bin", candidates);
    CollectCalculixCandidates("/usr/local/bin", candidates);
    CollectCalculixCandidates("/opt/homebrew/Cellar/calculix-ccx", candidates);
    CollectCalculixCandidates("/usr/local/Cellar/calculix-ccx", candidates);

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    for (auto iterator = candidates.rbegin(); iterator != candidates.rend(); ++iterator) {
        if (CanLaunchCalculix(*iterator)) {
            return *iterator;
        }
    }

    return std::nullopt;
}

}  // namespace

TEST(CalculixComposipyBenchmarkTest, DriverOptimisesRealCompositePlateBenchmark) {
    const std::optional<std::filesystem::path> executablePath = DiscoverCalculixExecutable();
    if (!executablePath.has_value()) {
        GTEST_SKIP() << "CalculiX executable was not found. Set LAMOPT_CCX_EXECUTABLE or CCX to enable this test.";
    }

    const std::filesystem::path tempDirectory = TempPath("calculix_composipy_benchmark");
    const std::filesystem::path templatePath = tempDirectory / "composipy_benchmark_template.inp";
    const lamopt::CalculixComposipyBenchmarkSpec spec;
    lamopt::WriteCalculixComposipyBenchmarkTemplate(templatePath, spec);

    lamopt::CalculixComposipyBenchmarkBackend backend(
        lamopt::MakeCalculixComposipyBenchmarkJobSetup(templatePath, tempDirectory, executablePath),
        spec
    );

    lamopt::AnalysisRequest baselineRequest;
    baselineRequest.designVariables = Eigen::VectorXd::Constant(1, 0.30);
    baselineRequest.workDirectory = tempDirectory / "baseline";
    baselineRequest.lowerBounds = Eigen::VectorXd::Constant(1, 0.12);
    baselineRequest.upperBounds = Eigen::VectorXd::Constant(1, 0.30);
    baselineRequest.requestSensitivities = true;

    const lamopt::AnalysisResult baselineResult = backend.evaluate(baselineRequest);
    ASSERT_TRUE(baselineResult.isSuccessful()) << baselineResult.diagnostics.message;
    ASSERT_EQ(baselineResult.objectives.size(), 1);
    ASSERT_EQ(baselineResult.constraints.size(), 2);
    ASSERT_TRUE(baselineResult.sensitivityPolicy.has_value());
    ASSERT_EQ(baselineResult.sensitivityPolicy->objectiveSources.size(), 1U);
    ASSERT_EQ(baselineResult.sensitivityPolicy->constraintSources.size(), 2U);
    EXPECT_EQ(baselineResult.sensitivityPolicy->objectiveSources[0],
              lamopt::ResponseSensitivitySource::BackendNative);
    EXPECT_EQ(baselineResult.sensitivityPolicy->constraintSources[0],
              lamopt::ResponseSensitivitySource::FiniteDifference);
    EXPECT_EQ(baselineResult.sensitivityPolicy->constraintSources[1],
              lamopt::ResponseSensitivitySource::FiniteDifference);
    EXPECT_LT(baselineResult.constraints.maxCoeff(), 0.0);

    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::GradientPenaltySubproblemSolver subproblemSolver({0.5, 10.0, 1.0e-12});

    lamopt::DriverOptions options;
    options.maxOuterIterations = 8;
    options.maxSubIterations = 3;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;
    options.finiteDifferenceStep = 1.0e-3;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request = baselineRequest;
    request.workDirectory = tempDirectory / "driver";

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.analysis.isSuccessful()) << result.message << " " << result.analysis.diagnostics.message;
    ASSERT_FALSE(result.history.empty());
    EXPECT_TRUE(result.analysis.hasAllGradients());
    EXPECT_LT(result.design(0), baselineRequest.designVariables(0));
    EXPECT_GT(result.design(0), 0.18);
    EXPECT_LT(result.design(0), 0.24);
    EXPECT_LE(result.analysis.constraints.maxCoeff(), 5.0e-3);
    EXPECT_LT(result.analysis.objectives(0), baselineResult.objectives(0));
    EXPECT_NE(result.analysis.diagnostics.message.find("Composipy benchmark responses assembled from real CalculiX output."),
              std::string::npos);
    EXPECT_NE(result.analysis.diagnostics.message.find("Sensitivity policy: objective[0]=backend_native; constraint[0]=finite_difference; constraint[1]=finite_difference."),
              std::string::npos);

    EXPECT_TRUE(std::filesystem::exists(request.workDirectory / "job.inp"));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.runDirectory / "job.dat"));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardOutputPath));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardErrorPath));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
