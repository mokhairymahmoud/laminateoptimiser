#include "../src/OptimisationPipeline/calculixLaminateParameterBenchmark.hpp"
#include "../src/OptimisationPipeline/configuredSubproblemSolver.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
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

lamopt::AnalysisRequest MakeInitialRequest(const std::filesystem::path& workDirectory) {
    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd(5);
    request.designVariables << 0.0, 0.0, 0.0, 0.0, 0.24;
    request.lowerBounds = Eigen::VectorXd(5);
    request.upperBounds = Eigen::VectorXd(5);
    *request.lowerBounds << -0.75, -0.75, -0.75, -0.75, 0.18;
    *request.upperBounds << 0.75, 0.75, 0.75, 0.75, 0.30;
    request.workDirectory = workDirectory;
    request.requestSensitivities = true;

    lamopt::LaminateSectionState section;
    section.variableOffset = 0;
    section.isBalanced = true;
    section.isSymmetric = true;
    section.sublaminateCount = 6;
    section.thicknessLowerBound = 0.18;
    section.thicknessUpperBound = 0.30;
    section.thickness = request.designVariables(4);
    section.laminationParameters = Eigen::VectorXd::Zero(4);
    request.laminateSections.push_back(section);
    return request;
}

}  // namespace

TEST(CalculixLaminateParameterBenchmarkTest, DriverUsesCanonicalLaminateTemplateParametersInRealRun) {
    const std::optional<std::filesystem::path> executablePath = DiscoverCalculixExecutable();
    if (!executablePath.has_value()) {
        GTEST_SKIP() << "CalculiX executable was not found. Set LAMOPT_CCX_EXECUTABLE or CCX to enable this test.";
    }

    const std::filesystem::path tempDirectory = TempPath("calculix_laminate_parameter_benchmark");
    const std::filesystem::path templatePath = tempDirectory / "laminate_parameter_benchmark.inp";

    lamopt::CalculixLaminateParameterBenchmarkSpec baselineSpec;
    baselineSpec.tipDisplacementLimit = 1.0;
    lamopt::WriteCalculixLaminateParameterBenchmarkTemplate(templatePath, baselineSpec);

    lamopt::CalculixLaminateParameterBenchmarkBackend baselineBackend(
        lamopt::MakeCalculixLaminateParameterBenchmarkJobSetup(templatePath, tempDirectory, executablePath),
        baselineSpec
    );

    const lamopt::AnalysisRequest baselineRequest = MakeInitialRequest(tempDirectory / "baseline");
    const lamopt::AnalysisResult baselineResult = baselineBackend.evaluate(baselineRequest);
    ASSERT_TRUE(baselineResult.isSuccessful()) << baselineResult.diagnostics.message;
    const double baselineTip = std::abs(baselineResult.extractedScalarValues.at("tip_u3"));
    ASSERT_GT(baselineTip, 0.0);

    lamopt::CalculixLaminateParameterBenchmarkSpec spec = baselineSpec;
    spec.tipDisplacementLimit = 0.95 * baselineTip;

    lamopt::CalculixLaminateParameterBenchmarkBackend backend(
        lamopt::MakeCalculixLaminateParameterBenchmarkJobSetup(templatePath, tempDirectory, executablePath),
        spec
    );
    auto derivativeProvider = lamopt::MakeCalculixLaminateParameterBenchmarkDerivativeProvider(
        MakeInitialRequest(tempDirectory / "driver"),
        spec
    );

    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::ConfiguredSubproblemSolver subproblemSolver;

    lamopt::DriverOptions options;
    options.maxOuterIterations = 10;
    options.maxSubIterations = 3;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;
    options.finiteDifferenceStep = 1.0e-3;
    options.stagnationTolerance = 1.0e-6;

    lamopt::AnalysisRequest request = MakeInitialRequest(tempDirectory / "driver");
    lamopt::GlobalOptimisationDriver driver(backend,
                                            approximationBuilder,
                                            subproblemSolver,
                                            options,
                                            &derivativeProvider);
    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.analysis.isSuccessful()) << result.message << " " << result.analysis.diagnostics.message;
    ASSERT_FALSE(result.history.empty());
    EXPECT_TRUE(result.analysis.hasAllGradients());
    ASSERT_TRUE(result.analysis.objectiveGradients.has_value());
    EXPECT_TRUE(std::any_of(result.history.begin(),
                            result.history.end(),
                            [](const lamopt::IterationRecord& record) { return record.accepted; }));
    EXPECT_NE(result.history.front().message.find("route: direct_core_laminate"), std::string::npos);
    EXPECT_LE(result.analysis.constraints.maxCoeff(), 5.0e-3);
    EXPECT_LT(result.design(4), request.designVariables(4));
    EXPECT_GT(result.design.head(4).norm(), 1.0e-3);
    EXPECT_NE(result.analysis.diagnostics.message.find("Canonical laminate template parameters and response-schema values assembled for real CalculiX output."),
              std::string::npos);
    EXPECT_NE(result.analysis.diagnostics.message.find("extracted FE quantities"), std::string::npos);
    EXPECT_EQ(result.analysis.sensitivityPolicy->objectiveSources[0],
              lamopt::ResponseSensitivitySource::OptimiserSideLaminate);
    EXPECT_EQ(result.analysis.sensitivityPolicy->constraintSources[0],
              lamopt::ResponseSensitivitySource::FiniteDifference);

    Eigen::MatrixXd expectedObjectiveGradient = Eigen::MatrixXd::Zero(request.designVariables.size(), 1);
    expectedObjectiveGradient(4, 0) = 1.0;
    EXPECT_TRUE(result.analysis.objectiveGradients->isApprox(expectedObjectiveGradient, 1.0e-12));

    std::ifstream renderedInput(request.workDirectory / "job.inp");
    std::stringstream renderedBuffer;
    renderedBuffer << renderedInput.rdbuf();
    const std::string renderedContents = renderedBuffer.str();
    EXPECT_EQ(renderedContents.find("{{SEC0_A_LP0}}"), std::string::npos);
    EXPECT_EQ(renderedContents.find("{{SEC0_THICKNESS}}"), std::string::npos);

    EXPECT_TRUE(std::filesystem::exists(request.workDirectory / "job.dat"));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardOutputPath));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardErrorPath));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
