#include "../src/OptimisationPipeline/defaultAnalysisBackend.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
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

std::filesystem::path WriteThicknessProbeTemplate(const std::filesystem::path& directory) {
    const std::filesystem::path templatePath = directory / "thickness_probe_template.inp";
    std::ofstream output(templatePath);
    output << "*HEADING\n";
    output << "LaminateOptimiser CalculiX end-to-end thickness probe\n";
    output << "*NODE, NSET=NALL\n";
    output << "1, 0., 0., 0.\n";
    output << "2, 1., 0., 0.\n";
    output << "3, 1., 1., 0.\n";
    output << "4, 0., 1., 0.\n";
    output << "*ELEMENT, TYPE=CPS4, ELSET=EALL\n";
    output << "1, 1, 2, 3, 4\n";
    output << "*MATERIAL, NAME=MAT\n";
    output << "*ELASTIC\n";
    output << "1000., 0.3\n";
    output << "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT\n";
    output << "{{THICKNESS}}\n";
    output << "*BOUNDARY\n";
    output << "1, 1, 2\n";
    output << "4, 1, 2\n";
    output << "2, 2, 2\n";
    output << "*STEP\n";
    output << "*STATIC\n";
    output << "*CLOAD\n";
    output << "2, 1, 1.0\n";
    output << "*NODE PRINT, NSET=NALL\n";
    output << "U\n";
    output << "*END STEP\n";
    return templatePath;
}

bool CanLaunchCalculix(const std::filesystem::path& executablePath) {
    if (executablePath.empty()) {
        return false;
    }

    const std::string command =
        lamopt::ShellQuoteForCommand(executablePath) + " -v >/dev/null 2>&1";
    return std::system(command.c_str()) == 0;
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

TEST(CalculixE2ETest, DriverImprovesPlateThicknessUsingRealCalculix) {
    const std::optional<std::filesystem::path> executablePath = DiscoverCalculixExecutable();
    if (!executablePath.has_value()) {
        GTEST_SKIP() << "CalculiX executable was not found. Set LAMOPT_CCX_EXECUTABLE or CCX to enable this test.";
    }

    const std::filesystem::path tempDirectory = TempPath("calculix_e2e_driver");
    const std::filesystem::path templatePath = WriteThicknessProbeTemplate(tempDirectory);

    const auto backend = lamopt::MakeDefaultAnalysisBackend({
        templatePath,
        tempDirectory,
        {{"{{THICKNESS}}", 0}},
        "job.inp",
        "analysis_results.txt",
        *executablePath,
        std::filesystem::path("ccx.stdout.log"),
        std::filesystem::path("ccx.stderr.log"),
        {},
        {},
        {},
        {
            {
                std::filesystem::path("{job_name}.dat"),
                "\\n\\s*2\\s+([-+0-9.Ee]+)\\s+[-+0-9.Ee]+\\s+[-+0-9.Ee]+",
                1,
                lamopt::CalculixMatchSelection::Last,
                1.0,
                0.0,
                "node_2_x_displacement"
            }
        },
        {}
    });

    ASSERT_TRUE(backend != nullptr);

    lamopt::AnalysisRequest baselineRequest;
    baselineRequest.designVariables = Eigen::VectorXd::Constant(1, 0.5);
    baselineRequest.workDirectory = tempDirectory / "baseline";
    baselineRequest.lowerBounds = Eigen::VectorXd::Constant(1, 0.5);
    baselineRequest.upperBounds = Eigen::VectorXd::Constant(1, 2.0);

    const lamopt::AnalysisResult baselineResult = backend->evaluate(baselineRequest);
    ASSERT_TRUE(baselineResult.isSuccessful()) << baselineResult.diagnostics.message;
    ASSERT_EQ(baselineResult.objectives.size(), 1);

    lamopt::LinearApproximationBuilder approximationBuilder;
    lamopt::GradientPenaltySubproblemSolver subproblemSolver({0.75, 5.0, 1.0e-12});

    lamopt::DriverOptions options;
    options.maxOuterIterations = 6;
    options.maxSubIterations = 2;
    options.requestSensitivities = true;
    options.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;
    options.finiteDifferenceStep = 1.0e-2;
    options.stagnationTolerance = 1.0e-6;

    lamopt::GlobalOptimisationDriver driver(*backend, approximationBuilder, subproblemSolver, options);

    lamopt::AnalysisRequest request = baselineRequest;
    request.workDirectory = tempDirectory / "driver";

    const lamopt::GlobalOptimisationResult result = driver.optimise(request);

    ASSERT_TRUE(result.analysis.isSuccessful()) << result.message << " " << result.analysis.diagnostics.message;
    ASSERT_FALSE(result.history.empty());
    EXPECT_TRUE(std::any_of(result.history.begin(),
                            result.history.end(),
                            [](const lamopt::IterationRecord& record) { return record.accepted; }));
    EXPECT_TRUE(result.analysis.hasAllGradients());
    EXPECT_GT(result.design(0), baselineRequest.designVariables(0));
    EXPECT_LE(result.design(0), 2.0 + 1.0e-12);
    EXPECT_GT(result.design(0), 1.5);
    EXPECT_LT(result.analysis.objectives(0), 0.5 * baselineResult.objectives(0));
    EXPECT_NE(result.analysis.diagnostics.message.find("Finite-difference gradients attached."),
              std::string::npos);

    EXPECT_TRUE(std::filesystem::exists(request.workDirectory / "job.inp"));
    EXPECT_TRUE(std::filesystem::exists(request.workDirectory / "job.dat"));
    EXPECT_TRUE(std::filesystem::exists(request.workDirectory / "fd_0" / "job.dat"));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.runDirectory / "job.dat"));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardOutputPath));
    EXPECT_TRUE(std::filesystem::exists(result.analysis.diagnostics.standardErrorPath));
}
