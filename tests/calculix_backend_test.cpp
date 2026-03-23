#include "../src/OptimisationPipeline/calculixJobBackend.hpp"
#include "../src/OptimisationPipeline/defaultAnalysisBackend.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

namespace {

std::filesystem::path TempPath(const std::string& name) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("laminateoptimiser_" + name);
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
}

std::filesystem::path WriteTemplate(const std::filesystem::path& directory) {
    const std::filesystem::path templatePath = directory / "template.inp";
    std::ofstream output(templatePath);
    output << "*Heading\n";
    output << "** X0={{X0}}\n";
    output << "** X1={{X1}}\n";
    return templatePath;
}

class ScopedEnvironmentVariable {
public:
    ScopedEnvironmentVariable(const char* name, const char* value)
        : m_name(name) {
        const char* current = std::getenv(name);
        if (current != nullptr) {
            m_hadValue = true;
            m_originalValue = current;
        }

        if (value != nullptr) {
            setenv(name, value, 1);
        } else {
            unsetenv(name);
        }
    }

    ~ScopedEnvironmentVariable() {
        if (m_hadValue) {
            setenv(m_name.c_str(), m_originalValue.c_str(), 1);
        } else {
            unsetenv(m_name.c_str());
        }
    }

private:
    std::string m_name;
    bool m_hadValue = false;
    std::string m_originalValue;
};

}  // namespace

TEST(CalculixBackendTest, BackendRendersParametersAndParsesResults) {
    const std::filesystem::path tempDirectory = TempPath("calculix_success");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.renderedInputFilename = "job.inp";
    config.resultFilename = "analysis_results.txt";
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}, {"{{X1}}", 1}};
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.5, -2.5);
    request.workDirectory = tempDirectory / "run";
    request.requestSensitivities = true;

    const lamopt::AnalysisResult result = backend.evaluate(request);

    ASSERT_TRUE(result.isSuccessful());
    EXPECT_EQ(result.diagnostics.backendName, "CalculixJobBackend");
    EXPECT_EQ(result.objectives.size(), 2);
    EXPECT_EQ(result.constraints.size(), 2);
    ASSERT_TRUE(result.objectiveGradients.has_value());
    ASSERT_TRUE(result.constraintGradients.has_value());
    EXPECT_NEAR((*result.objectiveGradients)(0, 0), 1.0, 1.0e-12);
    EXPECT_NEAR((*result.constraintGradients)(1, 1), -0.75, 1.0e-12);

    std::ifstream renderedInput(request.workDirectory / "job.inp");
    std::stringstream buffer;
    buffer << renderedInput.rdbuf();
    const std::string contents = buffer.str();
    EXPECT_NE(contents.find("1.5000000000000000e+00"), std::string::npos);
    EXPECT_NE(contents.find("-2.5000000000000000e+00"), std::string::npos);
}

TEST(CalculixBackendTest, DefaultConfigUsesResolvedExecutableAndLogCapture) {
    ScopedEnvironmentVariable lamoptExecutable("LAMOPT_CCX_EXECUTABLE", "/opt/calculix/bin/ccx-custom");
    ScopedEnvironmentVariable ccxExecutable("CCX", "/opt/calculix/bin/ccx-env");

    const std::filesystem::path tempDirectory = TempPath("calculix_default_config");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);

    const lamopt::CalculixJobConfig config = lamopt::MakeDefaultCalculixJobConfig({
        templatePath,
        tempDirectory,
        {{"{{X0}}", 0}, {"{{X1}}", 1}}
    });

    EXPECT_EQ(config.templateInputPath, templatePath);
    EXPECT_EQ(config.scratchRoot, tempDirectory);
    EXPECT_EQ(config.launchCommandTemplate, "'/opt/calculix/bin/ccx-custom' -i {job_name}");
    EXPECT_TRUE(config.launchInRunDirectory);
    ASSERT_TRUE(config.standardOutputFilename.has_value());
    ASSERT_TRUE(config.standardErrorFilename.has_value());
    EXPECT_EQ(*config.standardOutputFilename, std::filesystem::path("ccx.stdout.log"));
    EXPECT_EQ(*config.standardErrorFilename, std::filesystem::path("ccx.stderr.log"));
}

TEST(CalculixBackendTest, BackendUsesCalculixDefaultRunDirectoryWhenUnset) {
    const std::filesystem::path tempDirectory = TempPath("calculix_default_run_dir");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}, {"{{X1}}", 1}};
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 2.0);

    const lamopt::AnalysisResult result = backend.evaluate(request);

    ASSERT_TRUE(result.isSuccessful());
    EXPECT_EQ(result.diagnostics.runDirectory, tempDirectory / "calculix_job");
    EXPECT_TRUE(std::filesystem::exists(tempDirectory / "calculix_job" / "job.inp"));
}

TEST(CalculixBackendTest, BackendRunsInsideRunDirectoryAndCapturesLogs) {
    const std::filesystem::path tempDirectory = TempPath("calculix_logs");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate =
        "pwd; cat missing.stderr.probe; cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}, {"{{X1}}", 1}};
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 2.0);

    const lamopt::AnalysisResult result = backend.evaluate(request);

    ASSERT_TRUE(result.isSuccessful());
    ASSERT_FALSE(result.diagnostics.standardOutputPath.empty());
    ASSERT_FALSE(result.diagnostics.standardErrorPath.empty());

    std::ifstream stdoutFile(result.diagnostics.standardOutputPath);
    std::stringstream stdoutBuffer;
    stdoutBuffer << stdoutFile.rdbuf();
    EXPECT_NE(stdoutBuffer.str().find((tempDirectory / "calculix_job").string()), std::string::npos);

    std::ifstream stderrFile(result.diagnostics.standardErrorPath);
    std::stringstream stderrBuffer;
    stderrBuffer << stderrFile.rdbuf();
    EXPECT_NE(stderrBuffer.str().find("missing.stderr.probe"), std::string::npos);
}

TEST(CalculixBackendTest, DefaultFactoryBuildsCalculixBackend) {
    const std::filesystem::path tempDirectory = TempPath("calculix_default_factory");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);

    const auto backend = lamopt::MakeDefaultAnalysisBackend({
        templatePath,
        tempDirectory,
        {{"{{X0}}", 0}, {"{{X1}}", 1}},
        "job.inp",
        "analysis_results.txt",
        std::filesystem::path("/bin/sh"),
        std::filesystem::path("ccx.stdout.log"),
        std::filesystem::path("ccx.stderr.log")
    });

    ASSERT_TRUE(backend != nullptr);
    EXPECT_NE(dynamic_cast<lamopt::CalculixJobBackend*>(backend.get()), nullptr);
}

TEST(CalculixBackendTest, ConfiguredFactoryCanUseAbaqusCompatibilityBackend) {
    const std::filesystem::path tempDirectory = TempPath("abaqus_compat_factory");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::AbaqusJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}, {"{{X1}}", 1}};
    config.scratchRoot = tempDirectory;

    const auto backend = lamopt::MakeConfiguredAnalysisBackend({
        lamopt::AnalysisBackendKind::AbaqusCompatibility,
        std::nullopt,
        config
    });

    ASSERT_TRUE(backend != nullptr);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 2.0);
    request.workDirectory = tempDirectory / "run";

    const auto result = backend->evaluate(request);

    ASSERT_TRUE(result.isSuccessful());
    EXPECT_EQ(result.diagnostics.backendName, "AbaqusJobBackend");
}

TEST(CalculixBackendTest, BackendCanExtractResponsesFromCalculixTextOutput) {
    const std::filesystem::path tempDirectory = TempPath("calculix_text_extract");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_sample.dat";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {job_name}.dat";
    config.parameterMappings = {{"{{X0}}", 0}, {"{{X1}}", 1}};
    config.scratchRoot = tempDirectory;
    config.objectiveExtractions = {
        {std::filesystem::path("{job_name}.dat"), "TOTAL MASS\\s*=\\s*([-+0-9.Ee]+)", 1, lamopt::CalculixMatchSelection::Last, 1.0, 0.0, "mass"}
    };
    config.constraintExtractions = {
        {std::filesystem::path("{job_name}.dat"), "FAILURE INDEX\\s*=\\s*([-+0-9.Ee]+)", 1, lamopt::CalculixMatchSelection::Last, 1.0, -1.0, "failure_margin"},
        {std::filesystem::path("{job_name}.dat"), "TIP DISPLACEMENT\\s*=\\s*([-+0-9.Ee]+)", 1, lamopt::CalculixMatchSelection::Last, 1000.0, 0.0, "tip_displacement_mm"}
    };

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::Vector2d(1.0, 2.0);

    const lamopt::AnalysisResult result = backend.evaluate(request);

    ASSERT_TRUE(result.isSuccessful());
    ASSERT_EQ(result.objectives.size(), 1);
    ASSERT_EQ(result.constraints.size(), 2);
    EXPECT_NEAR(result.objectives(0), 12.75, 1.0e-12);
    EXPECT_NEAR(result.constraints(0), -0.08, 1.0e-12);
    EXPECT_NEAR(result.constraints(1), 3.1, 1.0e-12);
    EXPECT_EQ(result.diagnostics.message, "parsed from CalculiX output files");
    EXPECT_TRUE(std::filesystem::exists(tempDirectory / "calculix_job" / "analysis_results.txt"));
}

TEST(CalculixBackendTest, BackendReportsCommandFailure) {
    const std::filesystem::path tempDirectory = TempPath("calculix_failure");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate = "exit 3";
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, 1.0);
    request.workDirectory = tempDirectory / "run";

    const lamopt::AnalysisResult result = backend.evaluate(request);

    EXPECT_EQ(result.status, lamopt::AnalysisStatus::BackendFailure);
    EXPECT_EQ(result.diagnostics.exitCode, 3);
}

TEST(CalculixBackendTest, BackendReportsTimeout) {
    const std::filesystem::path tempDirectory = TempPath("calculix_timeout");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.launchCommandTemplate = "sleep 1";
    config.timeout = std::chrono::milliseconds(100);
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::Constant(1, 1.0);
    request.workDirectory = tempDirectory / "run";

    const lamopt::AnalysisResult result = backend.evaluate(request);

    EXPECT_EQ(result.status, lamopt::AnalysisStatus::Timeout);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
