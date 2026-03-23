#include "../src/OptimisationPipeline/calculixJobBackend.hpp"

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

}  // namespace

TEST(CalculixBackendTest, BackendRendersParametersAndParsesResults) {
    const std::filesystem::path tempDirectory = TempPath("calculix_success");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/abaqus/results_success.txt";

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

TEST(CalculixBackendTest, BackendUsesCalculixDefaultRunDirectoryWhenUnset) {
    const std::filesystem::path tempDirectory = TempPath("calculix_default_run_dir");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/abaqus/results_success.txt";

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
