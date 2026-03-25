#include "../src/OptimisationPipeline/calculixJobBackend.hpp"
#include "../src/OptimisationPipeline/laminateFeParameterisation.hpp"
#include "../src/OptimisationPipeline/laminateParameterisedBackend.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>

namespace {

std::filesystem::path TempPath(const std::string& name) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("laminateoptimiser_" + name);
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
    return path;
}

lamopt::LaminateSectionState MakeSection(const Eigen::Index variableOffset,
                                         const bool isBalanced,
                                         const bool isSymmetric,
                                         const int sublaminateCount = 6) {
    lamopt::LaminateSectionState section;
    section.variableOffset = variableOffset;
    section.isBalanced = isBalanced;
    section.isSymmetric = isSymmetric;
    section.sublaminateCount = sublaminateCount;
    section.thicknessLowerBound = 0.5;
    section.thicknessUpperBound = 2.0;
    section.thickness = 1.0;
    return section;
}

std::filesystem::path WriteTemplate(const std::filesystem::path& directory) {
    const std::filesystem::path templatePath = directory / "template.inp";
    std::ofstream output(templatePath);
    output << "*Heading\n";
    output << "** scalar={{X0}}\n";
    output << "** thickness={{SEC0_THICKNESS}}\n";
    output << "** a_lp0={{SEC0_A_LP0}}\n";
    output << "** d_lp1={{SEC0_D_LP1}}\n";
    output << "** offset={{SEC0_OFFSET}}\n";
    return templatePath;
}

}  // namespace

TEST(LaminateFeParameterisationTest, CanonicalLayoutForBalancedSymmetricSectionUsesADThenThickness) {
    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::LinSpaced(8, 0.1, 0.8);
    request.laminateSections.push_back(MakeSection(1, true, true));

    std::vector<lamopt::CanonicalLaminateSectionLayout> layouts;
    std::string message;
    ASSERT_TRUE(lamopt::BuildCanonicalLaminateSectionLayouts(request, layouts, message)) << message;
    ASSERT_EQ(layouts.size(), 1U);

    const auto& layout = layouts.front();
    EXPECT_EQ(layout.variableOffset, 1);
    EXPECT_EQ(layout.sectionSize, 5);
    EXPECT_EQ(layout.laminationParameterCount, 2);
    EXPECT_EQ(layout.thicknessIndex, 5);
    ASSERT_EQ(layout.laminationParameterSlots.size(), 4U);

    EXPECT_EQ(layout.laminationParameterSlots[0].globalIndex, 1);
    EXPECT_EQ(layout.laminationParameterSlots[0].blockKind, lamopt::LaminateTensorBlockKind::A);
    EXPECT_EQ(layout.laminationParameterSlots[0].laminationParameterIndex, 0);
    EXPECT_EQ(layout.laminationParameterSlots[1].blockKind, lamopt::LaminateTensorBlockKind::A);
    EXPECT_EQ(layout.laminationParameterSlots[1].laminationParameterIndex, 1);
    EXPECT_EQ(layout.laminationParameterSlots[2].blockKind, lamopt::LaminateTensorBlockKind::D);
    EXPECT_EQ(layout.laminationParameterSlots[2].laminationParameterIndex, 0);
    EXPECT_EQ(layout.laminationParameterSlots[3].blockKind, lamopt::LaminateTensorBlockKind::D);
    EXPECT_EQ(layout.laminationParameterSlots[3].laminationParameterIndex, 1);
}

TEST(LaminateFeParameterisationTest, CanonicalLayoutForGeneralSectionUsesADBThenThickness) {
    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::LinSpaced(20, 1.0, 20.0);
    request.laminateSections.push_back(MakeSection(2, false, false));

    std::vector<lamopt::CanonicalLaminateSectionLayout> layouts;
    std::string message;
    ASSERT_TRUE(lamopt::BuildCanonicalLaminateSectionLayouts(request, layouts, message)) << message;
    ASSERT_EQ(layouts.size(), 1U);

    const auto& layout = layouts.front();
    EXPECT_EQ(layout.sectionSize, 13);
    EXPECT_EQ(layout.laminationParameterCount, 4);
    EXPECT_EQ(layout.thicknessIndex, 14);
    ASSERT_EQ(layout.laminationParameterSlots.size(), 12U);

    for (int index = 0; index < 4; ++index) {
        EXPECT_EQ(layout.laminationParameterSlots[static_cast<std::size_t>(index)].blockKind,
                  lamopt::LaminateTensorBlockKind::A);
    }
    for (int index = 4; index < 8; ++index) {
        EXPECT_EQ(layout.laminationParameterSlots[static_cast<std::size_t>(index)].blockKind,
                  lamopt::LaminateTensorBlockKind::D);
    }
    for (int index = 8; index < 12; ++index) {
        EXPECT_EQ(layout.laminationParameterSlots[static_cast<std::size_t>(index)].blockKind,
                  lamopt::LaminateTensorBlockKind::B);
    }
}

TEST(LaminateFeParameterisationTest, JobBackendRendersResolvedLaminateTemplateParameters) {
    const std::filesystem::path tempDirectory = TempPath("laminate_fe_template_render");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.renderedInputFilename = "job.inp";
    config.resultFilename = "analysis_results.txt";
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}};
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend backend(config);

    lamopt::AnalysisRequest request;
    request.designVariables = Eigen::VectorXd::LinSpaced(6, 1.0, 6.0);
    request.workDirectory = tempDirectory / "run";
    request.laminateSections.push_back(MakeSection(1, true, true));
    request.templateParameters = lamopt::MakeLaminateTemplateParameters(request);

    const lamopt::AnalysisResult result = backend.evaluate(request);
    ASSERT_TRUE(result.isSuccessful()) << result.diagnostics.message;

    std::ifstream renderedInput(request.workDirectory / "job.inp");
    std::stringstream buffer;
    buffer << renderedInput.rdbuf();
    const std::string contents = buffer.str();

    EXPECT_NE(contents.find("scalar=1.0000000000000000e+00"), std::string::npos);
    EXPECT_NE(contents.find("thickness=6.0000000000000000e+00"), std::string::npos);
    EXPECT_NE(contents.find("a_lp0=2.0000000000000000e+00"), std::string::npos);
    EXPECT_NE(contents.find("d_lp1=5.0000000000000000e+00"), std::string::npos);
    EXPECT_NE(contents.find("offset=1"), std::string::npos);
}

TEST(LaminateFeParameterisationTest, ParameterisedBackendRebuildsLaminateTokensForEachEvaluation) {
    const std::filesystem::path tempDirectory = TempPath("laminate_parameterised_backend");
    const std::filesystem::path templatePath = WriteTemplate(tempDirectory);
    const std::filesystem::path fixturePath =
        std::filesystem::path(__FILE__).parent_path() / "fixtures/calculix/results_success.txt";

    lamopt::CalculixJobConfig config;
    config.templateInputPath = templatePath;
    config.renderedInputFilename = "job.inp";
    config.resultFilename = "analysis_results.txt";
    config.launchCommandTemplate = "cp " + fixturePath.string() + " {result_file}";
    config.parameterMappings = {{"{{X0}}", 0}};
    config.scratchRoot = tempDirectory;

    lamopt::CalculixJobBackend rawBackend(config);
    lamopt::LaminateParameterisedBackend backend(rawBackend);

    lamopt::AnalysisRequest request = {};
    request.designVariables = Eigen::VectorXd::LinSpaced(6, 1.0, 6.0);
    request.workDirectory = tempDirectory / "run_1";
    request.laminateSections.push_back(MakeSection(1, true, true));

    ASSERT_TRUE(backend.evaluate(request).isSuccessful());

    std::ifstream renderedInput1(request.workDirectory / "job.inp");
    std::stringstream buffer1;
    buffer1 << renderedInput1.rdbuf();
    EXPECT_NE(buffer1.str().find("thickness=6.0000000000000000e+00"), std::string::npos);

    request.designVariables(5) = 1.75;
    request.workDirectory = tempDirectory / "run_2";

    ASSERT_TRUE(backend.evaluate(request).isSuccessful());

    std::ifstream renderedInput2(request.workDirectory / "job.inp");
    std::stringstream buffer2;
    buffer2 << renderedInput2.rdbuf();
    EXPECT_NE(buffer2.str().find("thickness=1.7500000000000000e+00"), std::string::npos);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
