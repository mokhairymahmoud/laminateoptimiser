#include "../src/OptimisationPipeline/calculixJobBackend.hpp"
#include "../src/OptimisationPipeline/globalOptimisationDriver.hpp"

#include <Eigen/Dense>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::filesystem::path PrepareExampleDirectory() {
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / "lamopt_mock_calculix_e2e";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    return directory;
}

std::filesystem::path WriteTemplate(const std::filesystem::path& directory) {
    const std::filesystem::path templatePath = directory / "plate_template.inp";
    std::ofstream output(templatePath);
    if (!output) {
        throw std::runtime_error("Unable to write example template: " + templatePath.string());
    }

    output << "*HEADING\n";
    output << "Mock plate optimisation example\n";
    output << "** THICKNESS={{THICKNESS}}\n";
    output << "** Replace this file with your real CalculiX input deck.\n";
    return templatePath;
}

std::filesystem::path WriteMockSolverScript(const std::filesystem::path& directory) {
    const std::filesystem::path scriptPath = directory / "fake_ccx.sh";
    std::ofstream output(scriptPath);
    if (!output) {
        throw std::runtime_error("Unable to write mock solver script: " + scriptPath.string());
    }

    output << "#!/bin/sh\n";
    output << "set -eu\n";
    output << "thickness=$(awk -F= '/THICKNESS=/{print $2}' job.inp | tr -d '[:space:]')\n";
    output << "awk -v t=\"$thickness\" 'BEGIN {\n";
    output << "  mass = t;\n";
    output << "  tip = 0.18 / t;\n";
    output << "  printf(\"TOTAL MASS = %.10f\\n\", mass);\n";
    output << "  printf(\"TIP DISPLACEMENT = %.10f\\n\", tip);\n";
    output << "}' > job.dat\n";
    output << "printf 'mock solver completed for t=%s\\n' \"$thickness\"\n";

    output.close();
    std::filesystem::permissions(
        scriptPath,
        std::filesystem::perms::owner_exec | std::filesystem::perms::owner_read
            | std::filesystem::perms::owner_write | std::filesystem::perms::group_exec
            | std::filesystem::perms::group_read | std::filesystem::perms::others_exec
            | std::filesystem::perms::others_read,
        std::filesystem::perm_options::replace
    );
    return scriptPath;
}

void PrintVector(const std::string& label, const Eigen::VectorXd& values) {
    std::cout << label << ": [";
    for (Eigen::Index index = 0; index < values.size(); ++index) {
        if (index != 0) {
            std::cout << ", ";
        }
        std::cout << values(index);
    }
    std::cout << "]\n";
}

}  // namespace

int main() {
    try {
        const std::filesystem::path exampleDirectory = PrepareExampleDirectory();
        const std::filesystem::path templatePath = WriteTemplate(exampleDirectory);
        const std::filesystem::path mockSolverPath = WriteMockSolverScript(exampleDirectory);

        lamopt::CalculixJobConfig backendConfig;
        backendConfig.templateInputPath = templatePath;
        backendConfig.renderedInputFilename = "job.inp";
        backendConfig.resultFilename = "analysis_results.txt";
        backendConfig.launchCommandTemplate =
            lamopt::ShellQuoteForCommand(mockSolverPath) + " {job_name}";
        backendConfig.launchInRunDirectory = true;
        backendConfig.parameterMappings = {{"{{THICKNESS}}", 0}};
        backendConfig.scratchRoot = exampleDirectory / "runs";
        backendConfig.standardOutputFilename = std::filesystem::path("ccx.stdout.log");
        backendConfig.standardErrorFilename = std::filesystem::path("ccx.stderr.log");
        backendConfig.rawScalarExtractions = {
            {
                "mass",
                {
                    std::filesystem::path("{job_name}.dat"),
                    "TOTAL MASS\\s*=\\s*([-+0-9.Ee]+)",
                    1,
                    lamopt::CalculixMatchSelection::Last,
                    1.0,
                    0.0,
                    "mass"
                }
            },
            {
                "tip_displacement",
                {
                    std::filesystem::path("{job_name}.dat"),
                    "TIP DISPLACEMENT\\s*=\\s*([-+0-9.Ee]+)",
                    1,
                    lamopt::CalculixMatchSelection::Last,
                    1.0,
                    0.0,
                    "tip_displacement"
                }
            }
        };
        backendConfig.objectiveResponses = {
            {"mass", lamopt::DerivedScalarResponseTransform::Identity, 1.0, 0.0, "mass_objective"}
        };
        backendConfig.constraintResponses = {
            {
                "tip_displacement",
                lamopt::DerivedScalarResponseTransform::Affine,
                1.0,
                -0.50,
                "tip_displacement_minus_limit"
            }
        };

        lamopt::CalculixJobBackend backend(backendConfig);
        lamopt::LinearApproximationBuilder approximationBuilder;
        lamopt::GradientPenaltySubproblemSolver subproblemSolver({0.08, 10.0, 1.0e-12});

        lamopt::DriverOptions options;
        options.maxOuterIterations = 12;
        options.maxSubIterations = 4;
        options.requestSensitivities = true;
        options.gradientFallbackMode = lamopt::GradientFallbackMode::FiniteDifference;
        options.finiteDifferenceStep = 1.0e-3;
        options.stagnationTolerance = 5.0e-2;

        lamopt::GlobalOptimisationDriver driver(backend, approximationBuilder, subproblemSolver, options);

        lamopt::AnalysisRequest request;
        request.designVariables = Eigen::VectorXd::Constant(1, 0.40);
        request.lowerBounds = Eigen::VectorXd::Constant(1, 0.10);
        request.upperBounds = Eigen::VectorXd::Constant(1, 0.60);
        request.workDirectory = exampleDirectory / "optimisation_run";
        request.requestSensitivities = true;

        const lamopt::AnalysisResult baseline = backend.evaluate(request);
        if (!baseline.isSuccessful()) {
            std::cerr << "Baseline analysis failed: " << baseline.diagnostics.message << '\n';
            return EXIT_FAILURE;
        }

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Baseline design\n";
        PrintVector("  design", request.designVariables);
        PrintVector("  objectives", baseline.objectives);
        PrintVector("  constraints", baseline.constraints);
        std::cout << '\n';

        const lamopt::GlobalOptimisationResult result = driver.optimise(request);
        if (!result.analysis.isSuccessful()) {
            std::cerr << "Optimisation failed: " << result.message << '\n';
            std::cerr << "Analysis diagnostics: " << result.analysis.diagnostics.message << '\n';
            return EXIT_FAILURE;
        }

        const std::filesystem::path artifactDirectory = exampleDirectory / "artifacts";
        const lamopt::OptimisationLogPaths artifactPaths =
            lamopt::WriteOptimisationArtifacts(artifactDirectory, result);
        const std::filesystem::path checkpointPath = artifactDirectory / "driver.chk";
        lamopt::WriteCheckpoint(checkpointPath, result);

        std::cout << "Optimised design\n";
        PrintVector("  design", result.design);
        PrintVector("  objectives", result.analysis.objectives);
        PrintVector("  constraints", result.analysis.constraints);
        std::cout << "  converged: " << (result.converged ? "true" : "false") << '\n';
        std::cout << "  iterations: " << result.history.size() << '\n';
        std::cout << "  route message: " << result.message << '\n';
        std::cout << "  analysis diagnostics: " << result.analysis.diagnostics.message << '\n';
        std::cout << '\n';
        std::cout << "Generated files\n";
        std::cout << "  example directory: " << exampleDirectory << '\n';
        std::cout << "  final run directory: " << result.analysis.diagnostics.runDirectory << '\n';
        std::cout << "  summary json: " << artifactPaths.summaryJsonPath << '\n';
        std::cout << "  history csv: " << artifactPaths.historyCsvPath << '\n';
        std::cout << "  checkpoint: " << checkpointPath << '\n';
        std::cout << '\n';
        std::cout << "Interpretation\n";
        std::cout << "  objective = mass = thickness\n";
        std::cout << "  constraint = tip_displacement - 0.50 <= 0, with tip_displacement = 0.18 / thickness\n";
        std::cout << "  expected optimum is near thickness = 0.36\n";
        std::cout << '\n';
        std::cout << "To turn this into a real CalculiX run, keep the same driver setup and replace:\n";
        std::cout << "  1. plate_template.inp with your real deck\n";
        std::cout << "  2. fake_ccx.sh with the real CalculiX executable\n";
        std::cout << "  3. the regex extraction rules with the quantities you need from .dat/.frd/output files\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << "Example failed: " << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
