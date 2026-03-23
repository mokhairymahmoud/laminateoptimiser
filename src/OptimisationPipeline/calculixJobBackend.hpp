#pragma once

#include "jobBackend.hpp"

#include <cstdlib>

namespace lamopt {

using CalculixParameterMapping = JobParameterMapping;

struct CalculixJobConfig : JobBackendConfig {
    CalculixJobConfig() {
        renderedInputFilename = "job.inp";
        resultFilename = "analysis_results.txt";
        launchCommandTemplate = "ccx -i {job_name}";
        launchInRunDirectory = true;
        standardOutputFilename = "ccx.stdout.log";
        standardErrorFilename = "ccx.stderr.log";
        backendName = "CalculixJobBackend";
        defaultRunDirectoryName = "calculix_job";
    }
};

struct CalculixJobSetup {
    std::filesystem::path templateInputPath;
    std::filesystem::path scratchRoot;
    std::vector<CalculixParameterMapping> parameterMappings;
    std::string renderedInputFilename = "job.inp";
    std::filesystem::path resultFilename = "analysis_results.txt";
    std::optional<std::filesystem::path> executablePath;
    std::optional<std::filesystem::path> standardOutputFilename = std::filesystem::path("ccx.stdout.log");
    std::optional<std::filesystem::path> standardErrorFilename = std::filesystem::path("ccx.stderr.log");
};

[[nodiscard]] inline std::filesystem::path ResolveCalculixExecutable(
    const std::optional<std::filesystem::path>& explicitExecutable = std::nullopt) {
    if (explicitExecutable.has_value()) {
        return *explicitExecutable;
    }

    if (const char* envExecutable = std::getenv("LAMOPT_CCX_EXECUTABLE")) {
        if (*envExecutable != '\0') {
            return std::filesystem::path(envExecutable);
        }
    }
    if (const char* envExecutable = std::getenv("CCX")) {
        if (*envExecutable != '\0') {
            return std::filesystem::path(envExecutable);
        }
    }

    return std::filesystem::path("ccx");
}

[[nodiscard]] inline std::string ShellQuoteForCommand(const std::filesystem::path& path) {
    std::string quoted = "'";
    const std::string value = path.string();
    for (const char character : value) {
        if (character == '\'') {
            quoted += "'\\''";
        } else {
            quoted += character;
        }
    }
    quoted += "'";
    return quoted;
}

[[nodiscard]] inline CalculixJobConfig MakeDefaultCalculixJobConfig(const CalculixJobSetup& setup) {
    CalculixJobConfig config;
    config.templateInputPath = setup.templateInputPath;
    config.scratchRoot = setup.scratchRoot;
    config.parameterMappings = setup.parameterMappings;
    config.renderedInputFilename = setup.renderedInputFilename;
    config.resultFilename = setup.resultFilename;
    config.standardOutputFilename = setup.standardOutputFilename;
    config.standardErrorFilename = setup.standardErrorFilename;
    config.launchCommandTemplate =
        ShellQuoteForCommand(ResolveCalculixExecutable(setup.executablePath)) + " -i {job_name}";
    return config;
}

class CalculixJobBackend : public JobBackend {
public:
    explicit CalculixJobBackend(CalculixJobConfig config)
        : JobBackend(std::move(config)) {}
};

}  // namespace lamopt
