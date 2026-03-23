#pragma once

#include "jobBackend.hpp"

#include <cstdlib>
#include <regex>

namespace lamopt {

using CalculixParameterMapping = JobParameterMapping;

enum class CalculixMatchSelection {
    First,
    Last
};

struct CalculixTextScalarRule {
    std::filesystem::path sourceFilename;
    std::string pattern;
    int captureGroup = 1;
    CalculixMatchSelection matchSelection = CalculixMatchSelection::Last;
    double scale = 1.0;
    double offset = 0.0;
    std::string label;
};

struct CalculixNamedTextScalarRule {
    std::string id;
    CalculixTextScalarRule extraction;
};

enum class CalculixDerivedResponseTransform {
    Identity,
    Affine,
    AbsoluteAffine,
    InverseAffine
};

struct CalculixDerivedResponseRule {
    std::string sourceId;
    CalculixDerivedResponseTransform transform = CalculixDerivedResponseTransform::Identity;
    double scale = 1.0;
    double offset = 0.0;
    std::string label;
};

struct CalculixJobConfig : JobBackendConfig {
    std::vector<CalculixNamedTextScalarRule> rawScalarExtractions;
    std::vector<CalculixDerivedResponseRule> objectiveResponses;
    std::vector<CalculixDerivedResponseRule> constraintResponses;
    std::vector<CalculixTextScalarRule> objectiveExtractions;
    std::vector<CalculixTextScalarRule> constraintExtractions;

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
    std::vector<CalculixNamedTextScalarRule> rawScalarExtractions;
    std::vector<CalculixDerivedResponseRule> objectiveResponses;
    std::vector<CalculixDerivedResponseRule> constraintResponses;
    std::vector<CalculixTextScalarRule> objectiveExtractions;
    std::vector<CalculixTextScalarRule> constraintExtractions;
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
    config.rawScalarExtractions = setup.rawScalarExtractions;
    config.objectiveResponses = setup.objectiveResponses;
    config.constraintResponses = setup.constraintResponses;
    config.objectiveExtractions = setup.objectiveExtractions;
    config.constraintExtractions = setup.constraintExtractions;
    config.launchCommandTemplate =
        ShellQuoteForCommand(ResolveCalculixExecutable(setup.executablePath)) + " -i {job_name}";
    return config;
}

class CalculixJobBackend : public JobBackend {
public:
    explicit CalculixJobBackend(CalculixJobConfig config)
        : JobBackend(config)
        , m_config(std::move(config)) {}

protected:
    AnalysisResult parseSuccessfulRun(const ExecutionPaths& paths,
                                      const AnalysisDiagnostics& diagnostics) const override {
        if (usesResponseSchema()) {
            try {
                const std::unordered_map<std::string, double> rawValues =
                    extractNamedValues(m_config.rawScalarExtractions, paths);
                const Eigen::VectorXd objectives =
                    assembleDerivedValues(m_config.objectiveResponses, rawValues);
                const Eigen::VectorXd constraints =
                    assembleDerivedValues(m_config.constraintResponses, rawValues);
                writeResultFile(
                    paths.resultPath,
                    objectives,
                    constraints,
                    std::nullopt,
                    std::nullopt,
                    std::nullopt,
                    std::string("parsed from CalculiX response schema")
                );
            } catch (const std::exception& exception) {
                AnalysisResult result;
                result.status = AnalysisStatus::InvalidOutput;
                result.diagnostics = diagnostics;
                result.diagnostics.message = exception.what();
                return result;
            }

            return parseResultFile(paths.resultPath, diagnostics, m_config.backendName);
        }

        if (m_config.objectiveExtractions.empty() && m_config.constraintExtractions.empty()) {
            return JobBackend::parseSuccessfulRun(paths, diagnostics);
        }

        try {
            const Eigen::VectorXd objectives = extractValues(m_config.objectiveExtractions, paths);
            const Eigen::VectorXd constraints = extractValues(m_config.constraintExtractions, paths);
            writeResultFile(
                paths.resultPath,
                objectives,
                constraints,
                std::nullopt,
                std::nullopt,
                std::nullopt,
                std::string("parsed from CalculiX output files")
            );
        } catch (const std::exception& exception) {
            AnalysisResult result;
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics = diagnostics;
            result.diagnostics.message = exception.what();
            return result;
        }

        return parseResultFile(paths.resultPath, diagnostics, m_config.backendName);
    }

private:
    CalculixJobConfig m_config;

    [[nodiscard]] bool usesResponseSchema() const {
        return !m_config.rawScalarExtractions.empty()
            || !m_config.objectiveResponses.empty()
            || !m_config.constraintResponses.empty();
    }

    [[nodiscard]] static Eigen::VectorXd extractValues(const std::vector<CalculixTextScalarRule>& rules,
                                                       const ExecutionPaths& paths) {
        Eigen::VectorXd values(rules.size());
        for (std::size_t index = 0; index < rules.size(); ++index) {
            values(static_cast<Eigen::Index>(index)) = extractScalar(rules[index], paths);
        }
        return values;
    }

    [[nodiscard]] static std::unordered_map<std::string, double>
    extractNamedValues(const std::vector<CalculixNamedTextScalarRule>& rules,
                       const ExecutionPaths& paths) {
        std::unordered_map<std::string, double> values;
        values.reserve(rules.size());

        for (const CalculixNamedTextScalarRule& rule : rules) {
            if (rule.id.empty()) {
                throw std::runtime_error("CalculiX response schema extraction id must not be empty.");
            }
            if (values.find(rule.id) != values.end()) {
                throw std::runtime_error("CalculiX response schema extraction id is duplicated: " + rule.id);
            }

            values.emplace(rule.id, extractScalar(rule.extraction, paths));
        }

        return values;
    }

    [[nodiscard]] static Eigen::VectorXd
    assembleDerivedValues(const std::vector<CalculixDerivedResponseRule>& rules,
                          const std::unordered_map<std::string, double>& rawValues) {
        Eigen::VectorXd values(static_cast<Eigen::Index>(rules.size()));
        for (std::size_t index = 0; index < rules.size(); ++index) {
            values(static_cast<Eigen::Index>(index)) = evaluateDerivedRule(rules[index], rawValues);
        }
        return values;
    }

    [[nodiscard]] static double extractScalar(const CalculixTextScalarRule& rule,
                                              const ExecutionPaths& paths) {
        const std::filesystem::path sourcePath = resolveSourcePath(rule.sourceFilename, paths);
        std::ifstream sourceStream(sourcePath);
        if (!sourceStream) {
            throw std::runtime_error("CalculiX extraction source file not found: " + sourcePath.string());
        }

        std::ostringstream buffer;
        buffer << sourceStream.rdbuf();
        const std::string content = buffer.str();

        const std::regex expression(rule.pattern);
        std::sregex_iterator iterator(content.begin(), content.end(), expression);
        const std::sregex_iterator end;

        if (iterator == end) {
            throw std::runtime_error("CalculiX extraction pattern did not match for "
                                     + describeRule(rule, sourcePath) + ".");
        }

        std::smatch match = *iterator;
        if (rule.matchSelection == CalculixMatchSelection::Last) {
            for (; iterator != end; ++iterator) {
                match = *iterator;
            }
        }

        if (rule.captureGroup < 0 || rule.captureGroup >= static_cast<int>(match.size())) {
            throw std::runtime_error("CalculiX extraction capture group is out of range for "
                                     + describeRule(rule, sourcePath) + ".");
        }

        const double rawValue = std::stod(match[rule.captureGroup].str());
        return rule.scale * rawValue + rule.offset;
    }

    [[nodiscard]] static double
    evaluateDerivedRule(const CalculixDerivedResponseRule& rule,
                        const std::unordered_map<std::string, double>& rawValues) {
        const auto iterator = rawValues.find(rule.sourceId);
        if (iterator == rawValues.end()) {
            throw std::runtime_error("CalculiX response schema source id was not extracted: " + rule.sourceId);
        }

        const double value = iterator->second;
        switch (rule.transform) {
            case CalculixDerivedResponseTransform::Identity:
                return value;
            case CalculixDerivedResponseTransform::Affine:
                return rule.scale * value + rule.offset;
            case CalculixDerivedResponseTransform::AbsoluteAffine:
                return rule.scale * std::abs(value) + rule.offset;
            case CalculixDerivedResponseTransform::InverseAffine:
                if (std::abs(value) <= 1.0e-16) {
                    throw std::runtime_error("CalculiX response schema inverse transform encountered a zero source value for "
                                             + describeDerivedRule(rule) + ".");
                }
                return rule.scale / value + rule.offset;
        }

        throw std::runtime_error("Unsupported CalculiX response schema transform for "
                                 + describeDerivedRule(rule) + ".");
    }

    [[nodiscard]] static std::filesystem::path resolveSourcePath(const std::filesystem::path& sourceFilename,
                                                                 const ExecutionPaths& paths) {
        if (sourceFilename.is_absolute()) {
            return sourceFilename;
        }

        std::string resolved = sourceFilename.string();
        replaceToken(resolved, "{job_name}", paths.renderedInputPath.stem().string());
        replaceToken(resolved, "{run_dir}", paths.runDirectory.string());
        replaceToken(resolved, "{input_file}", paths.renderedInputPath.string());
        replaceToken(resolved, "{result_file}", paths.resultPath.string());
        return paths.runDirectory / resolved;
    }

    static void replaceToken(std::string& input,
                             const std::string& needle,
                             const std::string& replacement) {
        if (needle.empty()) {
            return;
        }

        std::size_t position = 0;
        while ((position = input.find(needle, position)) != std::string::npos) {
            input.replace(position, needle.size(), replacement);
            position += replacement.size();
        }
    }

    [[nodiscard]] static std::string describeRule(const CalculixTextScalarRule& rule,
                                                  const std::filesystem::path& sourcePath) {
        if (!rule.label.empty()) {
            return rule.label + " in " + sourcePath.string();
        }
        return sourcePath.string();
    }

    [[nodiscard]] static std::string describeDerivedRule(const CalculixDerivedResponseRule& rule) {
        if (!rule.label.empty()) {
            return rule.label;
        }
        return rule.sourceId;
    }
};

}  // namespace lamopt
