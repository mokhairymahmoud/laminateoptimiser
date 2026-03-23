#pragma once

#include "analysis.hpp"

#include <Eigen/Dense>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#if defined(__APPLE__) || defined(__unix__)
#include <csignal>
#include <sys/types.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#endif

namespace lamopt {

struct AbaqusParameterMapping {
    std::string token;
    std::size_t designIndex = 0;
};

struct AbaqusJobConfig {
    std::filesystem::path templateInputPath;
    std::string renderedInputFilename = "job.inp";
    std::filesystem::path resultFilename = "analysis_results.txt";
    std::string launchCommandTemplate;
    std::vector<AbaqusParameterMapping> parameterMappings;
    std::filesystem::path scratchRoot;
    std::chrono::milliseconds timeout{0};
    int maxRetries = 0;
    std::string backendName = "AbaqusJobBackend";
};

class AbaqusJobBackend : public AnalysisBackend {
public:
    explicit AbaqusJobBackend(AbaqusJobConfig config)
        : m_config(std::move(config)) {}

    AnalysisResult evaluate(const AnalysisRequest& request) override {
        AnalysisResult result;
        result.diagnostics.backendName = m_config.backendName;

        if (request.designVariables.size() == 0) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = "Design vector is empty.";
            return result;
        }

        std::filesystem::path runDirectory = request.workDirectory;
        if (runDirectory.empty()) {
            runDirectory = m_config.scratchRoot / "abaqus_job";
        }
        std::filesystem::create_directories(runDirectory);

        const std::filesystem::path renderedInputPath = runDirectory / m_config.renderedInputFilename;
        const std::filesystem::path resultPath = runDirectory / m_config.resultFilename;

        try {
            renderTemplate(request.designVariables, renderedInputPath);
        } catch (const std::exception& exception) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.runDirectory = runDirectory;
            result.diagnostics.message = exception.what();
            return result;
        }

        int attempts = 0;
        for (; attempts <= m_config.maxRetries; ++attempts) {
            const std::string command = fillCommandTemplate(renderedInputPath, runDirectory, resultPath);
            result.diagnostics.command = command;
            result.diagnostics.runDirectory = runDirectory;

            const CommandExecution commandExecution = runCommand(command, m_config.timeout);
            result.diagnostics.exitCode = commandExecution.exitCode;
            result.diagnostics.attempts = attempts + 1;

            if (commandExecution.status == AnalysisStatus::Success) {
                return parseResultFile(resultPath, result.diagnostics);
            }

            if (commandExecution.status == AnalysisStatus::Timeout) {
                result.status = AnalysisStatus::Timeout;
                result.diagnostics.message = commandExecution.message;
                return result;
            }

            result.status = AnalysisStatus::BackendFailure;
            result.diagnostics.message = commandExecution.message;
        }

        return result;
    }

private:
    struct CommandExecution {
        AnalysisStatus status = AnalysisStatus::BackendFailure;
        int exitCode = 0;
        std::string message;
    };

    AbaqusJobConfig m_config;

    void renderTemplate(const Eigen::VectorXd& designVariables,
                        const std::filesystem::path& renderedInputPath) const {
        std::ifstream templateInput(m_config.templateInputPath);
        if (!templateInput) {
            throw std::runtime_error("Unable to open Abaqus template input: "
                                     + m_config.templateInputPath.string());
        }

        std::ostringstream buffer;
        buffer << templateInput.rdbuf();
        std::string rendered = buffer.str();

        for (const AbaqusParameterMapping& mapping : m_config.parameterMappings) {
            if (mapping.designIndex >= static_cast<std::size_t>(designVariables.size())) {
                throw std::runtime_error("Design index out of bounds in Abaqus parameter mapping.");
            }
            replaceAll(rendered, mapping.token, formatDouble(designVariables(static_cast<Eigen::Index>(mapping.designIndex))));
        }

        std::ofstream output(renderedInputPath);
        output << rendered;
    }

    [[nodiscard]] std::string fillCommandTemplate(const std::filesystem::path& inputPath,
                                                  const std::filesystem::path& runDirectory,
                                                  const std::filesystem::path& resultPath) const {
        std::string command = m_config.launchCommandTemplate;
        replaceAll(command, "{input_file}", inputPath.string());
        replaceAll(command, "{run_dir}", runDirectory.string());
        replaceAll(command, "{result_file}", resultPath.string());
        replaceAll(command, "{job_name}", inputPath.stem().string());
        return command;
    }

    [[nodiscard]] static std::string formatDouble(double value) {
        std::ostringstream stream;
        stream.setf(std::ios::scientific);
        stream.precision(16);
        stream << value;
        return stream.str();
    }

    static void replaceAll(std::string& input,
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

    [[nodiscard]] static AnalysisResult parseResultFile(const std::filesystem::path& resultPath,
                                                        const AnalysisDiagnostics& diagnostics) {
        AnalysisResult result;
        result.diagnostics = diagnostics;

        std::ifstream resultStream(resultPath);
        if (!resultStream) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = "Abaqus result file not found: " + resultPath.string();
            return result;
        }

        std::unordered_map<std::string, std::string> fields;
        std::string line;
        while (std::getline(resultStream, line)) {
            const std::size_t separator = line.find('=');
            if (separator == std::string::npos) {
                continue;
            }
            fields[line.substr(0, separator)] = line.substr(separator + 1);
        }

        if (fields.count("objectives") == 0 || fields.count("constraints") == 0) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = "Abaqus result file is missing objectives or constraints.";
            return result;
        }

        result.objectives = parseVector(fields["objectives"]);
        result.constraints = parseVector(fields["constraints"]);
        if (fields.count("objective_gradients") != 0) {
            result.objectiveGradients = parseMatrix(fields["objective_gradients"]);
        }
        if (fields.count("constraint_gradients") != 0) {
            result.constraintGradients = parseMatrix(fields["constraint_gradients"]);
        }
        if (fields.count("objective_curvature") != 0) {
            result.objectiveCurvature = parseMatrix(fields["objective_curvature"]);
        }

        result.status = AnalysisStatus::Success;
        if (fields.count("message") != 0) {
            result.diagnostics.message = fields["message"];
        }
        return result;
    }

    [[nodiscard]] static Eigen::VectorXd parseVector(const std::string& value) {
        std::vector<double> entries;
        std::stringstream stream(value);
        std::string item;
        while (std::getline(stream, item, ',')) {
            if (!item.empty()) {
                entries.push_back(std::stod(item));
            }
        }

        Eigen::VectorXd vector(entries.size());
        for (std::size_t index = 0; index < entries.size(); ++index) {
            vector(static_cast<Eigen::Index>(index)) = entries[index];
        }
        return vector;
    }

    [[nodiscard]] static Eigen::MatrixXd parseMatrix(const std::string& value) {
        std::vector<Eigen::VectorXd> rows;
        std::stringstream stream(value);
        std::string row;
        while (std::getline(stream, row, ';')) {
            if (!row.empty()) {
                rows.push_back(parseVector(row));
            }
        }

        if (rows.empty()) {
            return Eigen::MatrixXd();
        }

        Eigen::MatrixXd matrix(rows[0].size(), static_cast<Eigen::Index>(rows.size()));
        for (std::size_t column = 0; column < rows.size(); ++column) {
            matrix.col(static_cast<Eigen::Index>(column)) = rows[column];
        }
        return matrix;
    }

    [[nodiscard]] static CommandExecution runCommand(const std::string& command,
                                                     std::chrono::milliseconds timeout) {
#if defined(__APPLE__) || defined(__unix__)
        CommandExecution result;
        pid_t pid = fork();
        if (pid < 0) {
            result.status = AnalysisStatus::BackendFailure;
            result.message = "Failed to fork Abaqus job process.";
            return result;
        }

        if (pid == 0) {
            setpgid(0, 0);
            execl("/bin/sh", "sh", "-lc", command.c_str(), static_cast<char*>(nullptr));
            _exit(127);
        }

        setpgid(pid, pid);
        int status = 0;
        const auto deadline = timeout.count() > 0
            ? std::optional<std::chrono::steady_clock::time_point>(std::chrono::steady_clock::now() + timeout)
            : std::nullopt;

        while (true) {
            const pid_t waitResult = waitpid(pid, &status, WNOHANG);
            if (waitResult == pid) {
                break;
            }
            if (waitResult < 0) {
                result.status = AnalysisStatus::BackendFailure;
                result.message = "Failed while waiting for Abaqus job process.";
                return result;
            }
            if (deadline.has_value() && std::chrono::steady_clock::now() >= *deadline) {
                killpg(pid, SIGKILL);
                waitpid(pid, &status, 0);
                result.status = AnalysisStatus::Timeout;
                result.message = "Abaqus job exceeded the configured timeout.";
                return result;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
            result.status = AnalysisStatus::Success;
            result.exitCode = 0;
            result.message = "Abaqus job completed successfully.";
            return result;
        }

        result.status = AnalysisStatus::BackendFailure;
        result.exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : status;
        result.message = "Abaqus job failed with exit code " + std::to_string(result.exitCode) + ".";
        return result;
#else
        (void)timeout;
        const int exitCode = std::system(command.c_str());
        CommandExecution result;
        result.exitCode = exitCode;
        result.status = exitCode == 0 ? AnalysisStatus::Success : AnalysisStatus::BackendFailure;
        result.message = exitCode == 0
            ? "Abaqus job completed successfully."
            : "Abaqus job failed with exit code " + std::to_string(exitCode) + ".";
        return result;
#endif
    }
};

}  // namespace lamopt
