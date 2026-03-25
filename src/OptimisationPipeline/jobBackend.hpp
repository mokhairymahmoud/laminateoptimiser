#pragma once

#include "analysis.hpp"

#include <Eigen/Dense>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#if defined(__APPLE__) || defined(__unix__)
#include <csignal>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#endif

namespace lamopt {

struct JobParameterMapping {
    std::string token;
    std::size_t designIndex = 0;
};

struct JobBackendConfig {
    std::filesystem::path templateInputPath;
    std::string renderedInputFilename = "job.inp";
    std::filesystem::path resultFilename = "analysis_results.txt";
    std::string launchCommandTemplate;
    std::vector<JobParameterMapping> parameterMappings;
    std::filesystem::path scratchRoot;
    std::chrono::milliseconds timeout{0};
    int maxRetries = 0;
    bool launchInRunDirectory = false;
    std::optional<std::filesystem::path> standardOutputFilename;
    std::optional<std::filesystem::path> standardErrorFilename;
    std::string backendName = "ExternalJobBackend";
    std::string defaultRunDirectoryName = "external_job";
};

class JobBackend : public AnalysisBackend {
public:
    explicit JobBackend(JobBackendConfig config)
        : m_config(std::move(config)) {}

    AnalysisResult evaluate(const AnalysisRequest& request) override {
        AnalysisResult result;
        result.diagnostics.backendName = backendName();

        if (request.designVariables.size() == 0) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = "Design vector is empty.";
            return result;
        }

        std::filesystem::path runDirectory = request.workDirectory;
        if (runDirectory.empty()) {
            runDirectory = m_config.scratchRoot / defaultRunDirectoryName();
        }
        std::filesystem::create_directories(runDirectory);

        const ExecutionPaths paths{
            runDirectory,
            runDirectory / m_config.renderedInputFilename,
            runDirectory / m_config.resultFilename,
            m_config.standardOutputFilename.has_value()
            ? std::optional<std::filesystem::path>(runDirectory / *m_config.standardOutputFilename)
            : std::nullopt,
            m_config.standardErrorFilename.has_value()
            ? std::optional<std::filesystem::path>(runDirectory / *m_config.standardErrorFilename)
            : std::nullopt
        };
        result.diagnostics.standardOutputPath = paths.standardOutputPath.value_or(std::filesystem::path());
        result.diagnostics.standardErrorPath = paths.standardErrorPath.value_or(std::filesystem::path());

        try {
            renderTemplate(request, paths.renderedInputPath);
        } catch (const std::exception& exception) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.runDirectory = runDirectory;
            result.diagnostics.message = exception.what();
            return result;
        }

        int attempts = 0;
        for (; attempts <= m_config.maxRetries; ++attempts) {
            const std::string command = fillCommandTemplate(paths.renderedInputPath, runDirectory, paths.resultPath);
            result.diagnostics.command = command;
            result.diagnostics.runDirectory = runDirectory;

            const CommandExecution commandExecution = runCommand(
                command,
                m_config.timeout,
                backendName(),
                m_config.launchInRunDirectory ? runDirectory : std::filesystem::path(),
                paths.standardOutputPath,
                paths.standardErrorPath
            );
            result.diagnostics.exitCode = commandExecution.exitCode;
            result.diagnostics.attempts = attempts + 1;

            if (commandExecution.status == AnalysisStatus::Success) {
                return parseSuccessfulRun(paths, result.diagnostics);
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

protected:
    struct ExecutionPaths {
        std::filesystem::path runDirectory;
        std::filesystem::path renderedInputPath;
        std::filesystem::path resultPath;
        std::optional<std::filesystem::path> standardOutputPath;
        std::optional<std::filesystem::path> standardErrorPath;
    };

    [[nodiscard]] const JobBackendConfig& config() const {
        return m_config;
    }

    [[nodiscard]] virtual AnalysisResult parseSuccessfulRun(const ExecutionPaths& paths,
                                                            const AnalysisDiagnostics& diagnostics) const {
        return parseResultFile(paths.resultPath, diagnostics, backendName());
    }

    [[nodiscard]] static AnalysisResult parseResultFile(const std::filesystem::path& resultPath,
                                                        const AnalysisDiagnostics& diagnostics,
                                                        const std::string& backendName) {
        AnalysisResult result;
        result.diagnostics = diagnostics;

        std::ifstream resultStream(resultPath);
        if (!resultStream) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = backendName + " result file not found: " + resultPath.string();
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
            result.diagnostics.message = backendName
                + " result file is missing objectives or constraints.";
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
        if (fields.count("constraint_curvature") != 0) {
            result.constraintCurvature = parseMatrix(fields["constraint_curvature"]);
        }

        result.status = AnalysisStatus::Success;
        if (fields.count("message") != 0) {
            result.diagnostics.message = fields["message"];
        }
        return result;
    }

    static void writeResultFile(const std::filesystem::path& resultPath,
                                const Eigen::VectorXd& objectives,
                                const Eigen::VectorXd& constraints,
                                const std::optional<Eigen::MatrixXd>& objectiveGradients = std::nullopt,
                                const std::optional<Eigen::MatrixXd>& constraintGradients = std::nullopt,
                                const std::optional<Eigen::MatrixXd>& objectiveCurvature = std::nullopt,
                                const std::optional<Eigen::MatrixXd>& constraintCurvature = std::nullopt,
                                const std::optional<std::string>& message = std::nullopt) {
        std::ofstream resultStream(resultPath);
        if (!resultStream) {
            throw std::runtime_error("Unable to write job result file: " + resultPath.string());
        }

        resultStream << "objectives=" << formatVector(objectives) << '\n';
        resultStream << "constraints=" << formatVector(constraints) << '\n';
        if (objectiveGradients.has_value()) {
            resultStream << "objective_gradients=" << formatMatrix(*objectiveGradients) << '\n';
        }
        if (constraintGradients.has_value()) {
            resultStream << "constraint_gradients=" << formatMatrix(*constraintGradients) << '\n';
        }
        if (objectiveCurvature.has_value()) {
            resultStream << "objective_curvature=" << formatMatrix(*objectiveCurvature) << '\n';
        }
        if (constraintCurvature.has_value()) {
            resultStream << "constraint_curvature=" << formatMatrix(*constraintCurvature) << '\n';
        }
        if (message.has_value()) {
            resultStream << "message=" << *message << '\n';
        }
    }

private:
    struct CommandExecution {
        AnalysisStatus status = AnalysisStatus::BackendFailure;
        int exitCode = 0;
        std::string message;
    };

    JobBackendConfig m_config;

    [[nodiscard]] const std::string& backendName() const {
        return m_config.backendName;
    }

    [[nodiscard]] const std::string& defaultRunDirectoryName() const {
        return m_config.defaultRunDirectoryName;
    }

    void renderTemplate(const AnalysisRequest& request,
                        const std::filesystem::path& renderedInputPath) const {
        std::ifstream templateInput(m_config.templateInputPath);
        if (!templateInput) {
            throw std::runtime_error("Unable to open " + backendName()
                                     + " template input: "
                                     + m_config.templateInputPath.string());
        }

        std::ostringstream buffer;
        buffer << templateInput.rdbuf();
        std::string rendered = buffer.str();

        for (const JobParameterMapping& mapping : m_config.parameterMappings) {
            if (mapping.designIndex >= static_cast<std::size_t>(request.designVariables.size())) {
                throw std::runtime_error("Design index out of bounds in "
                                         + backendName()
                                         + " parameter mapping.");
            }
            replaceAll(rendered,
                       mapping.token,
                       formatDouble(request.designVariables(static_cast<Eigen::Index>(mapping.designIndex))));
        }

        for (const ResolvedTemplateParameter& parameter : request.templateParameters) {
            if (parameter.token.empty()) {
                throw std::runtime_error("Resolved template parameter token must not be empty in "
                                         + backendName() + ".");
            }
            replaceAll(rendered, parameter.token, parameter.value);
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

    [[nodiscard]] static std::string shellQuote(const std::string& value) {
        std::string quoted = "'";
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

    [[nodiscard]] static std::string formatVector(const Eigen::VectorXd& value) {
        std::ostringstream stream;
        for (Eigen::Index index = 0; index < value.size(); ++index) {
            if (index != 0) {
                stream << ",";
            }
            stream << formatDouble(value(index));
        }
        return stream.str();
    }

    [[nodiscard]] static std::string formatMatrix(const Eigen::MatrixXd& value) {
        std::ostringstream stream;
        for (Eigen::Index column = 0; column < value.cols(); ++column) {
            if (column != 0) {
                stream << ";";
            }
            stream << formatVector(value.col(column));
        }
        return stream.str();
    }

    [[nodiscard]] static CommandExecution runCommand(const std::string& command,
                                                     std::chrono::milliseconds timeout,
                                                     const std::string& backendName,
                                                     const std::filesystem::path& workingDirectory,
                                                     const std::optional<std::filesystem::path>& standardOutputPath,
                                                     const std::optional<std::filesystem::path>& standardErrorPath) {
#if defined(__APPLE__) || defined(__unix__)
        CommandExecution result;
        pid_t pid = fork();
        if (pid < 0) {
            result.status = AnalysisStatus::BackendFailure;
            result.message = "Failed to fork " + backendName + " job process.";
            return result;
        }

        if (pid == 0) {
            setpgid(0, 0);
            if (!workingDirectory.empty()) {
                chdir(workingDirectory.c_str());
            }
            if (standardOutputPath.has_value()) {
                const int stdoutFd = open(standardOutputPath->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
                if (stdoutFd < 0) {
                    _exit(126);
                }
                dup2(stdoutFd, STDOUT_FILENO);
                close(stdoutFd);
            }
            if (standardErrorPath.has_value()) {
                const int stderrFd = open(standardErrorPath->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
                if (stderrFd < 0) {
                    _exit(126);
                }
                dup2(stderrFd, STDERR_FILENO);
                close(stderrFd);
            }
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
                result.message = "Failed while waiting for " + backendName + " job process.";
                return result;
            }
            if (deadline.has_value() && std::chrono::steady_clock::now() >= *deadline) {
                killpg(pid, SIGKILL);
                waitpid(pid, &status, 0);
                result.status = AnalysisStatus::Timeout;
                result.message = backendName + " job exceeded the configured timeout.";
                return result;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
            result.status = AnalysisStatus::Success;
            result.exitCode = 0;
            result.message = backendName + " job completed successfully.";
            return result;
        }

        result.status = AnalysisStatus::BackendFailure;
        result.exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : status;
        result.message = backendName + " job failed with exit code " + std::to_string(result.exitCode) + ".";
        return result;
#else
        (void)timeout;
        std::string shellCommand;
        if (!workingDirectory.empty()) {
            shellCommand += "cd " + shellQuote(workingDirectory.string()) + " && ";
        }
        shellCommand += command;
        if (standardOutputPath.has_value()) {
            shellCommand += " > " + shellQuote(standardOutputPath->string());
        }
        if (standardErrorPath.has_value()) {
            shellCommand += " 2> " + shellQuote(standardErrorPath->string());
        }

        const int exitCode = std::system(shellCommand.c_str());
        CommandExecution result;
        result.exitCode = exitCode;
        result.status = exitCode == 0 ? AnalysisStatus::Success : AnalysisStatus::BackendFailure;
        result.message = exitCode == 0
            ? backendName + " job completed successfully."
            : backendName + " job failed with exit code " + std::to_string(exitCode) + ".";
        return result;
#endif
    }
};

}  // namespace lamopt
