#pragma once

#include "jobBackend.hpp"

namespace lamopt {

using CalculixParameterMapping = JobParameterMapping;

struct CalculixJobConfig : JobBackendConfig {
    CalculixJobConfig() {
        renderedInputFilename = "job.inp";
        resultFilename = "analysis_results.txt";
        launchCommandTemplate = "ccx -i {job_name}";
        backendName = "CalculixJobBackend";
        defaultRunDirectoryName = "calculix_job";
    }
};

class CalculixJobBackend : public JobBackend {
public:
    explicit CalculixJobBackend(CalculixJobConfig config)
        : JobBackend(std::move(config)) {}
};

}  // namespace lamopt
