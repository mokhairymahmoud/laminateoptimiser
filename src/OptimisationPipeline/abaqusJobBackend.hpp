#pragma once

#include "jobBackend.hpp"

namespace lamopt {

using AbaqusParameterMapping = JobParameterMapping;

struct AbaqusJobConfig : JobBackendConfig {
    AbaqusJobConfig() {
        renderedInputFilename = "job.inp";
        resultFilename = "analysis_results.txt";
        backendName = "AbaqusJobBackend";
        defaultRunDirectoryName = "abaqus_job";
    }
};

class AbaqusJobBackend : public JobBackend {
public:
    explicit AbaqusJobBackend(AbaqusJobConfig config)
        : JobBackend(std::move(config)) {}
};

}  // namespace lamopt
