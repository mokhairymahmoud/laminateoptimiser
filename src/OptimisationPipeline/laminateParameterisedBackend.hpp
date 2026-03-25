#pragma once

#include "laminateFeParameterisation.hpp"

namespace lamopt {

class LaminateParameterisedBackend final : public AnalysisBackend {
public:
    explicit LaminateParameterisedBackend(AnalysisBackend& innerBackend,
                                          LaminateTemplateParameterisationOptions options = {})
        : m_innerBackend(innerBackend)
        , m_options(std::move(options)) {}

    AnalysisResult evaluate(const AnalysisRequest& request) override {
        AnalysisRequest parameterisedRequest = request;

        std::vector<ResolvedTemplateParameter> laminateParameters;
        std::string message;
        if (!BuildLaminateTemplateParameters(request, laminateParameters, message, m_options)) {
            AnalysisResult result;
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = message;
            return result;
        }

        parameterisedRequest.templateParameters.reserve(request.templateParameters.size() + laminateParameters.size());
        parameterisedRequest.templateParameters.insert(parameterisedRequest.templateParameters.end(),
                                                      laminateParameters.begin(),
                                                      laminateParameters.end());
        return m_innerBackend.evaluate(parameterisedRequest);
    }

private:
    AnalysisBackend& m_innerBackend;
    LaminateTemplateParameterisationOptions m_options;
};

}  // namespace lamopt
