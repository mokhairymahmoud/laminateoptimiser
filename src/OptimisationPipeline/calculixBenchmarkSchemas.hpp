#pragma once

#include "calculixJobBackend.hpp"

namespace lamopt {

struct CalculixPlateBucklingBenchmarkSchemaOptions {
    std::filesystem::path sourceFilename = std::filesystem::path("{job_name}.dat");
    double requiredBucklingFactor = 1.0;
    double tipDisplacementLimit = 1.0;
    double displacementScale = 1.0;
    std::string massPattern = "TOTAL MASS\\s*=\\s*([-+0-9.Ee]+)";
    std::string bucklingPattern = "BUCKLING FACTOR(?:\\s+1)?\\s*=\\s*([-+0-9.Ee]+)";
    std::string tipDisplacementPattern = "TIP DISPLACEMENT(?:\\s+U3)?\\s*=\\s*([-+0-9.Ee]+)";
    CalculixMatchSelection matchSelection = CalculixMatchSelection::Last;
};

inline void ConfigureCalculixPlateBucklingBenchmarkResponses(
    CalculixJobConfig& config,
    const CalculixPlateBucklingBenchmarkSchemaOptions& options = {}) {
    config.rawScalarExtractions = {
        {"mass",
         {options.sourceFilename,
          options.massPattern,
          1,
          options.matchSelection,
          1.0,
          0.0,
          "mass"}},
        {"buckling_lambda_1",
         {options.sourceFilename,
          options.bucklingPattern,
          1,
          options.matchSelection,
          1.0,
          0.0,
          "buckling_lambda_1"}},
        {"tip_u3",
         {options.sourceFilename,
          options.tipDisplacementPattern,
          1,
          options.matchSelection,
          options.displacementScale,
          0.0,
          "tip_u3"}}
    };

    config.objectiveResponses = {
        {"mass",
         CalculixDerivedResponseTransform::Identity,
         1.0,
         0.0,
         "objective_mass"}
    };

    config.constraintResponses = {
        {"buckling_lambda_1",
         CalculixDerivedResponseTransform::InverseAffine,
         options.requiredBucklingFactor,
         -1.0,
         "buckling_margin"},
        {"tip_u3",
         CalculixDerivedResponseTransform::AbsoluteAffine,
         1.0 / options.tipDisplacementLimit,
         -1.0,
         "tip_displacement_margin"}
    };
}

inline void ConfigureCalculixPlateBucklingBenchmarkResponses(
    CalculixJobSetup& setup,
    const CalculixPlateBucklingBenchmarkSchemaOptions& options = {}) {
    setup.rawScalarExtractions = {
        {"mass",
         {options.sourceFilename,
          options.massPattern,
          1,
          options.matchSelection,
          1.0,
          0.0,
          "mass"}},
        {"buckling_lambda_1",
         {options.sourceFilename,
          options.bucklingPattern,
          1,
          options.matchSelection,
          1.0,
          0.0,
          "buckling_lambda_1"}},
        {"tip_u3",
         {options.sourceFilename,
          options.tipDisplacementPattern,
          1,
          options.matchSelection,
          options.displacementScale,
          0.0,
          "tip_u3"}}
    };

    setup.objectiveResponses = {
        {"mass",
         CalculixDerivedResponseTransform::Identity,
         1.0,
         0.0,
         "objective_mass"}
    };

    setup.constraintResponses = {
        {"buckling_lambda_1",
         CalculixDerivedResponseTransform::InverseAffine,
         options.requiredBucklingFactor,
         -1.0,
         "buckling_margin"},
        {"tip_u3",
         CalculixDerivedResponseTransform::AbsoluteAffine,
         1.0 / options.tipDisplacementLimit,
         -1.0,
         "tip_displacement_margin"}
    };
}

}  // namespace lamopt
