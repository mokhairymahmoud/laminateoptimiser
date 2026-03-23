#pragma once

#include "calculixJobBackend.hpp"
#include "responseSchema.hpp"

#include <fstream>

namespace lamopt {

struct CalculixComposipyBenchmarkSpec {
    double plateLength = 360.0;
    double plateWidth = 360.0;
    std::vector<double> plyAngles = {-45.0, 45.0, 90.0, 0.0, 0.0, 0.0, 0.0, 90.0, 45.0, -45.0};
    double referencePlyThickness = 0.21;
    double e1 = 60800.0;
    double e2 = 58250.0;
    double e3 = 58250.0;
    double nu12 = 0.07;
    double nu13 = 0.07;
    double nu23 = 0.07;
    double g12 = 4550.0;
    double g13 = 4550.0;
    double g23 = 4550.0;
    double transverseLoad = -1.0;
    double edgeCompressionCornerLoad = -60.0;
    double edgeCompressionMidLoad = -120.0;
    double requiredBucklingFactor = 2.0;
    double tipDisplacementLimit = 0.75;
};

inline std::vector<DerivedScalarResponseRule>
MakeCalculixComposipyBenchmarkObjectiveRules(const CalculixComposipyBenchmarkSpec& spec = {}) {
    (void) spec;
    return {{
        "mass_per_area",
        DerivedScalarResponseTransform::Identity,
        1.0,
        0.0,
        "mass_per_area"
    }};
}

inline std::vector<DerivedScalarResponseRule>
MakeCalculixComposipyBenchmarkConstraintRules(const CalculixComposipyBenchmarkSpec& spec = {}) {
    return {{
                "buckling_lambda_1",
                DerivedScalarResponseTransform::InverseAffine,
                spec.requiredBucklingFactor,
                -1.0,
                "buckling_margin"
            },
            {
                "tip_u3",
                DerivedScalarResponseTransform::AbsoluteAffine,
                1.0 / spec.tipDisplacementLimit,
                -1.0,
                "tip_displacement_margin"
            }};
}

inline ResponseSensitivityPolicy
MakeCalculixComposipyBenchmarkSensitivityPolicy(const CalculixComposipyBenchmarkSpec& spec = {}) {
    (void) spec;
    ResponseSensitivityPolicy policy;
    policy.objectiveSources = {ResponseSensitivitySource::BackendNative};
    policy.constraintSources = {ResponseSensitivitySource::FiniteDifference,
                                ResponseSensitivitySource::FiniteDifference};
    return policy;
}

inline std::string FormatBenchmarkDouble(const double value) {
    std::ostringstream stream;
    stream.setf(std::ios::scientific);
    stream.precision(16);
    stream << value;
    return stream.str();
}

inline std::string OrientationNameForAngle(const double angle) {
    if (std::abs(angle) < 1.0e-12) {
        return "OR0";
    }
    if (std::abs(angle - 45.0) < 1.0e-12) {
        return "OR45";
    }
    if (std::abs(angle + 45.0) < 1.0e-12) {
        return "ORM45";
    }
    if (std::abs(angle - 90.0) < 1.0e-12) {
        return "OR90";
    }
    throw std::runtime_error("Unsupported benchmark ply angle: " + std::to_string(angle));
}

inline void WriteCalculixComposipyBenchmarkTemplate(
    const std::filesystem::path& templatePath,
    const CalculixComposipyBenchmarkSpec& spec = {},
    const std::string& plyThicknessToken = "{{PLY_THICKNESS}}") {
    std::ofstream output(templatePath);
    if (!output) {
        throw std::runtime_error("Unable to write CalculiX Composipy benchmark template: "
                                 + templatePath.string());
    }

    const double halfLength = 0.5 * spec.plateLength;
    const double halfWidth = 0.5 * spec.plateWidth;

    output << "*HEADING\n";
    output << "LaminateOptimiser CalculiX Composipy benchmark plate\n";
    output << "*NODE,NSET=NALL\n";
    output << "1,0.,0.,0.\n";
    output << "2," << FormatBenchmarkDouble(spec.plateLength) << ",0.,0.\n";
    output << "3," << FormatBenchmarkDouble(spec.plateLength) << ","
           << FormatBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "4,0.," << FormatBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "5," << FormatBenchmarkDouble(halfLength) << ",0.,0.\n";
    output << "6," << FormatBenchmarkDouble(spec.plateLength) << ","
           << FormatBenchmarkDouble(halfWidth) << ",0.\n";
    output << "7," << FormatBenchmarkDouble(halfLength) << ","
           << FormatBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "8,0.," << FormatBenchmarkDouble(halfWidth) << ",0.\n";

    output << "*ELEMENT,TYPE=S8R,ELSET=PLATE\n";
    output << "1,1,2,3,4,5,6,7,8\n";

    output << "*NSET,NSET=X0\n";
    output << "1,8,4\n";
    output << "*NSET,NSET=XA\n";
    output << "2,6,3\n";
    output << "*NSET,NSET=Y0\n";
    output << "1,5,2\n";
    output << "*NSET,NSET=YB\n";
    output << "4,7,3\n";
    output << "*NSET,NSET=TIP\n";
    output << "5\n";

    output << "*MATERIAL,NAME=LAMINA\n";
    output << "*ELASTIC,TYPE=ENGINEERING CONSTANTS\n";
    output << FormatBenchmarkDouble(spec.e1) << ","
           << FormatBenchmarkDouble(spec.e2) << ","
           << FormatBenchmarkDouble(spec.e3) << ","
           << FormatBenchmarkDouble(spec.nu12) << ","
           << FormatBenchmarkDouble(spec.nu13) << ","
           << FormatBenchmarkDouble(spec.nu23) << ","
           << FormatBenchmarkDouble(spec.g12) << ","
           << FormatBenchmarkDouble(spec.g13) << "\n";
    output << FormatBenchmarkDouble(spec.g23) << "\n";

    output << "*ORIENTATION,NAME=OR0\n";
    output << "1.,0.,0.,0.,1.,0.\n";
    output << "*ORIENTATION,NAME=OR45\n";
    output << "1.,0.,0.,0.,1.,0.\n";
    output << "3,45.\n";
    output << "*ORIENTATION,NAME=ORM45\n";
    output << "1.,0.,0.,0.,1.,0.\n";
    output << "3,-45.\n";
    output << "*ORIENTATION,NAME=OR90\n";
    output << "1.,0.,0.,0.,1.,0.\n";
    output << "3,90.\n";

    output << "*SHELL SECTION,ELSET=PLATE,COMPOSITE\n";
    for (const double angle : spec.plyAngles) {
        output << plyThicknessToken << ",,LAMINA," << OrientationNameForAngle(angle) << "\n";
    }

    output << "*BOUNDARY\n";
    output << "X0,1,1\n";
    output << "X0,3,3\n";
    output << "YB,2,2\n";
    output << "YB,3,3\n";
    output << "4,2,2\n";

    output << "*STEP\n";
    output << "*STATIC\n";
    output << "*CLOAD\n";
    output << "5,3," << FormatBenchmarkDouble(spec.transverseLoad) << "\n";
    output << "*NODE PRINT,NSET=TIP,GLOBAL=YES\n";
    output << "U\n";
    output << "*END STEP\n";

    output << "*STEP\n";
    output << "*BUCKLE\n";
    output << "1\n";
    output << "*CLOAD\n";
    output << "2,1," << FormatBenchmarkDouble(spec.edgeCompressionCornerLoad) << "\n";
    output << "6,1," << FormatBenchmarkDouble(spec.edgeCompressionMidLoad) << "\n";
    output << "3,1," << FormatBenchmarkDouble(spec.edgeCompressionCornerLoad) << "\n";
    output << "*NODE PRINT,NSET=TIP,GLOBAL=YES\n";
    output << "U\n";
    output << "*END STEP\n";
}

inline CalculixJobSetup MakeCalculixComposipyBenchmarkJobSetup(
    const std::filesystem::path& templateInputPath,
    const std::filesystem::path& scratchRoot,
    const std::optional<std::filesystem::path>& executablePath = std::nullopt) {
    CalculixJobSetup setup;
    setup.templateInputPath = templateInputPath;
    setup.scratchRoot = scratchRoot;
    setup.parameterMappings = {{"{{PLY_THICKNESS}}", 0}};
    setup.executablePath = executablePath;
    setup.rawScalarExtractions = {
        {"buckling_lambda_1",
         {std::filesystem::path("{job_name}.dat"),
          "B\\s+U\\s+C\\s+K\\s+L\\s+I\\s+N\\s+G\\s+F\\s+A\\s+C\\s+T\\s+O\\s+R\\s+O\\s+U\\s+T\\s+P\\s+U\\s+T[\\s\\S]*?\\n\\s*1\\s+([-+0-9.Ee]+)",
          1,
          CalculixMatchSelection::First,
          1.0,
          0.0,
          "buckling_lambda_1"}},
        {"tip_u3",
         {std::filesystem::path("{job_name}.dat"),
          "S\\s+T\\s+E\\s+P\\s+1[\\s\\S]*?displacements \\(vx,vy,vz\\) for set TIP and time\\s+[-+0-9.Ee]+\\s*\\n\\s*5\\s+[-+0-9.Ee]+\\s+[-+0-9.Ee]+\\s+([-+0-9.Ee]+)",
          1,
          CalculixMatchSelection::First,
          1.0,
          0.0,
          "tip_u3"}}
    };
    return setup;
}

class CalculixComposipyBenchmarkBackend final : public AnalysisBackend {
public:
    explicit CalculixComposipyBenchmarkBackend(CalculixJobSetup setup,
                                               CalculixComposipyBenchmarkSpec spec = {})
        : m_backend(MakeDefaultCalculixJobConfig(setup))
        , m_spec(std::move(spec)) {}

    AnalysisResult evaluate(const AnalysisRequest& request) override {
        AnalysisResult result = m_backend.evaluate(request);
        if (!result.isSuccessful()) {
            return result;
        }

        if (request.designVariables.size() != 1) {
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message =
                "CalculiX Composipy benchmark backend expects exactly one design variable (ply thickness).";
            return result;
        }

        const double plyThickness = request.designVariables(0);
        const double massPerArea = static_cast<double>(m_spec.plyAngles.size()) * plyThickness;
        result.extractedScalarValues["mass_per_area"] = massPerArea;
        result.sensitivityPolicy = MakeCalculixComposipyBenchmarkSensitivityPolicy(m_spec);

        const auto objectiveRules = MakeCalculixComposipyBenchmarkObjectiveRules(m_spec);
        const auto constraintRules = MakeCalculixComposipyBenchmarkConstraintRules(m_spec);

        result.objectives = Eigen::VectorXd::Constant(
            static_cast<Eigen::Index>(objectiveRules.size()),
            0.0
        );
        for (std::size_t index = 0; index < objectiveRules.size(); ++index) {
            result.objectives(static_cast<Eigen::Index>(index)) =
                EvaluateDerivedScalarResponseRule(objectiveRules[index], result.extractedScalarValues);
        }
        result.constraints = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(constraintRules.size()));
        for (std::size_t index = 0; index < constraintRules.size(); ++index) {
            result.constraints(static_cast<Eigen::Index>(index)) =
                EvaluateDerivedScalarResponseRule(constraintRules[index], result.extractedScalarValues);
        }

        if (request.requestSensitivities) {
            result.objectiveGradients = Eigen::MatrixXd::Constant(1, 1, static_cast<double>(m_spec.plyAngles.size()));
            result.objectiveGradientMask = Eigen::MatrixXi::Ones(1, 1);
        }

        if (!result.diagnostics.message.empty() && result.diagnostics.message.back() != ' ') {
            result.diagnostics.message += ' ';
        }
        result.diagnostics.message += "Composipy benchmark responses assembled from real CalculiX output.";
        if (request.requestSensitivities && result.sensitivityPolicy.has_value()) {
            result.diagnostics.message += " Sensitivity policy: "
                                       + DescribeResponseSensitivityPolicy(*result.sensitivityPolicy) + ".";
        }
        return result;
    }

private:
    CalculixJobBackend m_backend;
    CalculixComposipyBenchmarkSpec m_spec;
};

}  // namespace lamopt
