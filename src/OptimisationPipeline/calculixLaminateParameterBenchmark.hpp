#pragma once

#include "calculixJobBackend.hpp"
#include "laminateFeParameterisation.hpp"
#include "laminateParameterisedBackend.hpp"
#include "laminationParameterDerivatives.hpp"
#include "responseSchema.hpp"

#include <fstream>
#include <sstream>

namespace lamopt {

struct CalculixLaminateParameterBenchmarkSpec {
    double plateLength = 240.0;
    double plateWidth = 160.0;
    double baseE1 = 9000.0;
    double baseE2 = 7000.0;
    double baseG12 = 3200.0;
    double nu12 = 0.24;
    double load = -1.0;
    double tipDisplacementLimit = 0.02;
    double membraneScale = 0.35;
    double bendingScale = 0.20;
    double shearScale = 0.20;
};

inline std::vector<DerivedScalarResponseRule>
MakeCalculixLaminateParameterBenchmarkObjectiveRules(const CalculixLaminateParameterBenchmarkSpec& spec = {}) {
    (void) spec;
    return {{
        "thickness",
        DerivedScalarResponseTransform::Identity,
        1.0,
        0.0,
        "thickness_objective"
    }};
}

inline std::vector<DerivedScalarResponseRule>
MakeCalculixLaminateParameterBenchmarkConstraintRules(const CalculixLaminateParameterBenchmarkSpec& spec = {}) {
    return {{
        "tip_u3",
        DerivedScalarResponseTransform::AbsoluteAffine,
        1.0 / spec.tipDisplacementLimit,
        -1.0,
        "tip_displacement_margin"
    }};
}

inline ResponseSensitivityPolicy
MakeCalculixLaminateParameterBenchmarkSensitivityPolicy(const CalculixLaminateParameterBenchmarkSpec& spec = {}) {
    (void) spec;
    return ResponseSensitivityPolicy{
        {ResponseSensitivitySource::OptimiserSideLaminate},
        {ResponseSensitivitySource::FiniteDifference}
    };
}

inline SeparableLaminateQuantityModel
MakeCalculixLaminateParameterBenchmarkQuantityModel(const AnalysisRequest& request) {
    if (request.laminateSections.size() != 1) {
        throw std::runtime_error("Laminate-parameter benchmark quantity model expects exactly one laminate section.");
    }

    const LaminateSectionState& sectionState = request.laminateSections.front();
    const Eigen::Index sectionSize = CanonicalLaminateSectionVariableCount(sectionState);

    SeparableLaminateSectionResponseTerm thicknessTerm;
    thicknessTerm.linearCoefficients = Eigen::VectorXd::Zero(sectionSize);
    thicknessTerm.quadraticCoefficients = Eigen::VectorXd::Zero(sectionSize);
    thicknessTerm.linearCoefficients(sectionSize - 1) = 1.0;

    SeparableLaminateQuantityModel quantityModel;
    quantityModel.quantities.push_back({
        "thickness",
        SeparableLaminateResponse{0.0, DenseDesignResponseTerm{}, {thicknessTerm}}
    });
    return quantityModel;
}

inline AssembledLaminationParameterDerivativeProvider
MakeCalculixLaminateParameterBenchmarkDerivativeProvider(
    const AnalysisRequest& request,
    const CalculixLaminateParameterBenchmarkSpec& spec = {}) {
    SeparableLaminationParameterDerivativeOptions options;
    options.allowPartialRuleCoverage = true;
    return AssembledLaminationParameterDerivativeProvider(
        MakeCalculixLaminateParameterBenchmarkQuantityModel(request),
        MakeCalculixLaminateParameterBenchmarkObjectiveRules(spec),
        MakeCalculixLaminateParameterBenchmarkConstraintRules(spec),
        options
    );
}

inline std::string FormatLaminateParameterBenchmarkDouble(const double value) {
    std::ostringstream stream;
    stream.setf(std::ios::scientific);
    stream.precision(16);
    stream << value;
    return stream.str();
}

inline void WriteCalculixLaminateParameterBenchmarkTemplate(
    const std::filesystem::path& templatePath,
    const CalculixLaminateParameterBenchmarkSpec& spec = {}) {
    std::ofstream output(templatePath);
    if (!output) {
        throw std::runtime_error("Unable to write CalculiX laminate-parameter benchmark template: "
                                 + templatePath.string());
    }

    const double halfLength = 0.5 * spec.plateLength;
    const double halfWidth = 0.5 * spec.plateWidth;

    output << "*HEADING\n";
    output << "LaminateOptimiser canonical laminate parameter benchmark\n";
    output << "** SEC0_A_LP0={{SEC0_A_LP0}}\n";
    output << "** SEC0_A_LP1={{SEC0_A_LP1}}\n";
    output << "** SEC0_D_LP0={{SEC0_D_LP0}}\n";
    output << "** SEC0_D_LP1={{SEC0_D_LP1}}\n";
    output << "** SEC0_THICKNESS={{SEC0_THICKNESS}}\n";
    output << "*NODE,NSET=NALL\n";
    output << "1,0.,0.,0.\n";
    output << "2," << FormatLaminateParameterBenchmarkDouble(spec.plateLength) << ",0.,0.\n";
    output << "3," << FormatLaminateParameterBenchmarkDouble(spec.plateLength) << ","
           << FormatLaminateParameterBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "4,0.," << FormatLaminateParameterBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "5," << FormatLaminateParameterBenchmarkDouble(halfLength) << ",0.,0.\n";
    output << "6," << FormatLaminateParameterBenchmarkDouble(spec.plateLength) << ","
           << FormatLaminateParameterBenchmarkDouble(halfWidth) << ",0.\n";
    output << "7," << FormatLaminateParameterBenchmarkDouble(halfLength) << ","
           << FormatLaminateParameterBenchmarkDouble(spec.plateWidth) << ",0.\n";
    output << "8,0.," << FormatLaminateParameterBenchmarkDouble(halfWidth) << ",0.\n";

    output << "*ELEMENT,TYPE=S8R,ELSET=PLATE\n";
    output << "1,1,2,3,4,5,6,7,8\n";

    output << "*NSET,NSET=X0\n";
    output << "1,8,4\n";
    output << "*NSET,NSET=Y0\n";
    output << "1,5,2\n";
    output << "*NSET,NSET=TIP\n";
    output << "5\n";

    output << "*MATERIAL,NAME=ORTHO\n";
    output << "*ELASTIC,TYPE=ENGINEERING CONSTANTS\n";
    output << "{{LAM_E1}},{{LAM_E2}},{{LAM_E2}},{{LAM_NU12}},{{LAM_NU12}},{{LAM_NU12}},{{LAM_G12}},{{LAM_G12}}\n";
    output << "{{LAM_G12}}\n";
    output << "*SHELL SECTION,ELSET=PLATE,MATERIAL=ORTHO\n";
    output << "{{SEC0_THICKNESS}}\n";

    output << "*BOUNDARY\n";
    output << "X0,1,1\n";
    output << "X0,3,3\n";
    output << "Y0,2,2\n";
    output << "Y0,3,3\n";
    output << "1,2,2\n";

    output << "*STEP\n";
    output << "*STATIC\n";
    output << "*CLOAD\n";
    output << "5,3," << FormatLaminateParameterBenchmarkDouble(spec.load) << "\n";
    output << "*NODE PRINT,NSET=TIP,GLOBAL=YES\n";
    output << "U\n";
    output << "*END STEP\n";
}

inline CalculixJobSetup MakeCalculixLaminateParameterBenchmarkJobSetup(
    const std::filesystem::path& templateInputPath,
    const std::filesystem::path& scratchRoot,
    const std::optional<std::filesystem::path>& executablePath = std::nullopt) {
    CalculixJobSetup setup;
    setup.templateInputPath = templateInputPath;
    setup.scratchRoot = scratchRoot;
    setup.executablePath = executablePath;
    setup.rawScalarExtractions = {
        {"tip_u3",
         {std::filesystem::path("{job_name}.dat"),
          "displacements \\(vx,vy,vz\\) for set TIP and time\\s+[-+0-9.Ee]+\\s*\\n\\s*5\\s+[-+0-9.Ee]+\\s+[-+0-9.Ee]+\\s+([-+0-9.Ee]+)",
          1,
          CalculixMatchSelection::Last,
          1.0,
          0.0,
          "tip_u3"}}
    };
    return setup;
}

class CalculixLaminateParameterBenchmarkBackend final : public AnalysisBackend {
public:
    explicit CalculixLaminateParameterBenchmarkBackend(CalculixJobSetup setup,
                                                       CalculixLaminateParameterBenchmarkSpec spec = {})
        : m_jobBackend(MakeDefaultCalculixJobConfig(setup))
        , m_parameterisedBackend(m_jobBackend)
        , m_spec(std::move(spec)) {}

    AnalysisResult evaluate(const AnalysisRequest& request) override {
        AnalysisRequest parameterisedRequest = request;

        CanonicalLaminateSectionLayout layout;
        std::string message;
        if (!BuildCanonicalLaminateSectionLayout(request, 0, layout, message)) {
            AnalysisResult result;
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message = message;
            return result;
        }
        if (request.laminateSections.size() != 1) {
            AnalysisResult result;
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message =
                "CalculiX laminate-parameter benchmark expects exactly one laminate section.";
            return result;
        }
        if (layout.laminationParameterCount < 2 || layout.laminationParameterSlots.size() < 4) {
            AnalysisResult result;
            result.status = AnalysisStatus::InvalidOutput;
            result.diagnostics.message =
                "CalculiX laminate-parameter benchmark expects a balanced symmetric laminate section layout.";
            return result;
        }

        const double a0 = request.designVariables(layout.laminationParameterSlots[0].globalIndex);
        const double a1 = request.designVariables(layout.laminationParameterSlots[1].globalIndex);
        const double d0 = request.designVariables(layout.laminationParameterSlots[2].globalIndex);
        const double d1 = request.designVariables(layout.laminationParameterSlots[3].globalIndex);

        parameterisedRequest.templateParameters.push_back(
            {"{{LAM_E1}}", FormatLaminateParameterBenchmarkDouble(
                               materialFactor(m_spec.baseE1, a0, d0, m_spec.membraneScale, m_spec.bendingScale))});
        parameterisedRequest.templateParameters.push_back(
            {"{{LAM_E2}}", FormatLaminateParameterBenchmarkDouble(
                               materialFactor(m_spec.baseE2, a1, d1, m_spec.membraneScale, m_spec.bendingScale))});
        parameterisedRequest.templateParameters.push_back(
            {"{{LAM_G12}}", FormatLaminateParameterBenchmarkDouble(
                                shearFactor(m_spec.baseG12, d0, d1, m_spec.shearScale))});
        parameterisedRequest.templateParameters.push_back(
            {"{{LAM_NU12}}", FormatLaminateParameterBenchmarkDouble(m_spec.nu12)});

        AnalysisResult result = m_parameterisedBackend.evaluate(parameterisedRequest);
        if (!result.isSuccessful()) {
            return result;
        }

        const double thickness = request.designVariables(layout.thicknessIndex);
        result.extractedScalarValues["thickness"] = thickness;
        const auto objectiveRules = MakeCalculixLaminateParameterBenchmarkObjectiveRules(m_spec);
        const auto constraintRules = MakeCalculixLaminateParameterBenchmarkConstraintRules(m_spec);
        result.objectives = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(objectiveRules.size()));
        result.constraints = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(constraintRules.size()));
        for (std::size_t index = 0; index < objectiveRules.size(); ++index) {
            result.objectives(static_cast<Eigen::Index>(index)) =
                EvaluateDerivedScalarResponseRule(objectiveRules[index], result.extractedScalarValues);
        }
        for (std::size_t index = 0; index < constraintRules.size(); ++index) {
            result.constraints(static_cast<Eigen::Index>(index)) =
                EvaluateDerivedScalarResponseRule(constraintRules[index], result.extractedScalarValues);
        }
        result.sensitivityPolicy = MakeCalculixLaminateParameterBenchmarkSensitivityPolicy(m_spec);

        if (!result.diagnostics.message.empty() && result.diagnostics.message.back() != ' ') {
            result.diagnostics.message += ' ';
        }
        result.diagnostics.message += "Canonical laminate template parameters and response-schema values assembled for real CalculiX output.";
        return result;
    }

private:
    [[nodiscard]] static double materialFactor(const double baseValue,
                                               const double membraneValue,
                                               const double bendingValue,
                                               const double membraneScale,
                                               const double bendingScale) {
        const double factor = 1.15 + membraneScale * membraneValue + bendingScale * bendingValue;
        return baseValue * std::max(factor, 0.20);
    }

    [[nodiscard]] static double shearFactor(const double baseValue,
                                            const double bending0,
                                            const double bending1,
                                            const double scale) {
        const double factor = 1.10 + scale * 0.5 * (bending0 + bending1);
        return baseValue * std::max(factor, 0.20);
    }

    CalculixJobBackend m_jobBackend;
    LaminateParameterisedBackend m_parameterisedBackend;
    CalculixLaminateParameterBenchmarkSpec m_spec;
};

}  // namespace lamopt
