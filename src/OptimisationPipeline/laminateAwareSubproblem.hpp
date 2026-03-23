#pragma once

#include "laminateSectionBinding.hpp"
#include "subproblem.hpp"

#include <Eigen/Cholesky>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace lamopt {

struct LaminateProjectionOptions {
    int maxIterations = 30;
    double tolerance = 1.0e-8;
    double thicknessInteriorMargin = 1.0e-6;
};

class LaminateAwareSubproblemSolver : public SubproblemSolver {
public:
    LaminateAwareSubproblemSolver(SubproblemSolver& innerSolver,
                                  LaminateProjectionOptions options = {})
        : m_innerSolver(innerSolver)
        , m_options(options) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        SubproblemResult result = m_innerSolver.solve(problem);
        if (!result.success || problem.laminateSections.empty()) {
            return result;
        }

        std::ostringstream message;
        message << result.message;

        for (const LaminateSectionState& sectionState : problem.laminateSections) {
            ProjectionResult projection = dispatchProjection(sectionState, result.candidateDesign);
            if (!projection.success) {
                result.success = false;
                result.message = projection.message;
                return result;
            }
            result.iterations += projection.iterations;
            message << " | projected laminate section at offset "
                    << sectionState.variableOffset
                    << " in " << projection.iterations << " iterations";
        }

        message << " | laminate-aware projection complete";
        result.message = message.str();
        return result;
    }

private:
    struct ProjectionResult {
        bool success = false;
        int iterations = 0;
        double finalGap = std::numeric_limits<double>::infinity();
        std::string message;
    };

    struct ProjectionDispatchVisitor {
        using Result = ProjectionResult;

        LaminateAwareSubproblemSolver* solver = nullptr;
        Eigen::VectorXd* candidate = nullptr;

        template<typename Section>
        Result operator()(const LaminateSectionState& sectionState) const {
            return solver->template projectSection<Section>(sectionState, *candidate);
        }
    };

    SubproblemSolver& m_innerSolver;
    LaminateProjectionOptions m_options;

    ProjectionResult dispatchProjection(const LaminateSectionState& state, Eigen::VectorXd& candidate) {
        return DispatchLaminateSectionType(state, ProjectionDispatchVisitor{this, &candidate});
    }

    template<typename Section>
    ProjectionResult projectSection(const LaminateSectionState& state, Eigen::VectorXd& candidate) {
        using Vector = typename Section::Vector_t;
        using Hessian = typename Section::Hessian_t;

        ProjectionResult result;

        constexpr int kSectionSize = Section::Size;
        if (candidate.size() < state.variableOffset + kSectionSize) {
            result.message = "Candidate design is too small for the laminate section slice.";
            return result;
        }

        Vector target = candidate.segment(state.variableOffset, kSectionSize);
        Vector current = target;
        const double lowerThickness = state.thicknessLowerBound;
        const double upperThickness = std::max(state.thicknessUpperBound, lowerThickness + 2.0 * m_options.thicknessInteriorMargin);
        const int thicknessIndex = kSectionSize - 1;
        const double projectedThickness = clampThickness(target(thicknessIndex), lowerThickness, upperThickness);

        current(thicknessIndex) = projectedThickness;

        Section section;
        section.setBoundThickness(lowerThickness, upperThickness);
        section.Initialise(1.0, current);

        double dualityGap = std::numeric_limits<double>::infinity();
        int iteration = 0;
        do {
            Vector residual = -(current - target);
            Hessian hessian = Hessian::Identity();

            dualityGap = section.DualityGap();
            const double penalty = dualityGap / static_cast<double>(kSectionSize);
            const double predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();

            section.HessianEval(hessian);
            section.CalculateResiduals(current, residual, predictor * penalty);

            Eigen::LLT<Hessian> factorisation(hessian);
            if (factorisation.info() != Eigen::Success) {
                result.message = "Laminate section Hessian factorisation failed during predictor step.";
                return result;
            }

            Vector increment = residual;
            factorisation.solveInPlace(increment);
            section.UpdateIncrements(increment);

            const double correctedGap = section.DualityGap();
            const double corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(dualityGap, correctedGap);

            residual = -(current - target);
            section.CalculateResiduals(current, residual, corrector * penalty);
            increment = residual;
            factorisation.solveInPlace(increment);
            section.UpdateIncrements(increment);

            double primalStep = 1.0;
            double dualStep = 1.0;
            section.StepSize(primalStep, dualStep);
            section.UpdateVariables(primalStep, dualStep);
            current += primalStep * increment;
            current(thicknessIndex) = projectedThickness;
        } while (dualityGap > m_options.tolerance && ++iteration < m_options.maxIterations);

        result.iterations = iteration;
        result.finalGap = dualityGap;
        result.success = dualityGap <= m_options.tolerance;
        result.message = result.success
            ? "Laminate section projection converged."
            : "Laminate section projection reached the iteration limit.";

        if (result.success) {
            current(thicknessIndex) = projectedThickness;
            candidate.segment(state.variableOffset, kSectionSize) = current;
        }

        return result;
    }

    [[nodiscard]] double clampThickness(const double thickness,
                                        const double lowerThickness,
                                        const double upperThickness) const {
        const double lower = lowerThickness + m_options.thicknessInteriorMargin;
        const double upper = upperThickness - m_options.thicknessInteriorMargin;
        if (lower > upper) {
            return 0.5 * (lowerThickness + upperThickness);
        }
        return std::clamp(thickness, lower, upper);
    }
};

}  // namespace lamopt
