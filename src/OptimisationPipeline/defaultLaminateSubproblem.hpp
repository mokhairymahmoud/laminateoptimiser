#pragma once

#include "laminateAwareSubproblem.hpp"

#include <string>

namespace lamopt {

class DefaultLaminateSubproblemSolver : public SubproblemSolver {
public:
    DefaultLaminateSubproblemSolver(SubproblemSolver& defaultSolver,
                                    SubproblemSolver& directLaminateSolver,
                                    LaminateProjectionOptions projectionOptions = {})
        : m_defaultSolver(defaultSolver)
        , m_directLaminateSolver(directLaminateSolver)
        , m_projectionFallback(defaultSolver, projectionOptions) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        if (problem.laminateSections.empty()) {
            return m_defaultSolver.solve(problem);
        }

        SubproblemResult directResult = m_directLaminateSolver.solve(problem);
        if (directResult.success) {
            directResult.message += " | route: direct_core_laminate";
            return directResult;
        }

        SubproblemResult fallbackResult = m_projectionFallback.solve(problem);
        if (fallbackResult.success) {
            fallbackResult.message =
                "Direct core laminate route unavailable: " + directResult.message
                + " | route: laminate_projection_fallback | " + fallbackResult.message;
            return fallbackResult;
        }

        fallbackResult.message =
            "Direct core laminate route unavailable: " + directResult.message
            + " | laminate projection fallback failed: " + fallbackResult.message;
        return fallbackResult;
    }

private:
    SubproblemSolver& m_defaultSolver;
    SubproblemSolver& m_directLaminateSolver;
    LaminateAwareSubproblemSolver m_projectionFallback;
};

}  // namespace lamopt
