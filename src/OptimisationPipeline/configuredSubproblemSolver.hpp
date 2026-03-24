#pragma once

#include "coreMinMaxSubproblem.hpp"
#include "defaultLaminateSubproblem.hpp"

#include <memory>

namespace lamopt {

enum class BaseSubproblemKind {
    GradientPenalty,
    CoreMinMax2Var4Resp
};

struct ConfiguredSubproblemOptions {
    BaseSubproblemKind baseSubproblemKind = BaseSubproblemKind::GradientPenalty;
    GradientPenaltySubproblemOptions gradientPenaltyOptions;
    CoreMinMax2Var4RespOptions coreMinMaxOptions;
    CoreLaminateSection1RespOptions coreLaminateOptions;
    LaminateProjectionOptions laminateProjectionOptions;
};

class ConfiguredSubproblemSolver final : public SubproblemSolver {
public:
    explicit ConfiguredSubproblemSolver(ConfiguredSubproblemOptions options = {})
        : m_gradientPenaltySolver(options.gradientPenaltyOptions)
        , m_coreMinMaxSolver(options.coreMinMaxOptions)
        , m_directLaminateSolver(options.coreLaminateOptions)
        , m_routedSolver(selectBaseSolver(options.baseSubproblemKind,
                                          m_gradientPenaltySolver,
                                          m_coreMinMaxSolver),
                         m_directLaminateSolver,
                         options.laminateProjectionOptions) {}

    SubproblemResult solve(const ApproximationProblem& problem) override {
        return m_routedSolver.solve(problem);
    }

private:
    [[nodiscard]] static SubproblemSolver&
    selectBaseSolver(const BaseSubproblemKind kind,
                     GradientPenaltySubproblemSolver& gradientPenaltySolver,
                     CoreMinMax2Var4RespSubproblemSolver& coreMinMaxSolver) {
        switch (kind) {
            case BaseSubproblemKind::GradientPenalty:
                return gradientPenaltySolver;
            case BaseSubproblemKind::CoreMinMax2Var4Resp:
                return coreMinMaxSolver;
        }
        return gradientPenaltySolver;
    }

    GradientPenaltySubproblemSolver m_gradientPenaltySolver;
    CoreMinMax2Var4RespSubproblemSolver m_coreMinMaxSolver;
    CoreLaminateSection1RespSubproblemSolver m_directLaminateSolver;
    DefaultLaminateSubproblemSolver m_routedSolver;
};

inline std::unique_ptr<SubproblemSolver>
MakeConfiguredSubproblemSolver(const ConfiguredSubproblemOptions& options = {}) {
    return std::make_unique<ConfiguredSubproblemSolver>(options);
}

}  // namespace lamopt
