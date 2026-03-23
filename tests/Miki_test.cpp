#include "Miki/Miki.hpp"

#include <gtest/gtest.h>
#include <Eigen/Cholesky>

namespace {

Eigen::Vector4d SolveMikiObjective() {
    lampar::Miki<lampar::SingleMaterial, false> feasibility(1.0);

    Eigen::Vector4d gradient;
    gradient << 1.0, 0.0, 0.0, 0.0;

    Eigen::Vector4d design = Eigen::Vector4d::Zero();
    Eigen::Matrix4d hessian;
    Eigen::LLT<Eigen::Matrix4d> factorisation;

    double dualityGap = 0.0;
    int iteration = 0;
    do {
        Eigen::Vector4d increment = -gradient;
        hessian.setZero();

        dualityGap = feasibility.DualityGap();
        const double penalty = dualityGap / 3.0;
        const double predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
        feasibility.CalculateResiduals(increment, predictor * penalty);
        feasibility.Hessian(hessian);
        factorisation.compute(hessian);
        factorisation.solveInPlace(increment);
        feasibility.UpdateIncrements(increment);

        const double correctedGap = feasibility.DualityGap();
        const double corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(dualityGap, correctedGap);

        increment = -gradient;
        feasibility.CalculateResiduals(increment, corrector * penalty);
        factorisation.solveInPlace(increment);
        feasibility.UpdateIncrements(increment);

        double primalStep = 1.0;
        double dualStep = 1.0;
        feasibility.StepSize(primalStep, dualStep);
        feasibility.UpdateVariables(primalStep, dualStep);
        design += primalStep * increment;
    } while (dualityGap > 1.0e-10 && ++iteration < 20);

    return design;
}

}  // namespace

TEST(MikiRegression, OptimisationOverTheMikiConeConvergesToKnownPoint) {
    const Eigen::Vector4d design = SolveMikiObjective();

    EXPECT_NEAR(design(0), -std::sqrt(2.0), 1.0e-5);
    EXPECT_NEAR(design(1), 1.0, 1.0e-5);
    EXPECT_NEAR(design(2), 0.0, 1.0e-7);
    EXPECT_NEAR(design(3), 0.0, 1.0e-7);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
