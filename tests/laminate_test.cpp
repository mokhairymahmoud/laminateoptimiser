#include "Laminate/lpfeasible.hpp"
#include "Miki/Miki.hpp"

#include <gtest/gtest.h>
#include <Eigen/Cholesky>

namespace {

constexpr int kSubLaminates = 6;

using Laminate = lampar::lpfeasible<lampar::SingleMaterial, false, true, kSubLaminates, double>;

Laminate::Vector_t SolveLaminateFeasibilityCase() {
    Laminate laminate;
    Laminate::Hessian_t hessian;
    Laminate::Cholesky_t factorisation;

    Laminate::Vector_t design = Laminate::Vector_t::Zero();
    Laminate::Vector_t increment;
    Laminate::Vector_t residual;
    Laminate::Vector_t gradient;
    gradient << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    laminate.Initialise(1.0);

    double dualityGap = 0.0;
    int iteration = 0;
    do {
        residual = -gradient;
        hessian.setZero();

        dualityGap = laminate.DualityGap();
        const double penalty = dualityGap / static_cast<double>(kSubLaminates * 3);
        const double predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
        laminate.HessianEval(hessian);
        laminate.CalculateResiduals(design, residual, predictor * penalty);
        factorisation.compute(hessian);
        increment = residual;
        factorisation.solveInPlace(increment);
        laminate.UpdateIncrements(increment);

        const double correctedGap = laminate.DualityGap();
        const double corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(dualityGap, correctedGap);

        residual = -gradient;
        laminate.CalculateResiduals(design, residual, corrector * penalty);
        increment = residual;
        factorisation.solveInPlace(increment);
        laminate.UpdateIncrements(increment);

        double primalStep = 1.0;
        double dualStep = 1.0;
        laminate.StepSize(primalStep, dualStep);
        laminate.UpdateVariables(primalStep, dualStep);
        design += primalStep * increment;
    } while (dualityGap > 1.0e-10 && ++iteration < 20);

    return design;
}

}  // namespace

TEST(LaminateRegression, FeasibilitySolveMatchesReferenceDirectionCase) {
    const Laminate::Vector_t design = SolveLaminateFeasibilityCase();

    EXPECT_NEAR(design(0), 0.0, 1.0e-9);
    EXPECT_NEAR(design(1), -1.0, 1.0e-6);
    EXPECT_NEAR(design(2), -std::sqrt(2.0), 1.0e-5);
    EXPECT_NEAR(design(3), 0.0, 1.0e-9);
    EXPECT_NEAR(design(4), 0.0, 1.0e-9);
    EXPECT_NEAR(design(5), -1.0, 1.0e-6);
    EXPECT_NEAR(design(6), -std::sqrt(2.0), 1.0e-5);
    EXPECT_NEAR(design(7), 0.0, 1.0e-8);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
