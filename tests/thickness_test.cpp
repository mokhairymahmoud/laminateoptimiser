#include "BoundSDP/boundSDP.hpp"

#include <gtest/gtest.h>

namespace {

double SolveThicknessBoundProblem() {
    lampar::boundSDP<lampar::Upper> upperBound(100.0);
    lampar::boundSDP<lampar::Lower> lowerBound(-100.0);

    double design = -40.0;
    upperBound.Initialise(1.0, design);
    lowerBound.Initialise(1.0, design);

    double dualityGap = 0.0;
    int iteration = 0;
    do {
        const double response = design * design - design - 2.0;
        double increment = -response;
        double hessian = 2.0 * design - 1.0;

        dualityGap = upperBound.DualityGap() + lowerBound.DualityGap();
        const double penalty = dualityGap / 2.0;
        const double predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
        upperBound.CalculateResiduals(increment, predictor * penalty);
        lowerBound.CalculateResiduals(increment, predictor * penalty);
        upperBound.Hessian(hessian);
        lowerBound.Hessian(hessian);
        increment /= hessian;
        upperBound.UpdateIncrements(increment);
        lowerBound.UpdateIncrements(increment);

        const double correctedGap = upperBound.DualityGap() + lowerBound.DualityGap();
        const double corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(dualityGap, correctedGap);

        increment = -response;
        upperBound.CalculateResiduals(increment, corrector * penalty);
        lowerBound.CalculateResiduals(increment, corrector * penalty);
        increment /= hessian;
        upperBound.UpdateIncrements(increment);
        lowerBound.UpdateIncrements(increment);

        double primalStep = 1.0;
        double dualStep = 1.0;
        upperBound.StepSize(primalStep, dualStep);
        lowerBound.StepSize(primalStep, dualStep);
        upperBound.UpdateVariables(primalStep, dualStep);
        lowerBound.UpdateVariables(primalStep, dualStep);
        design += primalStep * increment;
    } while (dualityGap > 1.0e-10 && ++iteration < 20);

    return design;
}

}  // namespace

TEST(ThicknessRegression, BoundConstrainedScalarProblemConvergesToKnownRoot) {
    EXPECT_NEAR(SolveThicknessBoundProblem(), -1.0, 1.0e-6);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
