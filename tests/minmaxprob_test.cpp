#include "BoundSDP/boundSDP.hpp"
#include "GlobalOptimiser/approxFunction.hpp"
#include "GlobalOptimiser/scminmaxProb.hpp"

#include <gtest/gtest.h>
#include <cmath>

namespace {

using ApproximationFunction = globopt::approxFunction<4, 2>;
using SideConstraint = lampar::boundSDP<lampar::Lower>;

ApproximationFunction BuildProblem1Approximation(const ApproximationFunction::Vector_v& design) {
    ApproximationFunction approximationFunction;
    ApproximationFunction::Vector_r responses;
    ApproximationFunction::Matrix_t gradients = ApproximationFunction::Matrix_t::Zero();
    ApproximationFunction::Vector_r objectiveMask;
    objectiveMask << 1.0, 0.0, 0.0, 0.0;

    const double x1 = design(0);
    const double x2 = design(1);
    const double sqrt2 = std::sqrt(2.0);
    const double denominator = 2.0 * x1 * x2 + sqrt2 * x1 * x1;
    const double aux = sqrt2 * x1 * x2 + x1 * x1;
    const double common = sqrt2 * x2 + x1;
    const double objectiveScale = 100.0 * (2.0 * sqrt2 * x1 + x2);

    responses << 1.0,
                 (x2 + sqrt2 * x1) / denominator - 1.0,
                 1.0 / (x1 + sqrt2 * x2) - 1.0,
                 4.0 / 3.0 * x2 / denominator - 1.0;

    gradients.col(0) << 200.0 * sqrt2 / objectiveScale, 100.0 / objectiveScale;
    gradients.col(1) << -(sqrt2 * x1 * x2 + x1 * x1 + x2 * x2) / std::pow(aux, 2),
                         -1.0 / (sqrt2 * std::pow(common, 2));
    gradients.col(2) << -1.0 / std::pow(common, 2),
                         -sqrt2 / std::pow(common, 2);
    gradients.col(3) << -4.0 / 3.0 * x2 * (x2 + sqrt2 * x1) / std::pow(aux, 2),
                         4.0 / 3.0 * 1.0 / (sqrt2 * std::pow(common, 2));

    approximationFunction.ConfigureConLinModel(design, responses, gradients, objectiveMask, 1);
    return approximationFunction;
}

ApproximationFunction::Vector_v SolveProblem1() {
    SideConstraint sideConstraints[2];
    sideConstraints[0].setBound(0.0);
    sideConstraints[1].setBound(0.0);

    ApproximationFunction::Vector_v design;
    design << 1.0, 1.0;

    globopt::responseInterface<ApproximationFunction, SideConstraint> solver;

    for (int iteration = 0; iteration < 20; ++iteration) {
        ApproximationFunction approximationFunction = BuildProblem1Approximation(design);
        solver.scMinMaxProb(&approximationFunction, sideConstraints);
        solver.Solver(design);
    }

    return design;
}

}  // namespace

TEST(MinMaxProblemRegression, Problem1ConvergesToKnownSolution) {
    const ApproximationFunction::Vector_v design = SolveProblem1();

    EXPECT_NEAR(design(0), 0.788648, 1.0e-4);
    EXPECT_NEAR(design(1), 0.408326, 1.0e-4);

    ApproximationFunction approximationFunction = BuildProblem1Approximation(design);
    ApproximationFunction::Vector_r responses;
    approximationFunction.Eval(design, responses);

    EXPECT_NEAR(responses(0), 1.0, 1.0e-6);
    EXPECT_LT(responses(1), 0.0);
    EXPECT_LT(responses(2), 0.0);
    EXPECT_LT(responses(3), 0.0);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
