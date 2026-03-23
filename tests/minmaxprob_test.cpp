#include "BoundSDP/boundSDP.hpp"
#include "GlobalOptimiser/approxFunction.hpp"
#include "GlobalOptimiser/scminmaxProb.hpp"

#include <gtest/gtest.h>

namespace {

using ApproximationFunction = globopt::approxFunction<4, 2>;
using SideConstraint = lampar::boundSDP<lampar::Lower>;

ApproximationFunction::Vector_v SolveProblem1() {
    ApproximationFunction approximationFunction;
    SideConstraint sideConstraints[2];
    sideConstraints[0].setBound(0.0);
    sideConstraints[1].setBound(0.0);

    ApproximationFunction::Vector_v design;
    design << 1.0, 1.0;

    globopt::responseInterface<ApproximationFunction, SideConstraint> solver;
    solver.scMinMaxProb(&approximationFunction, sideConstraints);

    for (int iteration = 0; iteration < 20; ++iteration) {
        approximationFunction.InitialiseProblem1(design);
        solver.Solver(design);
    }

    return design;
}

}  // namespace

TEST(MinMaxProblemRegression, Problem1ConvergesToKnownSolution) {
    const ApproximationFunction::Vector_v design = SolveProblem1();

    EXPECT_NEAR(design(0), 0.788648, 1.0e-4);
    EXPECT_NEAR(design(1), 0.408326, 1.0e-4);

    ApproximationFunction approximationFunction;
    ApproximationFunction::Vector_r responses;
    approximationFunction.InitialiseProblem1(design);
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
