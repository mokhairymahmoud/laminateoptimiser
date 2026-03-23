#include "fSDPCone/fSDPCone.hpp"

#include <gtest/gtest.h>
#include <Eigen/Cholesky>

namespace {

using Matrix3 = Eigen::Matrix<double, 3, 3>;

lampar::fSDPCone<3, 4> BuildCone() {
    Matrix3 freeTerm;
    freeTerm.setIdentity();

    std::vector<Matrix3> matrices;
    Matrix3 basis;

    basis << 0.0, 1.0, 0.0,
             1.0, 0.0, 0.0,
             0.0, 0.0, 0.0;
    matrices.push_back(basis);

    basis << 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, -1.0;
    matrices.push_back(basis);

    basis << 0.0, 0.0, 1.0,
             0.0, 0.0, 0.0,
             1.0, 0.0, 0.0;
    matrices.push_back(basis);

    basis << 0.0, 0.0, 0.0,
             0.0, 0.0, 1.0,
             0.0, 1.0, 0.0;
    matrices.push_back(basis);

    lampar::fSDPCone<3, 4> feasibility;
    feasibility.InitConvexSet(freeTerm, matrices);
    return feasibility;
}

Eigen::Vector4d SolveConeProblem() {
    lampar::fSDPCone<3, 4> feasibility = BuildCone();
    Eigen::Vector4d gradient;
    gradient << 1.0, 0.0, 0.0, 0.0;

    Eigen::Vector4d design = Eigen::Vector4d::Zero();
    feasibility.Initialise(1.0, design);

    Eigen::Matrix4d hessian;
    Eigen::LLT<Eigen::Matrix4d> factorisation;
    double dualityGap = 0.0;
    int iteration = 0;
    do {
        Eigen::Vector4d increment = -gradient;
        hessian.setZero();

        dualityGap = feasibility.DualityGap();
        const double penalty = dualityGap / 3.0;
        feasibility.CalculateResiduals(increment, SDPA::Parameter<>::Predictor_Duality_Reduction() * penalty);
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

TEST(FSDPConeRegression, GenericConeMatchesReferenceMikiSolution) {
    const Eigen::Vector4d design = SolveConeProblem();

    EXPECT_NEAR(design(0), -std::sqrt(2.0), 1.0e-5);
    EXPECT_NEAR(design(1), 1.0, 1.0e-5);
    EXPECT_NEAR(design(2), 0.0, 1.0e-7);
    EXPECT_NEAR(design(3), 0.0, 1.0e-7);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
