#include "../src/GlobalOptimiser/dampingUtils.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace lampar::GlobalOptimiser;

class DampingUtilsTest : public ::testing::Test {
protected:
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    
    void SetUp() override {
        // Set up test vectors
        f = Vector::Zero(3);
        f_tilde = Vector::Zero(3);
        rho = Vector::Ones(3);
        x = Vector::Zero(3);
        xref = Vector::Zero(3);
        d = Vector::Zero(3);
    }

    Vector f, f_tilde, rho, x, xref, d;
    const double tolerance = 1e-10;
};

TEST_F(DampingUtilsTest, DampingMultiply_NegativeXi) {
    // Test when xi < 0
    double xi = -1.0;
    double rhoL = 0.5;
    double expected = 1.0 + (1.0 - rhoL) * std::tanh(xi) / (1.0 - rhoL);
    
    double result = DampingMultiply(xi, rhoL);
    
    EXPECT_NEAR(result, expected, tolerance);
}

TEST_F(DampingUtilsTest, DampingMultiply_PositiveXi) {
    // Test when xi >= 0
    double xi = 1.0;
    double rhoL = 0.5;
    double expected = 1.0 + xi;
    
    double result = DampingMultiply(xi, rhoL);
    
    EXPECT_NEAR(result, expected, tolerance);
}

TEST_F(DampingUtilsTest, DampingMultiply_ZeroXi) {
    // Test when xi = 0
    double xi = 0.0;
    double rhoL = 0.5;
    double expected = 1.0; // Should be 1.0 + 0.0
    
    double result = DampingMultiply(xi, rhoL);
    
    EXPECT_NEAR(result, expected, tolerance);
}

TEST_F(DampingUtilsTest, UpdateDamping_AllPositive) {
    // Test with all positive differences
    f << 2.0, 3.0, 4.0;
    f_tilde << 1.0, 1.0, 1.0;
    double d_scalar = 1.0;
    
    Vector expected = Vector::Ones(3);
    for(int i = 0; i < 3; ++i) {
        double xi = (f[i] - f_tilde[i]) / d_scalar;
        expected[i] *= (1.0 + xi);
    }
    
    UpdateDamping(f, f_tilde, d_scalar, rho);
    
    for(int i = 0; i < 3; ++i) {
        EXPECT_NEAR(rho[i], expected[i], tolerance);
    }
}

TEST_F(DampingUtilsTest, UpdateDamping_MixedDifferences) {
    // Test with mixed positive and negative differences
    f << -1.0, 0.0, 1.0;
    f_tilde << 1.0, 0.0, -1.0;
    double d_scalar = 1.0;
    rho.setConstant(0.5); // Set all rho values to 0.5
    
    Vector expected = Vector::Constant(3, 0.5);
    for(int i = 0; i < 3; ++i) {
        double xi = (f[i] - f_tilde[i]) / d_scalar;
        double eta = DampingMultiply(xi, expected[i]);
        expected[i] *= eta;
    }
    
    UpdateDamping(f, f_tilde, d_scalar, rho);
    
    for(int i = 0; i < 3; ++i) {
        EXPECT_NEAR(rho[i], expected[i], tolerance);
    }
}

TEST_F(DampingUtilsTest, DampingVector_ZeroDistance) {
    // Test when x = xref
    x.setZero();
    xref.setZero();
    
    DampingVector(xref, x, d);
    
    EXPECT_TRUE(d.isZero(tolerance));
}

TEST_F(DampingUtilsTest, DampingVector_UnitDistance) {
    // Test with unit distances
    x << 1.0, 2.0, 3.0;
    xref << 0.0, 1.0, 2.0;
    Vector expected = Vector::Constant(3, 1.0);
    
    DampingVector(xref, x, d);
    
    for(int i = 0; i < 3; ++i) {
        EXPECT_NEAR(d[i], expected[i], tolerance);
    }
}

TEST_F(DampingUtilsTest, DampingVector_NegativeDistance) {
    // Test with negative differences
    x << -1.0, -2.0, -3.0;
    xref << 1.0, 2.0, 3.0;
    Vector expected(3);
    expected << 2.0, 4.0, 6.0;  // |(-1 - 1)|, |(-2 - 2)|, |(-3 - 3)|
    
    DampingVector(xref, x, d);
    
    for(int i = 0; i < 3; ++i) {
        EXPECT_NEAR(d[i], expected[i], tolerance);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
