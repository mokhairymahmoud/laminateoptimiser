#include "../src/Miki/MikiConeBase.hpp"
#include <gtest/gtest.h>
#include <Eigen/Dense>

class MockMikiCone : public lampar::Miki::MikiConeBase<double> {
public:
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    // Expose protected members for testing
    void setTestMatrices(const Matrix& X, const Matrix& Y, 
                        const Matrix& dX, const Matrix& dY,
                        const Matrix& M) {
        this->m_X = X;
        this->m_Y = Y;
        this->m_dX = dX;
        this->m_dY = dY;
        this->m_M = M;
    }

    const Matrix& getX() const { return this->m_X; }
    const Matrix& getY() const { return this->m_Y; }
    const Matrix& getdX() const { return this->m_dX; }
    const Matrix& getdY() const { return this->m_dY; }
    const Matrix& getM() const { return this->m_M; }
    const Vector& getR() const { return this->m_r; }
    
    using Base = lampar::SDPA::SDPABase<double>;
    friend class MikiConeTest;
};

class MikiConeTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test matrices
        const int n = 3;  // Size of square matrices X, Y
        const int m = 2;  // Number of constraints (rows in M)
        
        X = Matrix::Identity(n, n);
        Y = Matrix::Identity(n, n);
        dX = Matrix::Zero(n, n);
        dY = Matrix::Zero(n, n);
        M = Matrix::Random(n, m);  // M is n x m for m constraints
        
        cone.setTestMatrices(X, Y, dX, dY, M);
    }

    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    
    MockMikiCone cone;
    Matrix X, Y, dX, dY, M;
    const double tolerance = 1e-10;
};

TEST_F(MikiConeTest, DualityGap) {
    // Test duality gap calculation
    double gap = cone.DualityGap();
    double expected = (X + dX).cwiseProduct(Y + dY).sum();
    EXPECT_NEAR(gap, expected, tolerance);
}

TEST_F(MikiConeTest, CalculateResiduals) {
    // Test residuals calculation
    const double mu = 0.1;
    Vector r = Vector::Random(2);  // r has length m (number of constraints)
    Matrix G = Matrix::Random(3, 2);  // G is n x m
    
    cone.CalculateResiduals(mu, r, G);
    
    // The internal m_RZ should be calculated as per the formula
    Matrix expectedRZ = mu * Matrix::Identity(3, 3) - X * Y - dX * dY;
    Matrix expectedDY = X.inverse() * expectedRZ;
    
    // Calculate expected residual vector
    Vector expectedR = r;
    expectedR.noalias() += M.transpose() * (Y * Vector::Ones(Y.cols()));
    expectedR.noalias() += M.transpose() * (expectedDY * Vector::Ones(expectedDY.cols()));
    
    // Verify dY calculation and residual update
    EXPECT_TRUE(cone.getdY().isApprox(expectedDY, tolerance));
    EXPECT_TRUE(cone.getR().isApprox(expectedR, tolerance));
}

TEST_F(MikiConeTest, UpdateIncrements) {
    Matrix G = Matrix::Random(3, 2);  // G is n x m
    Vector dx = G.col(0).head(M.cols());  // Take only the first m components
    
    // First set up RZ
    const double mu = 0.1;
    Vector r = Vector::Random(2);  // r has length m
    cone.CalculateResiduals(mu, r, G);
    
    // Now test increment update
    cone.UpdateIncrements(G);
    
    // Verify dX calculation
    Matrix expectedDX = Matrix::Zero(X.rows(), X.cols());  // Same size as X
    expectedDX.col(0).noalias() = M * dx;  // First column is M * dx
    EXPECT_TRUE(cone.getdX().isApprox(expectedDX, tolerance));
}

TEST_F(MikiConeTest, CalculateHessian) {
    Matrix G = Matrix::Zero(3, 2);  // G should be n x m
    Matrix B = Matrix::Zero(3, 2);  // B should be n x m
    
    cone.CalculateHessian(G, B);
    
    // Verify G matrix matches M
    EXPECT_TRUE(G.isApprox(M, tolerance));
    
    // B should be X^{-1} * M
    Matrix expectedB = X.inverse() * M;
    EXPECT_TRUE(B.isApprox(expectedB, tolerance));
}

TEST_F(MikiConeTest, StepSize) {
    double alphaP = 0.0;
    double alphaD = 0.0;
    
    cone.StepSize(alphaP, alphaD);
    
    // Step sizes should be positive and less than or equal to 1
    EXPECT_GT(alphaP, 0.0);
    EXPECT_LE(alphaP, 1.0);
    EXPECT_GT(alphaD, 0.0);
    EXPECT_LE(alphaD, 1.0);
}

TEST_F(MikiConeTest, UpdateVariables) {
    double alphaP = 0.5;
    double alphaD = 0.3;
    
    Matrix oldX = cone.getX();
    Matrix oldY = cone.getY();
    Matrix olddX = cone.getdX();
    Matrix olddY = cone.getdY();
    
    cone.UpdateVariables(alphaP, alphaD);
    
    // Verify X update
    EXPECT_TRUE(cone.getX().isApprox(oldX + alphaP * olddX, tolerance));
    EXPECT_TRUE(cone.getdX().isZero(tolerance));
    
    // Verify Y update
    EXPECT_TRUE(cone.getY().isApprox(oldY + alphaD * olddY, tolerance));
    EXPECT_TRUE(cone.getdY().isZero(tolerance));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
