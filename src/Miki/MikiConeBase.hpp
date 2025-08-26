#pragma once

#include "../SDPA/SDPABase.hpp"
#include <Eigen/Dense>

namespace lampar {
namespace Miki {

template<typename ScalarType>
class MikiConeBase : public SDPA::SDPABase<ScalarType> {
public:
    MikiConeBase() = default;
    virtual ~MikiConeBase() = default;
    
    // Implementation of base class interface
    ScalarType DualityGap() override {
        return (m_X + m_dX).cwiseProduct(m_Y + m_dY).sum();
    }
    
    void CalculateResiduals(ScalarType mu, 
                          const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& r,
                          const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& RZ) override {
        // Implementation from Algorithm 2 internal procedures
        m_RZ = mu * Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Identity(m_X.rows(), m_X.cols())
             - m_X * m_Y - m_dX * m_dY;
           
        m_r = r;
        // Y should affect each constraint in M properly
        m_r.noalias() += m_M.transpose() * (m_Y * Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>::Ones(m_Y.cols()));
        m_dY = m_X.inverse() * m_RZ;
        // Update residual with the effect of dY on each constraint
        m_r.noalias() += m_M.transpose() * (m_dY * Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>::Ones(m_dY.cols()));
    }
    
    void UpdateIncrements(const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& G) override {
        // dx is the solution to the increment equation (m x 1)
        Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> dx = G.col(0).head(m_M.cols());
        
        // Initialize m_dX with zeros of correct size
        m_dX = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Zero(m_X.rows(), m_X.cols());
        
        // Update first column of m_dX with M * dx
        m_dX.col(0).noalias() = m_M * dx;
        
        // Calculate dY
        m_dY.noalias() = m_X.inverse() * (m_RZ - m_dX * m_Y);
    }
    
    void CalculateHessian(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& G,
                         Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& B) override {
        G = m_M;  // G is n x m
        // Compute B as X^{-1} * M to maintain n x m dimensions
        B.noalias() = m_X.inverse() * m_M;
    }
    
    void StepSize(ScalarType& alphaP, ScalarType& alphaD) override {
        // Compute minimum eigenvalue for step size
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> es;
        es.compute(m_dX);
        ScalarType lambda = es.eigenvalues().minCoeff();
        
        // Calculate step size based on Algorithm 2
        ScalarType alpha = std::max(-lambda, 1.0 - m_delta);
        alphaP = alphaD = (1.0 - m_delta) / alpha;
    }
    
    void UpdateVariables(ScalarType alphaP, ScalarType alphaD) override {
        m_X += alphaP * m_dX;
        m_dX.setZero();
        m_Y += alphaD * m_dY;
        m_dY.setZero();
    }
    
protected:
    // Problem matrices
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_X; // Primal matrix
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_Y; // Dual matrix
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_dX; // Primal increment
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_dY; // Dual increment
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_M; // Constraint matrix
    Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> m_RZ; // Residual matrix
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> m_r; // Residual vector
    
    // Parameters
    ScalarType m_delta = 0.1; // Step size parameter
};

} // namespace Miki
} // namespace lampar
