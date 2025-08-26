#pragma once

#include "SDPABase.hpp"
#include <memory>

namespace lampar {
namespace SDPA {

template<typename ScalarType>
class SDPASolver {
public:
    explicit SDPASolver(std::unique_ptr<SDPABase<ScalarType>> sdpa)
        : m_sdpa(std::move(sdpa))
        , m_tolerance(1e-6) {}
    
    void Solve(Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& x) {
        if (!m_sdpa) return;

        // Initialize dimensions
        n = x.size();
        m = n; // For this solver, constraints match variables
        
        // Get initial duality gap
        ScalarType delta = m_sdpa->DualityGap();
        
        // Main solver loop
        while(delta > m_tolerance) {
            // Get objective function values
            auto [f, G, B] = GetObjectiveFunctionValues(x);
            
            // Create objective vector for permutation
            Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> o = 
                Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>::Ones(n);
            
            // Calculate barrier parameter
            ScalarType mu = delta / (n + m);
            
            // Calculate Hessian
            m_sdpa->CalculateHessian(G, B);
            
            // Predictor step
            m_sdpa->CalculateResiduals(m_options.betaP * mu, f, G);
            m_sdpa->UpdateIncrements(G);
            
            // Get predictor duality gap
            ScalarType deltaTilde = m_sdpa->DualityGap();
            
            // Corrector step
            ScalarType betaC = deltaTilde * deltaTilde / (delta * delta);
            betaC = std::max(betaC, m_options.betaP);
            
            m_sdpa->CalculateResiduals(betaC * mu, f, G);
            m_sdpa->UpdateIncrements(G);
            
            // Calculate step sizes
            ScalarType alphaP, alphaD;
            m_sdpa->StepSize(alphaP, alphaD);
            
            // Update variables
            m_sdpa->UpdateVariables(alphaP, alphaD);
            
            // Update duality gap
            delta = m_sdpa->DualityGap();
        }
    }
    
private:
    // Helper function to get objective function values
    auto GetObjectiveFunctionValues(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& x) {
        struct Result {
            Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> f;
            Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> G, B;
        };
        
        Result result;
        result.f = x; // Simple objective for testing
        result.G = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Identity(n, n);
        result.B = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);
        
        return result;
    }
    
    void SetTolerance(ScalarType tol) {
        m_tolerance = tol;
    }
    
private:
    std::unique_ptr<SDPABase<ScalarType>> m_sdpa;
    ScalarType m_tolerance;
    
    struct Options {
        ScalarType betaP = 0.5; // Default predictor coefficient
    } m_options;
    
    // Problem dimensions
    int n; // Number of variables
    int m; // Number of constraints
};

} // namespace SDPA
} // namespace lampar
