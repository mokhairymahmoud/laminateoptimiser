#pragma once

#include <Eigen/Dense>
#include <vector>

namespace lampar {
namespace SDPA {

template<typename ScalarType>
class SDPABase {
public:
    virtual ~SDPABase() = default;
    
    // Required interface methods from pseudocode
    virtual ScalarType DualityGap() = 0;
    virtual void CalculateResiduals(ScalarType mu, const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& f,
                                  const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& G) = 0;
    virtual void UpdateIncrements(const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& G) = 0;
    virtual void CalculateHessian(Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& G,
                                Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& B) = 0;
    virtual void StepSize(ScalarType& alphaP, ScalarType& alphaD) = 0;
    virtual void UpdateVariables(ScalarType alphaP, ScalarType alphaD) = 0;
    
protected:
    // Common utilities from pseudocode
    ScalarType CorrectorMultiplier(ScalarType delta, ScalarType deltaTilde) {
        ScalarType rho = std::min(delta, deltaTilde) / delta;
        return std::max(rho * rho, options.betaC);
    }
    
    void CreatePermutationMatrix(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& f,
                               const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& o,
                               Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>& P) {
        const int n = f.size();
        std::vector<int> index(n);
        std::vector<int> indexPrime(n);
        
        // Initialize indices
        for(int i = 0; i < n; ++i) {
            indexPrime[i] = i;
        }
        
        // Sort indices based on f values
        shellSort(f, indexPrime);
        
        // Find critical index
        int ic = 0;
        while(o[indexPrime[ic]] != 1 && ic < n) {
            ++ic;
        }
        
        // Create final index array
        for(int i = 0; i < ic; ++i) {
            index[i] = indexPrime[i];
        }
        for(int i = ic + 1; i < n; ++i) {
            index[i-1] = indexPrime[i];
        }
        index[n-1] = indexPrime[ic];
        
        // Create permutation matrix
        P.setZero(n, n);
        for(int i = 0; i < n; ++i) {
            P(i, index[i]) = 1;
        }
    }
    
private:
    void shellSort(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& f,
                  std::vector<int>& index) {
        const int n = f.size();
        for(int gap = n/2; gap > 0; gap /= 2) {
            for(int i = gap; i < n; ++i) {
                int tempIndex = index[i];
                ScalarType tempValue = f[index[i]];
                int j;
                for(j = i; j >= gap && f[index[j-gap]] > tempValue; j -= gap) {
                    index[j] = index[j-gap];
                }
                index[j] = tempIndex;
            }
        }
    }
    
    struct Options {
        ScalarType betaC = 0.1; // Default corrector multiplier minimum
        // Add other solver options as needed
    } options;
};

} // namespace SDPA
} // namespace lampar
