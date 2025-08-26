#pragma once
#include <Eigen/Dense>

namespace lampar {
namespace GlobalOptimiser {

template<typename ScalarType>
ScalarType DampingMultiply(const ScalarType xi, const ScalarType rhoL) {
    ScalarType eta;
    if (xi < 0) {
        eta = 1.0 + (1.0 - rhoL) * std::tanh(xi) / (1.0 - rhoL);
    } else {
        eta = 1.0 + xi;
    }
    return eta;
}

template<typename ScalarType>
void UpdateDamping(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& f, 
                  const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& f_tilde,
                  const ScalarType d,
                  Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& rho) {
    const size_t n = f.size();
    for(size_t i = 0; i < n; ++i) {
        ScalarType xi = (f[i] - f_tilde[i]) / d;
        ScalarType eta = DampingMultiply(xi, rho[i]);
        rho[i] = eta * rho[i];
    }
}

template<typename ScalarType>
void DampingVector(const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& xref,
                  const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& x,
                  Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>& d) {
    d = (x - xref).cwiseAbs();
}

} // namespace GlobalOptimiser
} // namespace lampar
