/*
 * fSDPCone.hpp
 *
 *  Created on: 15 avr. 2015
 *      Author: Ramzi
 */

#ifndef FSDPCONE_FSDPCONE_HPP_
#define FSDPCONE_FSDPCONE_HPP_

#include <Eigen/Dense>
#include <vector>
#include "../SDPA/SDPA.hpp"

using namespace std;

namespace lampar {

template<int Dim, int Size, typename ScalarType = double>
class fSDPCone : public
      SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>
{
public:
	typedef typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>	Base_t;
	typedef typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>::Matrix_t Matrix_t;
	typedef typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>::Vector_t Vector_t;
	typedef typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>::Hessian_t Hessian_t;
private:
	Matrix_t M0;
	vector<Matrix_t> M;
public:
	fSDPCone() {};
	fSDPCone(const Matrix_t &convSetM0, const vector<Matrix_t> &convSetM) {
		InitConvexSet(convSetM0, convSetM);
	};
	void InitConvexSet(const Matrix_t &convSetM0, const vector<Matrix_t> &convSetM) {
		M0 = convSetM0;
		M = convSetM;
	}
	~fSDPCone() {
		M.clear();
	};
	void Initval(const Vector_t &var, Matrix_t &primal) {
		Eval(var,primal);
		primal += M0;
	}
	void Eval(const Vector_t &var, Matrix_t &primal) const {
		primal = Matrix_t::Zero();
		for(int i=0;i<Size;++i)
			primal += var(i)*M[i];
	}
	void AdEval(const Matrix_t &dual, Vector_t &residual) const {
		Matrix_t temp;
		for(int i=0;i<Size;++i){
			temp.noalias() = M[i]*dual;
			residual(i) += temp.trace();
		}
	}
	void HessianEval(const Matrix_t &primal, const Matrix_t &dual, Hessian_t &hessian) const {
		Matrix_t temp,temp_primal,temp_dual;
		for (int jSize=0;jSize<Size;++jSize)
			for (int iSize=jSize;iSize<Size;++iSize){
				temp_primal.noalias() = M[jSize]*primal;
				temp_dual.noalias() = M[iSize]*dual;
				temp.noalias() = temp_primal*temp_dual;
				hessian(jSize,iSize) = hessian(iSize,jSize) +=temp.trace();
			}
	}
};

/*template<int Dim, int Size, typename ScalarType>
typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>::Matrix_t lampar::fSDPCone<Dim, Size, ScalarType>::M0;
template<int Dim, int Size, typename ScalarType>
vector<typename SDPA::Constraint<fSDPCone<Dim, Size, ScalarType>, Dim, Size, ScalarType>::Matrix_t> lampar::fSDPCone<Dim, Size, ScalarType>::M;*/

}  /* namespace lampar */




#endif /* FSDPCONE_FSDPCONE_HPP_ */
