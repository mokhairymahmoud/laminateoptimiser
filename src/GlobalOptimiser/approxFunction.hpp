/*
 * approxFunction.hpp
 *
 *  Created on: 3 oct. 2015
 *      Author: Ramzi
 */

#ifndef GLOBALOPTIMISER_APPROXFUNCTION_HPP_
#define GLOBALOPTIMISER_APPROXFUNCTION_HPP_

#include <Eigen/Dense>
#include <math.h>

namespace globopt {
template <int nResp, int nVar, typename ScalarType = double>
class approxFunction {
public:
	int NVAR;
	int NRESP;
	typedef Eigen::Matrix<ScalarType,nResp,1> Vector_r;
	typedef Eigen::Matrix<ScalarType,nVar,1> Vector_v;
	typedef Eigen::Matrix<ScalarType,nVar,nResp> Matrix_t;
	typedef Eigen::Matrix<ScalarType,nVar,nVar> Hessian_t;
	typedef Eigen::Matrix<ScalarType,nResp,nResp> Hessian_r;
	typedef Eigen::Matrix<int,nResp,nResp> Matrix_r_i;
	typedef Eigen::Matrix<int,nResp,1> Vector_r_i;

	typedef Eigen::Matrix<ScalarType,nResp-1,nResp-1> reduced_Hessian_r;
	typedef Eigen::Matrix<ScalarType,nResp-1,1> reduced_Vector_r;

private:
	Vector_r free_terms;
	Matrix_t sensitivitiesLinear;
	Matrix_t sensitivitiesReciprocal;
	Vector_r booleanVector;
	int dimObjSet;

private:
	void LinearOrReciprocal(Vector_v primalVar, ScalarType resp, Vector_v grad, ScalarType &freeterm, Vector_v &gradLinear, Vector_v &gradReciprocal){
		int iVar;
		Vector_v gradTemp;

		gradTemp = grad.cwiseProduct(primalVar);
		freeterm = resp;
		for (iVar = 0;iVar<nVar;iVar++){
			gradLinear(iVar) = 0;
			gradReciprocal(iVar) = 0;
			if (gradTemp(iVar)<0){
				gradReciprocal(iVar) = -primalVar(iVar)*gradTemp(iVar);
				freeterm += gradTemp(iVar);
			}
			else {
				gradLinear(iVar) = grad(iVar);
				freeterm -= gradTemp(iVar);
			}
		}
	};

public:
	approxFunction() {
		free_terms = Vector_r::Zero();
		sensitivitiesLinear = Matrix_t::Zero();
		sensitivitiesReciprocal = Matrix_t::Zero();
		booleanVector = Vector_r::Zero();
		dimObjSet = 0;
		NVAR = nVar;
		NRESP = nResp;
	};
	virtual ~approxFunction() {};
	int getBooleanVector(Vector_r &booleanVect) {
		booleanVect = booleanVector;
		return dimObjSet;
	}
	void ConfigureLinearModel(const Vector_v& referenceDesign,
	                         const Vector_r& referenceResponses,
	                         const Matrix_t& linearGradients,
	                         const Vector_r& objectiveMask,
	                         const int objectiveCount) {
		int iResp;

		free_terms = referenceResponses;
		sensitivitiesLinear = linearGradients;
		sensitivitiesReciprocal = Matrix_t::Zero();
		booleanVector = objectiveMask;
		dimObjSet = objectiveCount;

		for (iResp = 0; iResp < nResp; ++iResp) {
			free_terms(iResp) -= linearGradients.col(iResp).transpose() * referenceDesign;
		}
	}
	void ConfigureConLinModel(const Vector_v& referenceDesign,
	                          const Vector_r& referenceResponses,
	                          const Matrix_t& responseGradients,
	                          const Vector_r& objectiveMask,
	                          const int objectiveCount) {
		int iResp;
		ScalarType freeTerm;
		Vector_v gradLinear, gradReciprocal;

		booleanVector = objectiveMask;
		dimObjSet = objectiveCount;
		for (iResp = 0; iResp < nResp; ++iResp) {
			LinearOrReciprocal(referenceDesign,
			                  referenceResponses(iResp),
			                  responseGradients.col(iResp),
			                  freeTerm,
			                  gradLinear,
			                  gradReciprocal);
			free_terms(iResp) = freeTerm;
			sensitivitiesLinear.col(iResp) = gradLinear;
			sensitivitiesReciprocal.col(iResp) = gradReciprocal;
		}
	}
	void ConfigureModel(const Vector_r& modelFreeTerms,
	                    const Matrix_t& linearGradients,
	                    const Matrix_t& reciprocalGradients,
	                    const Vector_r& objectiveMask,
	                    const int objectiveCount) {
		free_terms = modelFreeTerms;
		sensitivitiesLinear = linearGradients;
		sensitivitiesReciprocal = reciprocalGradients;
		booleanVector = objectiveMask;
		dimObjSet = objectiveCount;
	}
	void Eval(Vector_v primalVar, Vector_r dualVar, Vector_r &responses, Matrix_t &gradient, Hessian_t &hessian){
		Vector_v gradLinear, gradReciprocal, gradTemp;
		int iResp;

		responses = free_terms;
		for (iResp=0;iResp<nResp;iResp++){
			//ConLin
			// resp = free_terms + gradLin*x + gradRecip /x

			// linear part
			gradLinear = sensitivitiesLinear.col(iResp);
			gradient.col(iResp) = gradLinear;
			responses(iResp) += gradLinear.transpose()*primalVar;

			// Reciprocal part
			gradReciprocal = sensitivitiesReciprocal.col(iResp);
			if (!gradReciprocal.isZero((ScalarType)0.0)) {
				const Vector_v primalVarInv = primalVar.cwiseInverse();
				gradTemp = gradReciprocal.cwiseProduct(primalVarInv);
				responses(iResp) += gradTemp.sum();
				gradTemp = gradTemp.cwiseProduct(primalVarInv);
				gradient.col(iResp) -= gradTemp;
				gradTemp = 2*gradTemp.cwiseProduct(primalVarInv);
				hessian.diagonal() +=  gradTemp*dualVar(iResp);
			}

		}
	};
	void Eval(Vector_v primalVar, Vector_r &responses){
		Vector_v gradLinear, gradReciprocal, gradTemp;
		int iResp;

		responses = free_terms;
		for (iResp=0;iResp<nResp;iResp++){
			//ConLin
			// resp = free_terms + gradLin*x + gradRecip /x

			// linear part
			gradLinear = sensitivitiesLinear.col(iResp);
			responses(iResp) += gradLinear.transpose()*primalVar;

			// Reciprocal part
			gradReciprocal = sensitivitiesReciprocal.col(iResp);
			if (!gradReciprocal.isZero((ScalarType)0.0)) {
				const Vector_v primalVarInv = primalVar.cwiseInverse();
				gradTemp = gradReciprocal.cwiseProduct(primalVarInv);
				responses(iResp) += gradTemp.sum();
			}
		}
	}
};

} /* namespace optsection */

#endif /* GLOBALOPTIMISER_APPROXFUNCTION_HPP_ */
