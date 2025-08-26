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
		dimObjSet = 0;
		NVAR = nVar;
		NRESP = nResp;
	};
	virtual ~approxFunction() {};
	int getBooleanVector(Vector_r &booleanVect) {
		booleanVect = booleanVector;
		return dimObjSet;
	}
	void Eval(Vector_v primalVar, Vector_r dualVar, Vector_r &responses, Matrix_t &gradient, Hessian_t &hessian){
		Vector_v gradLinear, gradReciprocal, gradTemp, primalVarInv;
		int iResp;

		responses = free_terms;
		primalVarInv = primalVar.cwiseInverse();
		for (iResp=0;iResp<nResp;iResp++){
			//ConLin
			// resp = free_terms + gradLin*x + gradRecip /x

			// linear part
			gradLinear = sensitivitiesLinear.col(iResp);
			gradient.col(iResp) = gradLinear;
			responses(iResp) += gradLinear.transpose()*primalVar;

			// Reciprocal part
			gradReciprocal = sensitivitiesReciprocal.col(iResp);
			gradTemp = gradReciprocal.cwiseProduct(primalVarInv);
			responses(iResp) += gradTemp.sum();
			gradTemp = gradTemp.cwiseProduct(primalVarInv);
			gradient.col(iResp) -= gradTemp;
			gradTemp = 2*gradTemp.cwiseProduct(primalVarInv);
			hessian.diagonal() +=  gradTemp*dualVar(iResp);

		}
	};
	void Eval(Vector_v primalVar, Vector_r &responses){
		Vector_v gradLinear, gradReciprocal, gradTemp, primalVarInv;
		int iResp;

		responses = free_terms;
		primalVarInv = primalVar.cwiseInverse();
		for (iResp=0;iResp<nResp;iResp++){
			//ConLin
			// resp = free_terms + gradLin*x + gradRecip /x

			// linear part
			gradLinear = sensitivitiesLinear.col(iResp);
			responses(iResp) += gradLinear.transpose()*primalVar;

			// Reciprocal part
			gradReciprocal = sensitivitiesReciprocal.col(iResp);
			gradTemp = gradReciprocal.cwiseProduct(primalVarInv);
			responses(iResp) += gradTemp.sum();
		}
	}
	void InitialiseProblem1(Vector_v primalVar){
		ScalarType funct, funct0;
		Vector_v grad;
		ScalarType x1, x2;
		ScalarType freeterm;
		Vector_v gradLinear;
		Vector_v gradReciprocal;

		// a 3 truss bars problem under stress constraints
		/*
		 * min V = 100*(2*sqrt(2)*x1+x2)
		 * s.t.
		 *    g1 : (x2+sqrt(2)*x1)/(2*x1*x2+sqrt(2)*x1*x1)-1 <= 0
		 *    g2 : 1.0/(x1+sqrt(2)*x2)-1 <= 0
		 *    g3 : x2/(2*x1*x2+sqrt(2)*x1*x1)-3.0/4.0 <= 0
		 *    g4 : -x1 <= 0
		 *    g5 : -x2 <= 0
		 */

		booleanVector<<1,0,0,0;dimObjSet = 1;
		x1 = primalVar(0) ; x2 = primalVar(1);
		// objective
		funct = 100 * (2*sqrt(2)*x1+x2);funct0=funct;
		funct = funct/funct0;
		grad(0) = 200 * sqrt(2)/funct0;
		grad(1) = 100/funct0;
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(0) = freeterm;
		sensitivitiesLinear.col(0) = gradLinear;
		sensitivitiesReciprocal.col(0) = gradReciprocal;
		// constraint 1
		funct = (x2+sqrt(2)*x1)/(2*x1*x2+sqrt(2)*x1*x1)-1;
		grad(0) = -(sqrt(2)*x1*x2+x1*x1+x2*x2)/pow(sqrt(2)*x1*x2+x1*x1,2);
		grad(1) = -1.0/(sqrt(2)*pow(sqrt(2)*x2+x1,2));
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(1) = freeterm;
		sensitivitiesLinear.col(1) = gradLinear;
		sensitivitiesReciprocal.col(1) = gradReciprocal;
		// constraint 2
		funct = 1.0/(x1+sqrt(2)*x2)-1;
		grad(0) = -1.0/(pow(sqrt(2)*x2+x1,2));
		grad(1) = -sqrt(2)/(pow(sqrt(2)*x2+x1,2));
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(2) = freeterm;
		sensitivitiesLinear.col(2) = gradLinear;
		sensitivitiesReciprocal.col(2) = gradReciprocal;
		// constraint 3
		funct = 4.0/3.0*x2/(2*x1*x2+sqrt(2)*x1*x1)-1;
		grad(0) = -4.0/3.0*x2*(x2+sqrt(2)*x1)/pow(sqrt(2)*x1*x2+x1*x1,2);
		grad(1) = 4.0/3.0*1.0/(sqrt(2)*pow(sqrt(2)*x2+x1,2));
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(3) = freeterm;
		sensitivitiesLinear.col(3) = gradLinear;
		sensitivitiesReciprocal.col(3) = gradReciprocal;
	}
	void InitialiseProblem2(Vector_v primalVar){
		ScalarType funct;
		Vector_v grad;
		ScalarType x1, x2;
		ScalarType freeterm;
		Vector_v gradLinear;
		Vector_v gradReciprocal;

		// a linear problem
		/*
		 * min f = .2*(2*x1+x2)
		 * s.t.
		 *    g1 : -x1-2*x2+7.5 <= 0
		 *    g2 : x1 >= 2.5
		 *    g3 : x2 >= 1.25
		 */

		booleanVector<<1,0;dimObjSet = 1;
		x1 = primalVar(0) ; x2 = primalVar(1);
		// objective
		funct = .2 * (2*x1+x2);
		grad(0) = .4/funct;
		grad(1) = .2/funct;funct = 1.0;
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(0) = freeterm;
		sensitivitiesLinear.col(0) = gradLinear;
		sensitivitiesReciprocal.col(0) = gradReciprocal;
		// constraint 1
		funct = -x1/7.5-2*x2/7.5+1;
		grad(0) = -1.0/7.5;
		grad(1) = -2.0/7.5;
		LinearOrReciprocal(primalVar, funct, grad, freeterm, gradLinear, gradReciprocal);
		free_terms(1) = freeterm;
		sensitivitiesLinear.col(1) = gradLinear;
		sensitivitiesReciprocal.col(1) = gradReciprocal;

	}
};

} /* namespace optsection */

#endif /* GLOBALOPTIMISER_APPROXFUNCTION_HPP_ */
