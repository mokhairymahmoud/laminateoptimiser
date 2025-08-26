/*
 * solutionProblem1.hpp
 *
 *  Created on: 24 oct. 2016
 *      Author: Ramzi
 */

#ifndef TESTPROBLEMS_SOLUTIONPROBLEM1_HPP_
#define TESTPROBLEMS_SOLUTIONPROBLEM1_HPP_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>

template <int nResp, int nVar, typename ScalarType = double>
class solutionProblem1 {
private:
	int NVAR;
	int NRESP;
public:
	typedef Eigen::Matrix<ScalarType,nResp,1> Vector_r;
	typedef Eigen::Matrix<ScalarType,nVar,nResp> Matrix_t;
	typedef Eigen::Matrix<ScalarType,nVar,nVar> Hessian_t;

private:
	Vector_r resp;
	Matrix_t gradient;
	Hessian_t hessian;
public:
	solutionProblem1() {
		NVAR = nVar;
		NRESP = nResp;
	}
	virtual ~solutionProblem1(){

	}
	int getNumbResponse() {
		return this->NRESP;
	}
	int getNumbVariable() {
		return this->NVAR;
	}
	void setResponse(Vector_r resp) {
		this->resp = resp;
	}
	Vector_r getResponse() {
		return this->resp;
	}
	void setGradient(Matrix_t gradient) {
		this->gradient = gradient;
	}
	Matrix_t getGradient() {
		return this->gradient;
	}
	void setHessian(Hessian_t hessian) {
		this->hessian = hessian;
	}
	Hessian_t getHessian() {
		return this->hessian;
	}
};

#endif /* TESTPROBLEMS_SOLUTIONPROBLEM1_HPP_ */
