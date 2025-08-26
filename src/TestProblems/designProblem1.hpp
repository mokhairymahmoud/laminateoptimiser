/*
 * solutionDesign.hpp
 *
 *  Created on: 24 oct. 2016
 *      Author: Ramzi
 */

#ifndef TESTPROBLEMS_DESIGNPROBLEM1_HPP_
#define TESTPROBLEMS_DESIGNPROBLEM1_HPP_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>

template <int nVar, typename ScalarType = double>
class DesignProblem1 {
public:
	int NVAR;
	typedef Eigen::Matrix<ScalarType,nVar,1> Vector_v;


private:
	Vector_v designVariables;

public:
	DesignProblem1() {
		NVAR = nVar;
	}
	virtual ~DesignProblem1(){

	}
	int getNumbVariable() {
		return this->NVAR;
	}
	void setDesignVariable(Vector_v designVariables) {
		this->designVariables = designVariables;
	}
	Vector_v getDesignVariable() {
		return this->designVariables;
	}
};

#endif /* TESTPROBLEMS_DESIGNPROBLEM1_HPP_ */
