/*
 * responseInterface.hpp
 *
 *  Created on: 18 oct. 2016
 *      Author: Ramzi
 */

#ifndef INTERFACEOPTIMISER_RESPONSEINTERFACE_HPP_
#define INTERFACEOPTIMISER_RESPONSEINTERFACE_HPP_
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>
#include "../GlobalOptimiser/approxFunction.hpp"

using namespace std;
using namespace Eigen;

namespace globopt {


template<typename designVariable, typename responseProblem, typename approximationProblem, typename analysisSolver, typename ScalarType = double>
class responseInterface
{
private:
	approximationProblem *pApproximationProblem;
	analysisSolver *pAnalysisSolver;
	designVariable designVariables;
	responseProblem responses;


public:

	responseInterface(designVariable designVariables, approximationProblem *pApproximationProblem, analysisSolver *pAnalysisSolver) {
		this->designVariables = designVariables;
		this->pApproximationProblem = pApproximationProblem;
		this->pAnalysisSolver = pAnalysisSolver;
	}
	~responseInterface() {
	}
	void setDesignPoint(designVariable designVariables) {
		this->designVariables = designVariables;
	}
	designVariable getDesignPoint() {
		return this->designVariables;
	}
	void setApproximationProblem(approximationProblem *pApproximationProblem) {
		this->pApproximationProblem = pApproximationProblem;
	}
	approximationProblem * getApproximationProblem() {
		return this->pApproximationProblem;
	}
	responseProblem getResponseProblem() {
		return this->responses;
	}

	responseProblem femSolver(designVariable designVariables){ // will check later about the output
		return this->responses = this->pAnalysisSolver->solve(designVariables);
		//	this->designVariables = designVariables;
	}

	responseProblem evalApproximationProblem(designVariable designVariables, approximationProblem *pApproximationProblem) {
		this->responses = this->pApproximationProblem->eval(designVariables);
		pApproximationProblem = this->pApproximationProblem;
		//	this->designVariables = designVariables;
		return this->responses;
	}

	responseProblem evalApproximationProblem(designVariable designVariables) {
		this->responses = this->pApproximationProblem->eval(designVariables);
		//	this->designVariables = designVariables;
		return this->responses;
	}

	approximationProblem * buildApproximationProblem(designVariable designVariables) {
		//	this->designVariables = designVariables;
		return this->pApproximationProblem->build(designVariables);
	}

};

}  /* namespace globopt */



#endif /* INTERFACEOPTIMISER_RESPONSEINTERFACE_HPP_ */
