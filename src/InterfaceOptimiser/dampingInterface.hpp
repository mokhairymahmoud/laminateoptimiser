/*
 * dampingInterface.hpp
 *
 *  Created on: 19 oct. 2016
 *      Author: Ramzi
 */

#ifndef INTERFACEOPTIMISER_DAMPINGINTERFACE_HPP_
#define INTERFACEOPTIMISER_DAMPINGINTERFACE_HPP_

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>
#include "../GlobalOptimiser/approxFunction.hpp"

using namespace std;
using namespace Eigen;

namespace globopt {


template<typename designVariable, typename responseProblem, typename dampingProblem, typename ScalarType = double>
class dampingInterface
{
public:
	typedef typename dampingProblem::Vector_r Vector_r;

private:

	dampingProblem *pDampingProblem;
	designVariable designVariables;
	Vector_r dampingFactor;   // should change  its place
	responseProblem dampingResponses;


private:
	ScalarType dampingRatio(ScalarType gapTerm, ScalarType dampingLow) {
		if (gapTerm >= 0)
			return gapTerm + (ScalarType)1.0;
		else
			return (ScalarType)1.0 + ((ScalarType)1.0 - dampingLow)*tanh(gapTerm/((ScalarType)1.0-dampingLow));
	}

public:

	dampingInterface(dampingProblem *pDampingProblem, designVariable designVariables) {
		this->designVariables = designVariables;
		this->pDampingProblem = pDampingProblem;
	}
	~dampingInterface() {
	}
	void setDesignPoint(designVariable designVariables) {
		this->designVariables = designVariables;
	}
	designVariable getDesignPoint() {
		return this->designVariables;
	}
	void setDampingInterface(dampingProblem *pDampingProblem) {
		this->pDampingProblem = pDampingProblem;
	}
	dampingProblem * getDampingInterface() {
		return this->pDampingProblem;
	}
	void setDampingFactor(Vector_r dampingFactor) {
		this->dampingFactor = dampingFactor;
	}
	Vector_r getDampingFactor() {
		return this->dampingFactor;
	}
	responseProblem getResponseDamping() {
		return this->dampingResponses;
	}


	responseProblem evalDampingProblem(designVariable designVariables, designVariable candidateDesignVariables, dampingProblem *pDampingProblem) {
		this->pDampingProblem->eval(designVariables, candidateDesignVariables, this->dampingResponses);
	//	this->designVariables = designVariables;
		pDampingProblem = this->pDampingProblem;
		return this->dampingResponses;
	}
	responseProblem evalDampingProblem(designVariable designVariables, designVariable candidateDesignVariables) {
		this->pDampingProblem->eval(designVariables, candidateDesignVariables, this->dampingResponses);
	//	this->designVariables = designVariables;
		return this->dampingResponses;
	}

	dampingProblem * buildDampingProblem(designVariable designVariables, designVariable initialDesignVariables /*sensitivities can be used too*/) {
	//	this->designVariables = designVariables;
		return this->pDampingProblem->build(designVariables, initialDesignVariables);
	}
	void updateDampingFactor(designVariable designVariables, designVariable candidateDesignVariables, responseProblem responses, responseProblem candidateResponses) {
		Vector_r gapTerm, damping;
		int nResp = responses.getNumResponses();
		this->pDampingProblem->eval(designVariables, candidateDesignVariables, this->dampingResponses);
		//	this->designVariables = designVariables;

		gapTerm = (responses.getResponse() - candidateResponses.getResponse()).cwiseInverse(this->dampingResponses.getResponse().cwiseProduct(this->dampingFactor));
		for (int iResp=0;iResp<nResp;iResp++){
			dampingFactor(iResp) = dampingRatio(gapTerm(iResp), 0.1 /*will be changed with option */);
		}
		this->dampingResponses.getResponse() = this->dampingFactor.cwiseProduct(this->dampingResponses.getResponse());
	}

};

}  /* namespace globopt */


#endif /* INTERFACEOPTIMISER_DAMPINGINTERFACE_HPP_ */
