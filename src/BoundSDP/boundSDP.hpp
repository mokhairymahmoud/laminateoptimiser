/*
 * boundSDP.hpp
 *
 *  Created on: 2 mai 2015
 *      Author: Ramzi
 */

#ifndef BOUNDSDP_BOUNDSDP_HPP_
#define BOUNDSDP_BOUNDSDP_HPP_

#include "../SDPA/scalarSDP.hpp"

namespace lampar {

enum Direction_t{Upper, Lower};

template<Direction_t Type, typename ScalarType = double>
class boundSDP {};

template<typename ScalarType>
class boundSDP<Upper, ScalarType> : public
      scalarSDP::Constraint<boundSDP<Upper, ScalarType>, ScalarType>
{
private:
	ScalarType bound_;
public:
	typedef typename scalarSDP::Constraint<boundSDP<Upper, ScalarType>, ScalarType>  Base_t;

public:
	boundSDP() {};
	boundSDP(ScalarType bound) {
		bound_ = bound;
	};
	void setBound(ScalarType bound) {
		bound_ = bound;
	}
//	boundSDP(ScalarType dualscale) : Base_t(0,dualscale) {}
	void Initval(const ScalarType &var, ScalarType &slack) {
		Eval(var,slack);
		slack += bound_;
	}
	void Eval(const ScalarType &var, ScalarType &slack) const {
		slack = -var;
	}
	void AdEval(const ScalarType &dual,  ScalarType &residual) const {
		residual -= dual;
	}
	void HessianEval(const ScalarType &slack, const ScalarType &dual, ScalarType &hessian) {
		hessian += slack*dual;
	}

};

template<typename ScalarType>
class boundSDP<Lower, ScalarType> : public
      scalarSDP::Constraint<boundSDP<Lower, ScalarType>, ScalarType>
{
private:
	ScalarType bound_;
public:
	typedef typename scalarSDP::Constraint<boundSDP<Lower, ScalarType>, ScalarType>  Base_t;

public:
	boundSDP() {};
	boundSDP(ScalarType bound) {
		bound_ = bound;
	};
	void setBound(ScalarType bound) {
		bound_ = bound;
	}
//	boundSDP(ScalarType dualscale) : Base_t(0,dualscale) {}
	~boundSDP(){};
	void Initval(const ScalarType &var, ScalarType &slack) {
		Eval(var,slack);
		slack -= bound_;
	}
	void Eval(const ScalarType &var, ScalarType &slack) const {
		slack = var;
	}
	void AdEval(const ScalarType &dual,  ScalarType &residual) const {
		residual += dual;
	}
	void HessianEval(const ScalarType &slack, const ScalarType &dual, ScalarType &hessian) {
		hessian += slack*dual;
	}
};
}  /* namespace lampar */

#endif /* BOUNDSDP_BOUNDSDP_HPP_ */
