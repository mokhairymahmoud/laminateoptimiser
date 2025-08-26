/*
 * scalarSDP.hpp
 *
 *  Created on: 2 mai 2015
 *      Author: Ramzi
 */

#ifndef SCALARSDP_HPP_
#define SCALARSDP_HPP_

#include <iostream>
#include "../SDPA/parameter.hpp"

namespace scalarSDP {

template<typename SDCone, typename ScalarType = double>
class Constraint  {

private:

	ScalarType slack_, dslack_, dual_, ddual_, inv_slack_, duality_residual_;
	SDCone * scalarSDP_Cone() {
		return static_cast<SDCone *>(this);
	}
	ScalarType FindStep(const ScalarType &val, const ScalarType &dval) const {
		return dval/val;
	}
public:
	void Initialise(ScalarType dual_scale, const ScalarType &var)  {
		scalarSDP_Cone()->Initval(var,slack_);
		inv_slack_ = 1.0/slack_;
		dual_ = dual_scale*inv_slack_;
		dslack_ = 0;
		ddual_ = 0;
	}
	void Initialise(ScalarType dual_scale)  {
		this->Initialise(dual_scale, 0);
	}
	void Eval(const ScalarType &var, ScalarType &slack) {
		scalarSDP_Cone()->Eval(var,slack);
	}
	void AdEval(const ScalarType &dual, ScalarType &residual)  {
		scalarSDP_Cone()->AdEval(dual,residual);
	}
	void HessianEval(const ScalarType &slack, const ScalarType &dual, ScalarType &hessian)  {
		scalarSDP_Cone()->HessianEval(slack,dual,hessian);
	}
	Constraint() {};
	Constraint(const ScalarType & var0, ScalarType dual_scale) {
		this->Initialise(dual_scale,var0);
	}
	~Constraint(){};
	void CalculateResiduals(ScalarType &residual, const ScalarType penalty) {
		AdEval(dual_,residual);
		duality_residual_ = penalty-slack_*dual_-dslack_*ddual_;
		ddual_ = inv_slack_*duality_residual_;
		AdEval(ddual_,residual);
	}
	void Hessian(ScalarType &hessian) {
		HessianEval(inv_slack_,dual_,hessian);
	}
	ScalarType DualityGap() const {
		ScalarType temp;
		temp = (slack_+dslack_)*(dual_+ddual_);
		return temp;
	}
	void UpdateIncrements(const ScalarType &dvar)
	{
		Eval(dvar,dslack_) ;
		duality_residual_ -= dslack_*dual_;
		ddual_ = inv_slack_*duality_residual_;
	}
	void StepSize(ScalarType &slack_step, ScalarType &dual_step) const
	{
		ScalarType eval, step;
		eval = FindStep(slack_,dslack_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		slack_step = slack_step < step? slack_step: step;
		eval = FindStep(dual_,ddual_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		dual_step = dual_step < step? dual_step: step;
	}
	void UpdateVariables(const ScalarType slack_step,const ScalarType dual_step)
	{
		slack_ += slack_step*dslack_;
		dual_ += dual_step*ddual_;
		inv_slack_ = 1.0/slack_;
		dslack_ = 0;
		ddual_ = 0;
	}
};
} // end namespace scalarSDP
#endif /* SCALARSDP_HPP_ */
