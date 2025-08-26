/*
 * sectionHelper.hpp
 *
 *  Created on: 15 mai 2015
 *      Author: Ramzi
 */

#ifndef SECTION_SECTIONBASE_HPP_
#define SECTION_SECTIONBASE_HPP_

#include "../Section/section.hpp"
#include <Eigen/Dense>

namespace optsection {
template<class derivedSection, int nVar, typename ScalarType = double>
class sectionBase  :public section<ScalarType>{
private:
	static const int size = nVar;
	typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,1> Vector_dynamic;
	typedef Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> Hessian_dynamic;
public:
	typedef Eigen::Matrix<ScalarType,nVar,1> Vector_t;
	typedef Eigen::Matrix<ScalarType,nVar,nVar> Hessian_t;
	sectionBase() {

	};
	~sectionBase() {};
private:
	int getSize_impl(){
		return size;
	}
	void Initialise_impl(ScalarType dual_scale, const void *var){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		const Vector_dynamic * pVar = static_cast<const Vector_dynamic *> (var);
		pSection->Initialise(dual_scale,*pVar);
	}

	void StepSize_impl(ScalarType &primal_step, ScalarType &dual_step){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		pSection->StepSize(primal_step,  dual_step);
	}
	ScalarType DualityGap_impl(){
		ScalarType  pprimal_step;
		derivedSection * pSection = static_cast<derivedSection *> (this);
		pprimal_step = pSection->DualityGap();
		return pprimal_step;
	}
	void HessianEval_impl(void *hessian){
		derivedSection * pSection = dynamic_cast<derivedSection *> (this);
		Hessian_dynamic * pHessian = static_cast<Hessian_dynamic *> (hessian);
		pSection->HessianEval(*pHessian);
	}
	void CalculateResiduals_impl(const void *var, void *residual, const ScalarType penalty){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		const Vector_dynamic * pVar = static_cast<const Vector_dynamic *> (var);
		Vector_dynamic * pResidual = static_cast<Vector_dynamic *> (residual);
		pSection->CalculateResiduals(*pVar, *pResidual,penalty);
	}
	void UpdateIncrements_impl(const void *dvar){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		const Vector_dynamic * pDvar = static_cast<const Vector_dynamic *> (dvar);
		pSection->UpdateIncrements(*pDvar);
	}
	void UpdateVariables_impl(const ScalarType primal_step,const ScalarType dual_step){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		pSection->UpdateVariables(primal_step,  dual_step);
	}
	void setHessian_impl(void *hessian){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		Hessian_dynamic * pHessian = static_cast<Hessian_dynamic *> (hessian);
		pSection->setHessian(*pHessian);
	}
	Hessian_dynamic getHessian_impl(){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		Hessian_dynamic Hessian;
		Hessian = pSection->getHessian();
		return Hessian;
	}
	Vector_dynamic getResidual_impl(){
		derivedSection * pSection = static_cast<derivedSection *> (this);
		Hessian_dynamic Vector;
		Vector = pSection->getResidual();
		return Vector;
	}
};


} // end namespace optsection



#endif /* SECTION_SECTIONBASE_HPP_ */
