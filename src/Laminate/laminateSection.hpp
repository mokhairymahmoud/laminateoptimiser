/*
 * laminateSection.hpp
 *
 *  Created on: 15 mai 2015
 *      Author: Ramzi
 */

#ifndef LAMINATE_LAMINATESECTION_HPP_
#define LAMINATE_LAMINATESECTION_HPP_

#include <iostream>

#include <Eigen/Dense>
#include "../Section/sectionBase.hpp"
#include "../Laminate/lpfeasible.hpp"
#include "../Miki/Miki.hpp"
#include "../Laminate/varsize.hpp"
#include "../BoundSDP/boundSDP.hpp"

using namespace std;

namespace lampar {
template<SubLaminate_t Type, bool isBalanced, bool isSymmetric = true, int NSUBLAM = 5, typename ScalarType = double>
class laminateSection : public
	optsection::sectionBase<laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>, internal::varsize<isBalanced,isSymmetric>::SIZE+1, ScalarType>
{
public:
	typedef optsection::sectionBase<laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>, internal::varsize<isBalanced,isSymmetric>::SIZE+1, ScalarType> base_t;
	typedef typename base_t::Vector_t Vector_t;
	typedef typename base_t::Hessian_t Hessian_t;
	typedef lpfeasible<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType> base_lpfeasible;
	typedef typename base_lpfeasible::Vector_t Vector_lpfeasible;
	typedef typename base_lpfeasible::Hessian_t Hessian_lpfeasible;
	static const int Size = base_lpfeasible::Size+1;

private:
	base_lpfeasible lpfeasible_;
	boundSDP<Upper, ScalarType> upperBoundThickness_;
	boundSDP<Lower, ScalarType> lowerBoundThickness_;
	static Hessian_t hessian_;
	static Vector_t residual_;
public:
	laminateSection() {	};
	~laminateSection() { };
	void Initialise(ScalarType dual_scale, const Vector_t &var)  {
		lpfeasible_.Initialise(dual_scale);
		upperBoundThickness_.Initialise(dual_scale, var(base_lpfeasible::Size));
		lowerBoundThickness_.Initialise(dual_scale, var(base_lpfeasible::Size));
	}
	void setBoundThickness(ScalarType lowerBoundThickness, ScalarType upperBoundThickness)  {
		upperBoundThickness_.setBound(upperBoundThickness);
		lowerBoundThickness_.setBound(lowerBoundThickness);
	}
	void StepSize(ScalarType &primal_step, ScalarType &dual_step) {
		lpfeasible_.StepSize(primal_step,dual_step);
		upperBoundThickness_.StepSize(primal_step,dual_step);
		lowerBoundThickness_.StepSize(primal_step,dual_step);
	}
	ScalarType DualityGap() {
		ScalarType duality;
		duality = lpfeasible_.DualityGap();
		duality += upperBoundThickness_.DualityGap();
		duality += lowerBoundThickness_.DualityGap();
		return  duality;
	}
	void HessianEval(Hessian_t &hessian){
		Hessian_lpfeasible hessian_lp;

		hessian_lp = hessian.block(0,0,base_lpfeasible::Size,base_lpfeasible::Size);
		lpfeasible_.HessianEval(hessian_lp);
		hessian.block(0,0,base_lpfeasible::Size,base_lpfeasible::Size) = hessian_lp ;
		upperBoundThickness_.Hessian(hessian(base_lpfeasible::Size,base_lpfeasible::Size));
		lowerBoundThickness_.Hessian(hessian(base_lpfeasible::Size,base_lpfeasible::Size));
		hessian_ = hessian;
	}
	void HessianEval(Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> &hessian){
		Hessian_t fixed_hessian = hessian;
		HessianEval(fixed_hessian);
		hessian = fixed_hessian;
	}
	void CalculateResiduals(const Vector_t &var, Vector_t &residual, const ScalarType penalty){
		Vector_lpfeasible residual_lp;
		residual_lp = residual.block(0,0,base_lpfeasible::Size,1);
		lpfeasible_.CalculateResiduals(var.block(0,0,base_lpfeasible::Size,1),residual_lp,penalty);
		residual.block(0,0,base_lpfeasible::Size,1) = residual_lp;
		upperBoundThickness_.CalculateResiduals(residual(base_lpfeasible::Size),penalty);
		lowerBoundThickness_.CalculateResiduals(residual(base_lpfeasible::Size),penalty);
		residual_ = residual;
	}
	void CalculateResiduals(const Eigen::Matrix<ScalarType,Eigen::Dynamic,1> &var,
	                        Eigen::Matrix<ScalarType,Eigen::Dynamic,1> &residual,
	                        const ScalarType penalty){
		Vector_t fixed_var = var;
		Vector_t fixed_residual = residual;
		CalculateResiduals(fixed_var, fixed_residual, penalty);
		residual = fixed_residual;
	}
	void UpdateIncrements(const Vector_t &dvar){
		lpfeasible_.UpdateIncrements(dvar.block(0,0,base_lpfeasible::Size,1));
		upperBoundThickness_.UpdateIncrements(dvar(base_lpfeasible::Size));
		lowerBoundThickness_.UpdateIncrements(dvar(base_lpfeasible::Size));
	}
	void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step){
		lpfeasible_.UpdateVariables(primal_step,dual_step);
		upperBoundThickness_.UpdateVariables(primal_step,dual_step);
		lowerBoundThickness_.UpdateVariables(primal_step,dual_step);
	}
	Hessian_t getHessian(){
		return hessian_;
	}
	Vector_t getResidual(){
		return residual_;
	}
};
template<SubLaminate_t Type, bool isBalanced, bool isSymmetric, int NSUBLAM, typename ScalarType>
typename lampar::laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>::Hessian_t lampar::laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>::hessian_;
template<SubLaminate_t Type, bool isBalanced, bool isSymmetric, int NSUBLAM, typename ScalarType>
typename lampar::laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>::Vector_t lampar::laminateSection<Type, isBalanced, isSymmetric, NSUBLAM, ScalarType>::residual_;

}  /* namespace lampar */


#endif /* LAMINATE_LAMINATESECTION_HPP_ */
