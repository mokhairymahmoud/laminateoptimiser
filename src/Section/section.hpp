/*
 * section.hpp
 *
 *  Created on: 15 mai 2015
 *      Author: Ramzi
 */

#ifndef SECTION_SECTION_HPP_
#define SECTION_SECTION_HPP_
#include <Eigen/Dense>
#include <cstddef>
#include <iostream>

namespace optsection {
	template<typename ScalarType = double>
	class section {
	private:
		virtual ScalarType DualityGap_impl() = 0;
		virtual void StepSize_impl(ScalarType &primal_step, ScalarType &dual_step) =0;
		virtual void HessianEval_impl(void *hessian) =0;
		virtual void CalculateResiduals_impl(const void *var, void *residual, const ScalarType penalty) =0;
		virtual void UpdateIncrements_impl(const void *dvar) =0;
		virtual void UpdateVariables_impl(const ScalarType primal_step,const ScalarType dual_step) =0;
		virtual void setHessian_impl(void *hessian) =0;
		virtual void Initialise_impl(ScalarType dual_scale, const void *var) = 0;
		virtual int getSize_impl() = 0;
		virtual Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> getHessian_impl() =0;
		virtual Eigen::Matrix<ScalarType,Eigen::Dynamic,1> getResidual_impl() =0;

	public:
		section(){};
		virtual ~section(){};
		int getSize(){
			return this->getSize_impl();
		}
		template<int nVar>
		void Initialise(ScalarType dual_scale, const Eigen::Matrix<ScalarType,nVar,1> &var)  {
			this->Initialise_impl(dual_scale,&var);
		}

		template<int nVar>
		void HessianEval(Eigen::Matrix<ScalarType,nVar,nVar> &hessian){
			this->HessianEval_impl(&hessian);
		}
		template<int nVar>
		void CalculateResiduals(const Eigen::Matrix<ScalarType,nVar,1> &var, Eigen::Matrix<ScalarType,nVar,1> &residual, const ScalarType penalty){
			this->CalculateResiduals_impl(&var, &residual, penalty);
		}
		template<int nVar>
		void UpdateIncrements(const Eigen::Matrix<ScalarType,nVar,1> &dvar){
			this->UpdateIncrements_impl(&dvar);
		}

		void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step){
			this->UpdateVariables_impl(primal_step, dual_step);
		}

		void StepSize(ScalarType &primal_step, ScalarType &dual_step){
			this->StepSize_impl(primal_step, dual_step);
		}

		ScalarType DualityGap() {
			return this->DualityGap_impl();
		}

		template<int nVar>
		void setHessian(Eigen::Matrix<ScalarType,nVar,nVar> hessian){
			this->setHessian_impl(&hessian);
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> getHessian(){
			return this->getHessian_impl();
		}
		Eigen::Matrix<ScalarType,Eigen::Dynamic,1> getResidual(){
			return this->getResidual_impl();
		}
	};


} /* namespace optsection */

#endif /* SECTION_SECTION_HPP_ */
