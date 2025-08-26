/*
 * SDPA.hpp
 *
 *  Created on: 6 Jan 2015
 *      Author: root
 */

#ifndef SDPA_HPP_
#define SDPA_HPP_

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "../SDPA/parameter.hpp"

namespace SDPA {

template<typename SDCone, int Dim, int Size, typename ScalarType = double>
class Constraint  {
public:
	typedef Eigen::Matrix<ScalarType,Dim,Dim> Matrix_t;
	typedef Eigen::Matrix<ScalarType,Size,1> Vector_t;
	typedef Eigen::Matrix<ScalarType,Size,Size> Hessian_t;
private:
	typedef Eigen::LLT<Matrix_t> Cholesky_t;
	Matrix_t primal_, dprimal_, dual_, ddual_, inv_primal_, duality_residual_;
	SDCone * SDPA_Cone()
	{
		return static_cast<SDCone *>(this);
	}
	ScalarType FindStep(const Matrix_t &val, const Matrix_t &dval) const
	{
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> es(dval,val);
		return es.eigenvalues().minCoeff();
	}
public:
	void Initialise(ScalarType dual_scale, const Vector_t &var)  {
		SDPA_Cone()->Initval(var,primal_);
		Cholesky_t L = primal_.llt();
		inv_primal_ = Matrix_t::Identity();
		L.solveInPlace(inv_primal_);
		dual_ = dual_scale*inv_primal_; // dual_ = dual_scale * inv_primal
		dprimal_ = Matrix_t::Zero();
		ddual_ = Matrix_t::Zero();
	}
	void Initialise(ScalarType dual_scale)  {
		this->Initialise(dual_scale, Vector_t::Zero());
	}
	void Eval(const Vector_t &var, Matrix_t &primal) {
		SDPA_Cone()->Eval(var,primal);
	}
	void AdEval(const Matrix_t &dual, Vector_t &residual)  {
		SDPA_Cone()->AdEval(dual,residual);
	}
	void HessianEval(const Matrix_t &primal, const Matrix_t &dual, Hessian_t &hessian)  {
		SDPA_Cone()->HessianEval(primal,dual,hessian);
	}
	Constraint() {};
	Constraint(const Vector_t & var0, ScalarType dual_scale)
	{
		this->Initialise(dual_scale,var0);
	}
	~Constraint(){};
	void CalculateResiduals(Vector_t &residual, const ScalarType penalty)
	{
		AdEval(dual_,residual);
		duality_residual_.noalias() = penalty*Matrix_t::Identity()-primal_*dual_-dprimal_*ddual_;
		ddual_ = inv_primal_*duality_residual_;
		ddual_.noalias() = (ddual_ + ddual_.transpose())/2.0;
		AdEval(ddual_,residual);
	}
	void Hessian(Hessian_t &hessian)
	{
		HessianEval(inv_primal_,dual_,hessian);
	}
	ScalarType DualityGap() const
	{
		Matrix_t temp;
		temp.noalias() = (primal_+dprimal_)*(dual_+ddual_);
		return temp.trace();
	}
	void UpdateIncrements(const Vector_t &dvar)
	{
		Eval(dvar,dprimal_) ;
		duality_residual_ -= dprimal_*dual_;
		ddual_ = inv_primal_*duality_residual_;
		// symmetrisation (which is damn ugly)
		duality_residual_ = ddual_.transpose(); // duality_residual used as temp
		ddual_ = (ddual_ +duality_residual_)/2.0;
	}
	void StepSize(ScalarType &primal_step, ScalarType &dual_step) const
	{
		ScalarType eval, step;
		eval = FindStep(primal_,dprimal_);
		step = Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
		eval = FindStep(dual_,ddual_);
		step = Parameter<ScalarType>::Step_Size_Control(eval);
		dual_step = dual_step < step? dual_step: step;
	}
	void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step)
	{
		primal_ += primal_step*dprimal_;
		dual_ += dual_step*ddual_;
		Cholesky_t L = primal_.llt();
		inv_primal_ = Matrix_t::Identity();
		L.solveInPlace(inv_primal_);
		dprimal_ = Matrix_t::Zero();
		ddual_ = Matrix_t::Zero();
	}
};

/*
 *
template<typename SDCone, int Dim, int Size, typename ScalarType = double>
class Constraint {
private:
	typedef Eigen::Matrix<ScalarType,Dim,Dim> SDPA_MatrixType;
	typedef Eigen::Matrix<ScalarType,Size,1> SDPA_VectorType;
	typedef Eigen::LLT<SDPA_MatrixType> SDPA_Cholesky;
	SDPA_MatrixType primal_, dprimal_, dual_, ddual_, inv_primal_, duality_residual_;
	SDConeBase<SDCone,Dim,ScalarType> * SDPA_Cone_;
	ScalarType FindStep(const SDPA_MatrixType &val, const SDPA_MatrixType &dval)
	{
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d> es(dval,val);
		return es.eigenvalues().minCoeff();
	}
public:
	Constraint(SDConeBase<SDCone,Dim,ScalarType> * cone,const SDPA_VectorType & var0,ScalarType dual_scale)  {
		SDPA_Cone_ = cone;
		SDPA_Cone_->Initialise(primal_);
		SDPA_Cone_->Eval(var0, primal_);
		SDPA_Cholesky L = primal_.llt();
		inv_primal_ = SDPA_MatrixType::Identity();
		L.solveInPlace(inv_primal_);
		dual_ = dual_scale*inv_primal_; // dual_ = dual_scale * inv_primal
		dprimal_ = SDPA_MatrixType::Zero();
		ddual_ = SDPA_MatrixType::Zero();
	}
	~Constraint(){};
	void CalculateResiduals(SDPA_VectorType &residual, const ScalarType penalty)
	{
		SDPA_Cone_->AdEval(dual_,residual);
		duality_residual_.noalias() = penalty*SDPA_MatrixType::Identity()-primal_*dual_-dprimal_*ddual_;
		ddual_ = inv_primal_*duality_residual_;
		ddual_.noalias() = (ddual_ + ddual_.transpose())/2.0;
		SDPA_Cone_->AdEval(dual_,residual);
	}
	ScalarType DualityGap() const
	{
		SDPA_MatrixType temp;
		temp.noalias() = (primal_+dprimal_)*(dual_+ddual_);
		return temp.trace();
	}
	void UpdateIncrements(const SDPA_VectorType &dvar)
	{
		SDPA_Cone_->Eval(dvar,dprimal_) ;
		duality_residual_ -= dprimal_*dual_;
		ddual_ = inv_primal_*duality_residual_;
		ddual_.noalias() = (ddual_ + ddual_.transpose())/2.0;
	}
	void StepSize(ScalarType &primal_step, ScalarType &dual_step) const
	{
		ScalarType eval, step;
		eval = FindStep(primal_,dprimal_);
		step = Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
		eval = FindStep(dual_,ddual_);
		step = Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
	}
	void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step)
	{
		primal_ += primal_step*dprimal_;
		dual_ += primal_step*ddual_;
		SDPA_Cholesky L = primal_.llt();
		inv_primal_ = SDPA_MatrixType::Identity();
		L.solveInPlace(inv_primal_);
		dprimal_ = SDPA_MatrixType::Zero();
		ddual_ = SDPA_MatrixType::Zero();
	}
};
*/

} // end namespace SDPA

#endif /* SDPA_HPP_ */

