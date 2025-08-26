/*
 * lpfSDPCone.h
 *
 *  Created on: 2 Feb 2015
 *      Author: root
 */

#ifndef LPFSDPCONE_H_
#define LPFSDPCONE_H_
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <vector>
#include "../fSDPCone/fSDPCone.hpp"
#include "../Laminate/varsize.hpp"
#include "../Laminate/cltweight.hpp"
#include "../Miki/Miki.hpp"

namespace lampar {
template<SubLaminate_t Type, bool isBalanced, bool isSymmetric = true, int NSUBLAM = 5, typename ScalarType = double>
class lpfSDPCone {
public:
	static const int NLAMPAR = internal::varsize<isBalanced,isSymmetric>::NVAR;
    static const int NLAMPARTYPE = internal::varsize<isBalanced,isSymmetric>::NVARTYPE;
    static const int Size = NLAMPAR*NLAMPARTYPE;

    typedef typename fSDPCone<3, NLAMPAR, ScalarType>::Matrix_t MikiMatrix_t;
	typedef typename fSDPCone<3, NLAMPAR, ScalarType>::Vector_t MikiVector_t;
	typedef typename fSDPCone<3, NLAMPAR, ScalarType>::Hessian_t MikiHessian_t;
	typedef Eigen::LLT<MikiHessian_t> MikiCholesky_t;

	typedef typename Eigen::Matrix<ScalarType,Size,1> Vector_t;
	typedef typename Eigen::Matrix<ScalarType,Size,Size> Hessian_t;
	typedef Eigen::LLT<Hessian_t> Cholesky_t;

private:
	static lampar::internal::cltweight<isSymmetric,NSUBLAM,ScalarType> clt;

	fSDPCone<3, NLAMPAR, ScalarType> sublaminate_[NSUBLAM];
	MikiCholesky_t sublaminate_cholesky_hessian_[NSUBLAM];
	MikiVector_t sublaminate_var_[NSUBLAM], sublaminate_dvar_[NSUBLAM];
	Vector_t duality_residual_;
	Cholesky_t cholesky_dual_;
public:
	lpfSDPCone(void) {};
	~lpfSDPCone(void) {};
	void Initialise(const MikiMatrix_t &convSetM0, const vector<MikiMatrix_t> &convSetM, ScalarType dual_scale)  {
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			sublaminate_var_[iSubLam] = MikiVector_t::Zero();
			sublaminate_[iSubLam].InitConvexSet(convSetM0,convSetM);
			sublaminate_[iSubLam].Initialise(dual_scale,sublaminate_var_[iSubLam]);
		}
	}

	void Eval(Vector_t &var)  {
		int iBlock;
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++)  {
			for (int ilampar_type = 0; ilampar_type < NLAMPARTYPE; ilampar_type ++)  {
				iBlock = ilampar_type*NLAMPAR;
				var.block(iBlock,0,NLAMPAR,1) += clt.weight[ilampar_type][iSubLam]*sublaminate_var_[iSubLam];
			}
		}
	}

	ScalarType DualityGap() {
		ScalarType gap = 0.0;
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			gap += sublaminate_[iSubLam].DualityGap();
		}
		return gap;
	}

	void StepSize(ScalarType &primal_step, ScalarType &dual_step) {
		ScalarType temp_primal_step, temp_dual_step;
		sublaminate_[0].StepSize(primal_step, dual_step);
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			temp_primal_step = primal_step;temp_dual_step = dual_step;
			sublaminate_[iSubLam].StepSize(primal_step, dual_step);
			primal_step = primal_step < temp_primal_step? primal_step: temp_primal_step;
			dual_step = dual_step < temp_dual_step? dual_step: temp_dual_step;
		}
	}

	void HessianEval(Hessian_t &hessian)  {
		Hessian_t hessian_dual, hessian_dual_inv;
		MikiHessian_t sub_hessian,sub_hessian_inv;
		int iBlock, jBlock;

		hessian_dual = Hessian_t::Zero();
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			sub_hessian = MikiHessian_t::Zero();
			sublaminate_[iSubLam].Hessian(sub_hessian);
			sublaminate_cholesky_hessian_[iSubLam] = sub_hessian.llt();
			sub_hessian_inv = MikiHessian_t::Identity();
			sublaminate_cholesky_hessian_[iSubLam].solveInPlace(sub_hessian_inv);

			// should calculate only the symmetric part, will do it later
			for (int ilampar_type = 0; ilampar_type < NLAMPARTYPE; ++ilampar_type)
				for (int jlampar_type = 0; jlampar_type < NLAMPARTYPE; ++jlampar_type){
					iBlock = ilampar_type*NLAMPAR; jBlock = jlampar_type*NLAMPAR;
					hessian_dual.block(iBlock,jBlock,NLAMPAR,NLAMPAR) += sub_hessian_inv *clt.weight[ilampar_type][iSubLam]*clt.weight[jlampar_type][iSubLam];
				}
		}
		cholesky_dual_ = hessian_dual.llt();
		hessian_dual_inv = Hessian_t::Identity();
		cholesky_dual_.solveInPlace(hessian_dual_inv);
		hessian += hessian_dual_inv;
	}

	void CalulateResiduals(const Vector_t &var, Vector_t &residual, const ScalarType penalty){
		int iBlock;
		MikiVector_t miki_temp;
		Vector_t laminate_temp;

		duality_residual_ = -var;
		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			sublaminate_dvar_[iSubLam] = MikiVector_t::Zero();
			sublaminate_[iSubLam].CalculateResiduals(sublaminate_dvar_[iSubLam], penalty);
			miki_temp = sublaminate_dvar_[iSubLam];
			sublaminate_cholesky_hessian_[iSubLam].solveInPlace(miki_temp);
			for (int ilampar_type = 0; ilampar_type < NLAMPARTYPE ; ilampar_type++){
				iBlock = ilampar_type*NLAMPAR;
				duality_residual_.block(iBlock,0,NLAMPAR,1) += clt.weight[ilampar_type][iSubLam]*(sublaminate_var_[iSubLam]+miki_temp);
			}
		}
		laminate_temp = duality_residual_;
		cholesky_dual_.solveInPlace(laminate_temp);
		residual += laminate_temp;
	}

	void UpdateIncrements(const Vector_t &dvar)  {
		int iBlock;
		duality_residual_ = - duality_residual_ + dvar;
		cholesky_dual_.solveInPlace(duality_residual_);
		for (int iSubLam = 0; iSubLam < NSUBLAM; iSubLam++)  {
			for (int ilampar_type = 0; ilampar_type < NLAMPARTYPE; ++ilampar_type)  {
				iBlock = ilampar_type*NLAMPAR;
				sublaminate_dvar_[iSubLam] += clt.weight[ilampar_type][iSubLam]*duality_residual_.block(iBlock,0,NLAMPAR,1);
			}
			sublaminate_cholesky_hessian_[iSubLam].solveInPlace(sublaminate_dvar_[iSubLam]);
			sublaminate_[iSubLam].UpdateIncrements(sublaminate_dvar_[iSubLam]);
		}
	}

	void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step){

		for (int iSubLam = 0;iSubLam<NSUBLAM;iSubLam++){
			sublaminate_var_[iSubLam] += primal_step*sublaminate_dvar_[iSubLam];
			sublaminate_[iSubLam].UpdateVariables(primal_step,dual_step);
		}
	}
};

template<SubLaminate_t Type, bool isBalanced, bool isSymmetric, int NSUBLAM, typename ScalarType>
lampar::internal::cltweight<isSymmetric,NSUBLAM,ScalarType> lampar::lpfSDPCone<Type,isBalanced,isSymmetric,NSUBLAM,ScalarType> ::clt;
} /* namespace lampar */

#endif /* lpfSDPCone_H_ */
