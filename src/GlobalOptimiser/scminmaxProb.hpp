/*
 * scminmaxProb.hpp
 *
 *  Created on: 3 oct. 2015
 *      Author: Ramzi
 */

#ifndef GLOBALOPTIMISER_SCMINMAXPROB_HPP_
#define GLOBALOPTIMISER_SCMINMAXPROB_HPP_
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <vector>
#include "../GlobalOptimiser/approxFunction.hpp"
#include "../SDPA/parameter.hpp"

using namespace std;
using namespace Eigen;

namespace globopt {

template<typename matrixType>
class Compare
{
	matrixType* _values;
public:
	Compare(matrixType* values) : _values(values) {}
	~Compare() {}
public:
	bool operator() (const int& a, const int& b) const
		{ return (_values->data())[a] > (_values->data())[b]; }
};

template<typename approximationFunction, typename sideConstraint, typename ScalarType = double>
class responseInterface
{
private:
	int nVar;
	int nResp;
	int nRespReduced;
	approximationFunction *approximationObject;
	sideConstraint *sideObject;
public:
	typedef typename approximationFunction::Vector_r Vector_r;
	typedef typename approximationFunction::Vector_v Vector_v;
	typedef typename approximationFunction::Matrix_t Matrix_t;
	typedef typename approximationFunction::Hessian_t Hessian_t;
	typedef typename approximationFunction::Hessian_r Hessian_r;
	typedef Eigen::LLT<Hessian_t> Cholesky_t;

	typedef typename approximationFunction::reduced_Hessian_r reduced_Hessian_r;
	typedef typename approximationFunction::reduced_Vector_r reduced_Vector_r;
	typedef Eigen::LLT<reduced_Hessian_r> reduced_Cholesky_r;

	typedef typename approximationFunction::Matrix_r_i Matrix_r_i;
	typedef typename approximationFunction::Vector_r_i Vector_r_i;

private:
	Cholesky_t cholesky__primal_hessian_;
	reduced_Cholesky_r reduced_cholesky_duality_hessian_;
	reduced_Vector_r reduced_duality_hessian_vector_;
	ScalarType exact_system_matrix_11, exact_system_matrix_12, exact_system_matrix_22;
	Vector_r booleanVector_;
	PermutationMatrix<Dynamic,Dynamic> permutationMatrix_;
	reduced_Vector_r reduced_BooleanVector_;

	Vector_r dual_ , ddual_ ;
	Vector_r slack_, dslack_;
	Vector_v primal_, dprimal_;
	ScalarType objective_primal_, objective_dprimal_;

	reduced_Vector_r reduced_ddual_vector_;
	ScalarType reduced_ddual_scalar_;
	bool verbose_ = true;
	int lastIterationCount_ = 0;


private:

	ScalarType FindStep(const ScalarType &val, const ScalarType &dval) const {
		return dval/val;
	}
	template<typename matrixType>
	ScalarType FindStep(const matrixType &val, const matrixType &dval) const {
		return (dval.cwiseProduct(val.cwiseInverse())).minCoeff();
	}

	ScalarType maxResponses(Vector_r resp, Vector_r booleanVector){
		int iResp , iRespFirst =0;
		ScalarType maxResp;

		while ((int)booleanVector(iRespFirst) ==0)
			iRespFirst++;
		maxResp = resp(iRespFirst);
		for (iResp = iRespFirst;iResp<nResp;iResp++){
			if ((int)booleanVector(iResp) !=0){
				if (resp(iResp)>maxResp)
					maxResp = resp(iResp);
			}
		}
		return maxResp;
	}
	void creatPermutationMatrixOld(Vector_r resp, Vector_r booleanVector, Hessian_r &permutation){
		int i,j,inc,midIndex;
		int index[nResp], indexTemp[nResp];
		ScalarType midResp;

		for(inc = 0;inc<nResp;inc++)
			indexTemp[inc] = inc;
		permutation = Hessian_r::Zero();
		for(inc = nResp/2;inc>0;inc/=2){
			for(i = inc;i< nResp;i++){
				for(j=i-inc;j>=0;j-=inc){
					if(resp(j+inc)<=resp(j))
						break;
					else{
						midResp = resp(j);
						resp(j) = resp(j+inc);
						resp(j+inc) = midResp;
						midIndex = indexTemp[j];
						indexTemp[j] = indexTemp[j+inc];
						indexTemp[j+inc] = midIndex;
					}
				}
			}
		}
		inc=0;
		while (booleanVector(indexTemp[inc])!=1){
			inc++;
		}
		for (i=0;i<inc;i++){
			index[i] = indexTemp[i];
		}
		for (i=inc+1;i<nResp;i++){
			index[i-1] = indexTemp[i];
		}
		index[nResp-1] = indexTemp [inc];

		for(inc = 0;inc<nResp;inc++) {
			permutation(inc,index[inc]) = 1;
		}

	}
	void creatPermutationMatrix(Vector_r resp){
		int inc =0, indexObj;
		globopt::Compare<Vector_r> operatorCompare(&resp);

		permutationMatrix_.setIdentity();
		sort(permutationMatrix_.indices().data(), permutationMatrix_.indices().data()+permutationMatrix_.indices().size(),operatorCompare);
		while (booleanVector_(permutationMatrix_.indices().data()[inc])!=1) inc++;
		indexObj = permutationMatrix_.indices().data()[inc];
		permutationMatrix_.indices().block(inc,0,nRespReduced-inc,1)=permutationMatrix_.indices().block(inc+1,0,nRespReduced-inc,1);
		permutationMatrix_.indices().data()[nRespReduced]=indexObj;
		permutationMatrix_ = permutationMatrix_.transpose();
	}
public:
	void scMinMaxProb() {

	}
	void scMinMaxProb(approximationFunction *appObj, sideConstraint *sideObj) {
		approximationObject = appObj;
		sideObject = sideObj;
		nVar = approximationObject->NVAR;
		nResp = approximationObject->NRESP;
		nRespReduced = nResp-1;
		permutationMatrix_.resize(nResp);
	}

	~responseInterface() {
	}
	void setVerbose(const bool verbose) {
		verbose_ = verbose;
	}
	int getLastIterationCount() const {
		return lastIterationCount_;
	}
	void Initialise(ScalarType dual_scale, Vector_v &var)  {
		Vector_r booleanVector;
		ScalarType dimObjSet;
		int iVar;

		// temporary
		primal_ = var;  // <<.788,.41; // temporary
		// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].Initialise(dual_scale,var(iVar));

		approximationObject->Eval(var, slack_);
		dimObjSet = approximationObject->getBooleanVector(booleanVector);
		// end change
		objective_primal_ = maxResponses(slack_,booleanVector)*dimObjSet/dimObjSet;
		slack_ = -slack_+objective_primal_*booleanVector+1.0*booleanVector;

		objective_primal_ = 1;
/*		slack_ = Vector_r::Ones();
		dual_ = booleanVector/dimObjSet+
				(approximationFunction::Vector_r::Ones()-booleanVector)/(nResp-dimObjSet);
*/
		dual_ = dual_scale*slack_.cwiseInverse();
		dslack_ = approximationFunction::Vector_r::Zero();
		ddual_ = approximationFunction::Vector_r::Zero();
		dprimal_ = approximationFunction::Vector_v::Zero();
		objective_dprimal_ = 0;
	}
	void HessianEval(Matrix_t gradient_matrix, Hessian_t &primal_hessian)  {
		Hessian_r duality_hessian;
		Vector_r duality_hessian_tempVect, booleanVector;
		reduced_Hessian_r reduced_duality_hessian;
		reduced_Vector_r reduced_duality_hessian_vector_temp;
		reduced_Vector_r reduced_BooleanVector_temp;
		Matrix_t gradient_matrix_temp;
		int iVar;

		duality_hessian = approximationFunction::Hessian_r::Zero();
	// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].Hessian(primal_hessian(iVar,iVar));
	// end change

		cholesky__primal_hessian_ =  primal_hessian.llt();
		gradient_matrix_temp = gradient_matrix;
		cholesky__primal_hessian_.solveInPlace(gradient_matrix_temp);

		duality_hessian_tempVect = slack_.cwiseProduct(dual_.cwiseInverse());
		duality_hessian.diagonal() = duality_hessian_tempVect;

		duality_hessian += gradient_matrix.transpose()*gradient_matrix_temp;
		duality_hessian = (permutationMatrix_* duality_hessian).eval()*permutationMatrix_.transpose();
		booleanVector = permutationMatrix_*booleanVector_;
		reduced_duality_hessian = duality_hessian.block(0,0,nRespReduced,nRespReduced);
		reduced_duality_hessian_vector_ = duality_hessian.block(0,nRespReduced,nRespReduced,1);
		exact_system_matrix_11 = duality_hessian(nRespReduced,nRespReduced);
		reduced_BooleanVector_ = booleanVector.block(0,0,nRespReduced,1);

		reduced_cholesky_duality_hessian_ = reduced_duality_hessian.llt();

		reduced_duality_hessian_vector_temp = reduced_duality_hessian_vector_;
		reduced_BooleanVector_temp = reduced_BooleanVector_;

		reduced_cholesky_duality_hessian_.solveInPlace(reduced_duality_hessian_vector_temp);
		reduced_cholesky_duality_hessian_.solveInPlace(reduced_BooleanVector_temp);

		exact_system_matrix_11 -= reduced_duality_hessian_vector_.transpose()*reduced_duality_hessian_vector_temp;
		exact_system_matrix_12 = (ScalarType)1.0-reduced_duality_hessian_vector_.transpose()*reduced_BooleanVector_temp;
		exact_system_matrix_22 = reduced_BooleanVector_.transpose()*reduced_BooleanVector_temp;
	}
	void CalculateResiduals(const ScalarType penalty, Vector_r response, Matrix_t gradient_matrix){
		int iVar;
		Vector_v primal_residual_temp;
		reduced_Vector_r reduced_duality_residual_vector_temp;

		dslack_ = penalty*Vector_r::Ones()-slack_.cwiseProduct(dual_)-dslack_.cwiseProduct(ddual_);
		ddual_ = response- booleanVector_*objective_primal_+ slack_;
		ddual_ +=  dslack_.cwiseProduct(dual_.cwiseInverse());
		objective_dprimal_ = -(ScalarType)1.0+booleanVector_.transpose()*dual_;

		dprimal_ = -gradient_matrix*dual_;
	// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].CalculateResiduals(dprimal_(iVar),penalty);
	// end change

		primal_residual_temp = dprimal_;
		cholesky__primal_hessian_.solveInPlace(primal_residual_temp);
		ddual_ += gradient_matrix.transpose()*primal_residual_temp;
		ddual_ = permutationMatrix_*ddual_;

		reduced_ddual_vector_ = ddual_.block(0,0,nRespReduced,1);
		reduced_ddual_scalar_ = ddual_(nRespReduced);

		reduced_duality_residual_vector_temp = reduced_ddual_vector_;
		reduced_cholesky_duality_hessian_.solveInPlace(reduced_duality_residual_vector_temp);
		reduced_ddual_scalar_ += reduced_duality_hessian_vector_.transpose()*reduced_duality_residual_vector_temp;
		objective_dprimal_ += reduced_BooleanVector_.transpose()*reduced_duality_residual_vector_temp;
	}
	void UpdateIncrements(Matrix_t gradient_matrix)	{
		ScalarType tempVal;
		int iVar;

		tempVal = (ScalarType)1.0/(exact_system_matrix_12*exact_system_matrix_12+exact_system_matrix_11*exact_system_matrix_22);
		objective_dprimal_= tempVal*(exact_system_matrix_11*objective_dprimal_+exact_system_matrix_12*reduced_ddual_scalar_);
		reduced_ddual_scalar_ = (-exact_system_matrix_12*objective_dprimal_+reduced_ddual_scalar_)/exact_system_matrix_11;
	//	dreduced_duality_scalar = tempVal*(exact_system_matrix_22*reduced_ddual_scalar_+exact_system_matrix_12*objective_dprimal_);


		reduced_ddual_vector_ -= reduced_duality_hessian_vector_*reduced_ddual_scalar_+reduced_BooleanVector_*objective_dprimal_;
		reduced_cholesky_duality_hessian_.solveInPlace(reduced_ddual_vector_);

		ddual_.block(0,0,nRespReduced,1) = reduced_ddual_vector_;
		ddual_(nRespReduced) = reduced_ddual_scalar_;
		ddual_ = permutationMatrix_.transpose()*ddual_;

		dprimal_ = dprimal_-gradient_matrix*ddual_;
		cholesky__primal_hessian_.solveInPlace(dprimal_);
	// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].UpdateIncrements(dprimal_(iVar));
	// end change
		dslack_ = dslack_-slack_.cwiseProduct(ddual_);
		dslack_ = dslack_.cwiseProduct(dual_.cwiseInverse());
	}
	ScalarType DualityGap() const {
		ScalarType temp;
		int iVar;

		temp = (slack_+dslack_).transpose()*(dual_+ddual_);
	// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			temp += sideObject[iVar].DualityGap();
	// end change
		return temp;
	}
	void StepSize(ScalarType &primal_step, ScalarType &dual_step) const{
		ScalarType eval, step;
		int iVar;

	// replace this once we add lamination
		eval = FindStep(primal_,dprimal_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].StepSize(primal_step,dual_step);
	// end replace
		eval = FindStep(slack_,dslack_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
		eval = FindStep(objective_primal_,objective_dprimal_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		primal_step = primal_step < step? primal_step: step;
		eval = FindStep(dual_,ddual_);
		step = SDPA::Parameter<ScalarType>::Step_Size_Control(eval);
		dual_step = dual_step < step? dual_step: step;
	}
	void UpdateVariables(const ScalarType primal_step,const ScalarType dual_step){
		int iVar;

		primal_ += primal_step*dprimal_;  // change this once we add lamination
		slack_ += primal_step*dslack_;
		objective_primal_ += primal_step*objective_dprimal_;
		dual_ += dual_step*ddual_;
	// change this once we add lamination
		for (iVar=0;iVar<nVar;iVar++)
			sideObject[iVar].UpdateVariables(primal_step,dual_step);
	// end change

		objective_dprimal_ = 0;
		ddual_ = Vector_r::Zero();
		dslack_ = Vector_r::Zero();
		dprimal_ = Vector_v::Zero(); // change this once we add lamination
	}
	void Solver(Vector_v &var){
		Matrix_t gradient_matrix;
		Hessian_t primal_hessian;
		Vector_r response;

		double d_g, penalty, d_g_c, predictor, corrector;
		double pstep = 1.0, dstep = 1.0;
	    int iter = 0;
	    lastIterationCount_ = 0;

	    primal_ = var;
	    Initialise(1.0, primal_);
	    do
	    {
			d_g = DualityGap();
			penalty = d_g/(nResp+nVar);
			predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
			primal_hessian = approximationFunction::Hessian_t::Zero();
			gradient_matrix = approximationFunction::Matrix_t::Zero();
		// change this once we add lamination
			approximationObject->Eval(primal_, dual_, response, gradient_matrix, primal_hessian);
			approximationObject->getBooleanVector(booleanVector_);
			creatPermutationMatrix(response);
		// end change
			HessianEval(gradient_matrix, primal_hessian);
			CalculateResiduals(predictor*penalty, response, gradient_matrix);

			UpdateIncrements(gradient_matrix);
			d_g_c = DualityGap();
			corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(d_g,d_g_c);

			CalculateResiduals(corrector*penalty, response, gradient_matrix);

			if (verbose_)
				std::cout << "correction factor: " << corrector << std::endl;

			UpdateIncrements(gradient_matrix);
			pstep = 1.0; dstep = 1.0;
			StepSize(pstep,dstep);
			UpdateVariables(pstep,dstep);

			lastIterationCount_ = iter + 1;
			if (verbose_)
				std::cout << "iteration: " << iter+1 << '\t' << pstep << '\t' << dstep << '\t' << d_g <<std::endl;
		} while(d_g > 1.0e-10 and ++iter<20);


		if (verbose_) {
			std::cout <<"x = "<<primal_.transpose();
			std::cout<<std::endl;
		}
		approximationObject->Eval(primal_, response);
		if (verbose_) {
			std::cout <<"response_ = "<<response.transpose();
			std::cout<<std::endl;
		}
		var = primal_;
	}
};

}  /* namespace globopt */



#endif /* GLOBALOPTIMISER_SCMINMAXPROB_HPP_ */
