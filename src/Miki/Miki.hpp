/*
 * Miki.h
 *
 *  Created on: 7 Jan 1005
 *      Author: root
 */

#ifndef MIKI_H_
#define MIKI_H_

#include <Eigen/Dense>
#include "../SDPA/SDPA.hpp"

namespace lampar {

enum SubLaminate_t{SingleMaterial, MultiMaterial};

template<SubLaminate_t Type, bool isBalanced, typename ScalarType = double>
class Miki {};

template<typename ScalarType>
class Miki<SingleMaterial, false, ScalarType> : public
      SDPA::Constraint<Miki<SingleMaterial, false, ScalarType>, 3, 4, ScalarType>
{
public:
	typedef typename SDPA::Constraint<Miki<SingleMaterial, false, ScalarType>, 3, 4, ScalarType>	Base_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, false, ScalarType>, 3, 4, ScalarType>::Matrix_t Matrix_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, false, ScalarType>, 3, 4, ScalarType>::Vector_t Vector_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, false, ScalarType>, 3, 4, ScalarType>::Hessian_t Hessian_t;
public:
	Miki() {};
	Miki(ScalarType dualscale) : Base_t(Vector_t::Zero(),dualscale) {}
	~Miki(){};
	void Initval(const Vector_t &var, Matrix_t &primal)
	{
		Eval(var,primal);
		primal += Matrix_t::Identity();
	}
	void Eval(const Vector_t &var, Matrix_t &primal) const
	{
		primal <<  0.0   , var(0), var(2),
				   var(0), var(1), var(3),
				   var(2), var(3),-var(1);
	}
	void AdEval(const Matrix_t &dual, Vector_t &residual) const
	{
		residual(0) += dual(1,0)+dual(0,1);
		residual(1) += dual(1,1)-dual(2,2);
		residual(2) += dual(2,0)+dual(0,2);
		residual(3) += dual(2,1)+dual(1,2);
	}
	void HessianEval(const Matrix_t &primal, const Matrix_t &dual, Hessian_t &hessian) const
	{
		hessian(0,0) += dual(0,0)*primal(1,1)+dual(1,1)*primal(0,0)+dual(0,1)*primal(0,1)+dual(1,0)*primal(1,0);
		hessian(0,1) =
		hessian(1,0) += dual(0,1)*primal(1,1)+dual(1,1)*primal(0,1)-dual(0,2)*primal(1,2)-dual(1,2)*primal(0,2);
		hessian(0,2) =
		hessian(2,0) += dual(0,0)*primal(1,2)+dual(1,2)*primal(0,0)+dual(0,1)*primal(0,2)+dual(0,2)*primal(0,1);
		hessian(0,3) =
		hessian(3,0) += dual(0,2)*primal(1,1)+dual(1,1)*primal(0,2)+dual(0,1)*primal(1,2)+dual(1,2)*primal(0,1);
		hessian(1,1) += dual(1,1)*primal(1,1)+dual(2,2)*primal(2,2)-dual(1,2)*primal(1,2)-dual(2,1)*primal(2,1);
		hessian(1,2) =
		hessian(2,1) += dual(0,1)*primal(1,2)+primal(1,0)*dual(2,1)-dual(0,2)*primal(2,2)-primal(2,0)*dual(2,2);
		hessian(1,3) =
		hessian(3,1) += primal(1,1)*dual(2,1)+dual(1,1)*primal(1,2)-dual(1,2)*primal(2,2)-primal(2,1)*dual(2,2);
        hessian(2,2) += dual(0,0)*primal(2,2)+primal(0,0)*dual(2,2)+dual(2,0)*primal(2,0)+dual(0,2)*primal(0,2);
        hessian(2,3) =
        hessian(3,2) += dual(1,0)*primal(2,2)+primal(0,1)*dual(2,2)+dual(2,0)*primal(2,1)+primal(0,2)*dual(1,2);
        hessian(3,3) += dual(1,1)*primal(2,2)+primal(1,1)*dual(2,2)+dual(2,1)*primal(2,1)+dual(1,2)*primal(1,2);
	}
};

template<typename ScalarType>
class Miki<SingleMaterial, true, ScalarType> : public
      SDPA::Constraint<Miki<SingleMaterial, true, ScalarType>, 3, 2, ScalarType>
{
public:
	typedef typename SDPA::Constraint<Miki<SingleMaterial, true, ScalarType>, 3, 2, ScalarType>	Base_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, true, ScalarType>, 3, 2, ScalarType>::Matrix_t Matrix_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, true, ScalarType>, 3, 2, ScalarType>::Vector_t Vector_t;
	typedef typename SDPA::Constraint<Miki<SingleMaterial, true, ScalarType>, 3, 2, ScalarType>::Hessian_t Hessian_t;
public:
	Miki() {};
	Miki(const ScalarType dualscale) : Base_t(Vector_t::Zero(),dualscale) {}
	~Miki(){};
	void Initval(const Vector_t &var, Matrix_t &primal)
	{
		Eval(var,primal);
		primal += Matrix_t::Identity();
	}
	void Eval(const Vector_t &var, Matrix_t &primal)
	{
		primal <<  0.0   , var(0), 0.0   ,
				   var(0), var(1), 0.0   ,
				   0.0   , 0.0   ,-var(1);
	}
	void AdEval(const Matrix_t &dual, Vector_t &residual)
	{
		residual(0) += dual(1,0)+dual(0,1);
		residual(1) += dual(1,1)-dual(2,2);
	}
	void HessianEval(const Matrix_t &primal, const Matrix_t &dual, Hessian_t &hessian)
	{
		hessian(0,0) += dual(0,0)*primal(1,1)+dual(1,1)*primal(0,0)+dual(0,1)*primal(0,1)+dual(1,0)*primal(1,0);
		hessian(0,1) =
		hessian(1,0) += dual(0,1)*primal(1,1)+dual(1,1)*primal(0,1);
		hessian(1,1) += dual(1,1)*primal(1,1)+dual(2,2)*primal(2,2);
	}
};

}  /* namespace lampar */

#endif /* MIKI_H_ */
