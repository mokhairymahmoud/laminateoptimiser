/*
 * Miki_test.cpp
 *
 *  Created on: 7 Jan 2015
 *      Author: root
 */

#include "Miki/Miki.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Cholesky>


int main()
{
	lampar::Miki<lampar::SingleMaterial,false> feasibility(1.0);
	Eigen::Matrix<double, 4, 1> g, dx, x;
	Eigen::Matrix<double, 4, 4> B;
	Eigen::Matrix<double, 3, 3> X;
	Eigen::LLT<Eigen::Matrix<double, 4, 4> > LB;
	double d_g, penalty, d_g_c, predictor, corrector;
	double pstep = 1.0, dstep = 1.0;

	g(0) = 1;
	g(1) = 0;
	g(2) = 0;
	g(3) = 0;
    x = Eigen::Matrix<double, 4, 1>::Zero();
    int iter = 0;
	//feasibility.Initialise(x, X);
    do
    {
		dx = -g; B = Eigen::Matrix<double, 4, 4>::Zero();
		d_g = feasibility.DualityGap();
		penalty = d_g/3.0;
		predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
		feasibility.CalculateResiduals(dx,predictor*penalty);
		feasibility.Hessian(B);
		//std::cout << "Hessian: "<< std::endl << B << std::endl;
		LB = B.llt();
		LB.solveInPlace(dx);
		//std::cout << "predictor increment: " << std::endl << dx.transpose() << std::endl;
		feasibility.UpdateIncrements(dx);
		d_g_c = feasibility.DualityGap();
		corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(d_g,d_g_c);
		dx = -g;
		feasibility.CalculateResiduals(dx,corrector*penalty);
		//std::cout << "corrector residual: " << std::endl << dx.transpose() << std::endl;
		LB.solveInPlace(dx);
		std::cout << "correction factor: " << corrector << std::endl;
		//std::cout << "corrector increment: " << std::endl << dx.transpose() << std::endl;
		feasibility.UpdateIncrements(dx);
		pstep = 1.0; dstep = 1.0;
		feasibility.StepSize(pstep,dstep);
		feasibility.UpdateVariables(pstep,dstep);
		x += pstep*dx;
		std::cout << "iteration: " << iter+1 << '\t' << pstep << '\t' << dstep << '\t' << d_g <<std::endl;
	} while(d_g > 1.0e-10 and ++iter<20);


	std::cout <<"x = "<<x<<std::endl;
	return 0;
}


