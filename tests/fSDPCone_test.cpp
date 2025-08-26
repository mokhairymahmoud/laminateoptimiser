/*
 * fSDPCone_test.cpp
 *
 *  Created on: 15 avr. 2015
 *      Author: Ramzi
 */
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include "fSDPCone/fSDPCone.hpp"


using namespace std;
typedef typename Eigen::Matrix<double, 3, 3> matrix3;
int main()
{
	matrix3 M0 , Mi;
	vector<matrix3> M;

	M0 = matrix3::Identity();
	Mi <<  0.0, 1.0, 0.0,
		   1.0, 0.0, 0.0,
		   0.0, 0.0, 0.0;
	M.push_back(Mi);
	Mi <<  0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0,-1.0;
	M.push_back(Mi);
	Mi <<  0.0, 0.0, 1.0,
		   0.0, 0.0, 0.0,
		   1.0, 0.0, 0.0;
	M.push_back(Mi);
	Mi <<  0.0, 0.0, 0.0,
		   0.0, 0.0, 1.0,
		   0.0, 1.0, 0.0;
	M.push_back(Mi);

	lampar::fSDPCone<3,4> feasibility;
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
    feasibility.InitConvexSet(M0, M);
	feasibility.Initialise(1.0,x);
    do
    {
		dx = -g; B = Eigen::Matrix<double, 4, 4>::Zero();
		d_g = feasibility.DualityGap();
		penalty = d_g/3.0;
		predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
		feasibility.CalculateResiduals(dx,predictor*penalty);
		feasibility.Hessian(B);
		std::cout << "Hessian: "<< std::endl << B << std::endl;
		LB = B.llt();
		LB.solveInPlace(dx);
		std::cout << "predictor increment: " << std::endl << dx.transpose();
		std::cout<< std::endl;
		feasibility.UpdateIncrements(dx);
		d_g_c = feasibility.DualityGap();
		corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(d_g,d_g_c);
		dx = -g;
		feasibility.CalculateResiduals(dx,corrector*penalty);
		std::cout << "corrector residual: " << std::endl << dx.transpose();
		std::cout<< std::endl;
		LB.solveInPlace(dx);
		std::cout << "correction factor: " << corrector << std::endl;
		std::cout << "corrector increment: " << std::endl << dx.transpose();
		std::cout<< std::endl;
		feasibility.UpdateIncrements(dx);
		pstep = 1.0; dstep = 1.0;
		feasibility.StepSize(pstep,dstep);
		feasibility.UpdateVariables(pstep,dstep);
		x += pstep*dx;
		std::cout << "iteration: " << iter+1 << '\t' << pstep << '\t' << dstep << '\t' << d_g <<std::endl;
	} while(d_g > 1.0e-10 and ++iter<20);
    std::cout <<"x = "<<x<<endl;
	return 0;
}




