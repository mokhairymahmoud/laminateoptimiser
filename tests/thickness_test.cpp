/*
 * thickness_test.cpp
 *
 *  Created on: 2 mai 2015
 *      Author: Ramzi
 */

#include <iostream>
#include "BoundSDP/boundSDP.hpp"


int main()
{
	lampar::boundSDP<lampar::Upper> feasibility_U(100);
	lampar::boundSDP<lampar::Lower> feasibility_L(-100);
	double g, dx, x;
	double B;


	double d_g, penalty, d_g_c, predictor, corrector;
	double pstep = 1.0, dstep = 1.0;

	x = -40;
    int iter = 0;
	feasibility_U.Initialise(1.0, x);
	feasibility_L.Initialise(1.0, x);
    do
    {
    	g = x*x-1*x-2;
		dx = -g; B = 2*x-1;
		d_g = feasibility_U.DualityGap()+feasibility_L.DualityGap();
		penalty = d_g/2.0;
		predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
		feasibility_U.CalculateResiduals(dx,predictor*penalty);
		feasibility_L.CalculateResiduals(dx,predictor*penalty);
		feasibility_U.Hessian(B);
		feasibility_L.Hessian(B);
		//std::cout << "Hessian: "<< std::endl << B << std::endl;
		dx = dx /B;
		//std::cout << "predictor increment: " << std::endl << dx.transpose() << std::endl;
		feasibility_U.UpdateIncrements(dx);
		feasibility_L.UpdateIncrements(dx);
		d_g_c = feasibility_U.DualityGap()+feasibility_L.DualityGap();
		corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(d_g,d_g_c);
		dx = -g;
		feasibility_U.CalculateResiduals(dx,corrector*penalty);
		feasibility_L.CalculateResiduals(dx,corrector*penalty);
		//std::cout << "corrector residual: " << std::endl << dx.transpose() << std::endl;
		dx = dx /B;
		std::cout << "correction factor: " << corrector << std::endl;
		//std::cout << "corrector increment: " << std::endl << dx.transpose() << std::endl;
		feasibility_U.UpdateIncrements(dx);
		feasibility_L.UpdateIncrements(dx);
		pstep = 1.0; dstep = 1.0;
		feasibility_U.StepSize(pstep,dstep);
		feasibility_L.StepSize(pstep,dstep);
		feasibility_U.UpdateVariables(pstep,dstep);
		feasibility_L.UpdateVariables(pstep,dstep);
		x += pstep*dx;
		std::cout << "iteration: " << iter+1 << '\t' << pstep << '\t' << dstep << '\t' << d_g <<std::endl;
	} while(d_g > 1.0e-10 and ++iter<20);


	std::cout <<"x = "<<x<<std::endl;
	return 0;
}


