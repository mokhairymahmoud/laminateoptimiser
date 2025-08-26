/*
 * minmaxprob_test.cpp
 *
 *  Created on: 3 oct. 2015
 *      Author: Ramzi
 */

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "GlobalOptimiser/approxFunction.hpp"
#include "GlobalOptimiser/scminmaxProb.hpp"
#include "BoundSDP/boundSDP.hpp"
#include <math.h>

typedef typename Eigen::Matrix<double, 3, 1> vector3;
typedef typename Eigen::Matrix<double, 8, 1> vector8;
typedef typename Eigen::Matrix<double, 8, 8> matrix8;
typedef typename Eigen::Matrix<double, 3, 3> matrix3;
typedef typename globopt::approxFunction<4,2> approxFunc1;
typedef typename lampar::boundSDP<lampar::Lower> sideConstraint;

int main()
{
	vector3 a;
	matrix3 A,LA;
	int i;

	a <<  5.3, 1.0, 3.0;
	A <<  4.0, -1.0, 2.0,
		  -1.0, 6, 0,
		  2.0, 0.0, 5.0;
	Eigen::LLT<matrix3> LAtemp;
	LAtemp = A.llt();
	LA = LAtemp.matrixL();
	std::cout<<"LA="<<std::endl;
	std::cout<<LA;
	std::cout<<std::endl;


	approxFunc1::Vector_v x0;
	approxFunc1 *approxObject;
	approxObject= new approxFunc1();
	sideConstraint *sideObject;
	sideObject = new sideConstraint[2];

	sideObject[0].setBound(0);sideObject[1].setBound(0);
//	x0<<.788,.41;
	x0<<1,1;
	globopt::responseInterface<approxFunc1,sideConstraint> mmprob;

	mmprob.scMinMaxProb(approxObject, sideObject);
	for (i=0;i<20;i++){
		std::cout<<"Solver iteration = "<<i<<std::endl;
		approxObject->InitialiseProblem1(x0);
		mmprob.Solver(x0);
	}

	return 0;
}








