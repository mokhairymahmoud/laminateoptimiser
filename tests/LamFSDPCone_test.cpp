/*
 * LamFSDPCone_test.cpp
 *
 *  Created on: 30 avr. 2015
 *      Author: Ramzi
 */
#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <Eigen/Core>

#include "LamFSDPCone/lpfSDPCone.hpp"
#include "Miki/Miki.hpp"

using namespace std;
typedef typename Eigen::Matrix<double, 3, 3> matrix3;

#define  nSubLam  10
#define isBalanced  false
#define isSymmetric  true
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

	Eigen::Matrix2d mat;
	  mat << 1, 2,
	         3, 4;
	  cout << "Here is mat.sum():       " << mat.sum()       << endl;
	  cout << "Here is mat.prod():      " << mat.prod()      << endl;
	  cout << "Here is mat.mean():      " << mat.mean()      << endl;
	  cout << "Here is mat.minCoeff():  " << Mi.minCoeff()  << endl;
	  cout << "Here is mat.maxCoeff():  " << mat.maxCoeff()  << endl;
	  cout << "Here is mat.trace():     " << mat.trace()     << endl;

/*
    lampar::lpfSDPCone<lampar::SingleMaterial, isBalanced,isSymmetric,nSubLam,double> lam;
    lampar::lpfSDPCone<lampar::SingleMaterial, isBalanced,isSymmetric,nSubLam,double>::Hessian_t H;
    lampar::lpfSDPCone<lampar::SingleMaterial, isBalanced,isSymmetric,nSubLam,double>::Cholesky_t LH;
    lampar::lpfSDPCone<lampar::SingleMaterial, isBalanced,isSymmetric,nSubLam,double>::Vector_t v,dv,r,g;
	double d_g, penalty, d_g_c, predictor, corrector;
	double pstep = 1.0, dstep = 1.0;
	int nGap = 1*nSubLam*3;
	int iter = 0;
	string fileName = "test1.txt";
	fstream os;
	double R= 1, coordX, coordY,theta, step;
	int nbr= 100;

	step = 2*M_PI/(double)nbr;
	os.open(fileName.c_str(), ios_base::out);
	os<<"M = [";
//	for (int i=0;i<=nbr;i++)
	int i = 1;
	{
		iter = 0;
		theta = step*i;
		coordX = R*cos(theta);
		coordY = R*sin(theta);

		v = v.Zero();

		g <<0,0,coordX,coordY,0,0,0,0;

		cout<<"g= "<<g.transpose();
		cout<<endl;
		os<<g.transpose()<<" ";
		lam.Initialise(M0,M,1.0);
		do
		{
			r = -g;
			H = H.Zero();
			d_g = lam.DualityGap() ;
			penalty = d_g/nGap;
			lam.HessianEval(H);
			predictor = SDPA::Parameter<>::Predictor_Duality_Reduction();
			lam.CalulateResiduals(v,r,predictor*penalty);
	//		std::cout << "Hessian: "<< std::endl << H << std::endl;
			LH = H.llt();
			dv = r;
			LH.solveInPlace(dv);
	//		std::cout << "predictor increment: " <<  dv.transpose();
	//		std::cout<< std::endl;
			lam.UpdateIncrements(dv);
			d_g_c = lam.DualityGap() ;
			corrector = SDPA::Parameter<>::Corrector_Duality_Reduction(d_g,d_g_c);
			r = -g;
			lam.CalulateResiduals(v,r,corrector*penalty);
			dv = r;
	//		std::cout << "corrector residual: " << dv.transpose();
	//		std::cout<< std::endl;
			LH.solveInPlace(dv);
	//		std::cout << "correction factor: " << corrector << std::endl;
	//		std::cout << "corrector increment: " <<  dv.transpose();
	//		std::cout<< std::endl;
			lam.UpdateIncrements(dv);
			pstep = 1.0; dstep = 1.0;
			lam.StepSize(pstep,dstep);
			lam.UpdateVariables(pstep,dstep);
			v += pstep*dv;
			std::cout << "iteration: " << iter+1 << '\t' << pstep << '\t' << dstep << '\t' << d_g <<std::endl;
		} while(d_g > 1.0e-10 and ++iter<20);

		std::cout << "v = "<<v<<std::endl;
		os<<v.transpose()<<";";
		os<<endl;
	}
	os<<"];"<<endl;
	os.close();*/
	return 0;
}








