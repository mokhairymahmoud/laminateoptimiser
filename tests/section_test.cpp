/*
 * section_test.cpp
 *
 *  Created on: 15 mai 2015
 *      Author: Ramzi
 */
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Section/section.hpp"
#include "Laminate/laminateSection.hpp"
using namespace Eigen;
//#include "Miki/Miki.hpp"
typedef Eigen::Matrix<double,Dynamic,Dynamic> Hessian;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> MyVector;
typedef Eigen::Matrix<double,4,4> Hessian_t;
typedef Eigen::Matrix<double,3,3> Hessian_t1;

void display (Eigen::Matrix<double,4,4> M){
	std::cout<<M;std::cout<<endl;
}
void test_helper(void * M){
	Hessian * pHessian = static_cast<Hessian *> (M);
	display(* pHessian);
}
template<int n>
void test(Eigen::Matrix<double,n,n> H){
	test_helper(&H);
}

int main()
{
	Hessian Z1;
	Hessian_t Z;
	Z=Hessian_t::Ones();
	Z1.resize(4,4);
	int n=Z1.size();
	for (int i=0;i<n;i++)
		Z1.data()[i] = 5;
	test(Z1);
cout<<"Ramzi";
/*
	vector <optsection::section *>::iterator pIS;
	vector <optsection::section *> pSection;
	pSection.resize(3);


	double thickness, width;
	double t[3]={5,10,20}, w[3]={5,10,20};
	int i=0;
*/
/*
	pSection[0] = new lampar::laminateSection<double>();
//	pSection[0]->setId(1);
	pSection[1] = new lampar::beamSection<double>;
	pSection[2] = new lampar::laminateSection<double>();
//	pSection[1]->setId(2);
//	pSection[0]->getSection()->volume(&thickness,&width);
//	pSection[1]->getSection()->volume(&thickness,&width);
	for (pIS=pSection.begin();pIS!=pSection.end();pIS++,i++){
		thickness =t[i]; width=w[i];
		(*pIS)->volume(&thickness,&width);
	}
*/
	typedef lampar::laminateSection<lampar::SingleMaterial,false> laminate1;
	typedef lampar::laminateSection<lampar::SingleMaterial,true> laminate2;
	optsection::section<> * pSection[6];
	laminate1 *pLaminate1;
	laminate2 *pLaminate2;
	int i;

	laminate1::Vector_t V1;
	double penalty1;
	double primal = 1, dual = 1;

	for (i=0;i<3;i++){
		pLaminate1 = new laminate1();
		V1=laminate1::Vector_t::Ones();V1(pLaminate1->Size-1) =1;
		pLaminate1->setBoundThickness(0,10);
		pLaminate1->Initialise(1,V1);
		pSection[i] = pLaminate1;
	}

	laminate2::Vector_t V2;

	for (i=3;i<6;i++){
		pLaminate2 = new laminate2();
		V2=laminate2::Vector_t::Ones();V2(pLaminate2->Size-1) =1;
		pLaminate2->setBoundThickness(-5,10);
		pLaminate2->Initialise(1,V2);
		pSection[i] = pLaminate2;
	}

	Hessian M;
	MyVector V,R;
	int Size, k;
	for (i=0;i<6;i++){
	//	M1=M1.Zero();
		Size = pSection[i]->getSize();
		M.resize(Size,Size);
		V.resize(Size);
		R.resize(Size);
		for (k=0;k<Size*Size;k++)
			M.data()[k] = 0;
		for (k=0;k<Size;k++){
			V.data()[k] = 5;
			R.data()[k] = 1;
		}
		pSection[i]->HessianEval(M);
		M = pSection[i]->getHessian();
		std::cout<<M;std::cout<<endl;
		penalty1 = pSection[i]->DualityGap();
		std::cout<<"penalty1 ="<<penalty1<<endl;
	//	R1=R1.Zero();
		pSection[i]->CalculateResiduals(V,R,0);
		R = pSection[i]->getResidual();
		std::cout<<"R1= "<<R;std::cout<<endl;
		pSection[i]->UpdateIncrements(R);
		primal = 1; dual = 1;
		pSection[i]->StepSize(primal,dual);
		std::cout<<"primal ="<<primal<<"   dual ="<<dual<<endl;
		pSection[i]->UpdateVariables(primal,dual);
	}

	return 0;
}




