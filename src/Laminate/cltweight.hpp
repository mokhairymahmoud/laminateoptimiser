/*
 * cltweight.hpp
 *
 *  Created on: 11 avr. 2015
 *      Author: Ramzi
 */

#ifndef LAMINATE_CLTWEIGHT_HPP_
#define LAMINATE_CLTWEIGHT_HPP_

#include <Eigen/Dense>

namespace lampar {
namespace internal {
// weight class
template<bool isSymmetric, int NSUBLAM, typename ScalarType = double>
class cltweight{
//public:
//	ScalarType weight[3][NSUBLAM];
};

template<int NSUBLAM, typename ScalarType>
class cltweight<true,NSUBLAM,ScalarType>{
public:
	ScalarType weight[2][NSUBLAM];
	cltweight(void){
		ScalarType z[NSUBLAM+1];
		ScalarType step = (ScalarType)1/(ScalarType)NSUBLAM;
		ScalarType z2_l,z2_lp1;
		int iSubLam;

		z[0] =0;
		for (int iz=1;iz<NSUBLAM+1;iz++){
			iSubLam = iz-1;
			z[iz] = z[iz-1] + step;
			z2_l = z[iz]*z[iz];
			z2_lp1 = z[iz-1]*z[iz-1];

			weight[0][iSubLam] = step;
			weight[1][iSubLam] = z2_l*z[iz]-z2_lp1*z[iz-1];
		}
	}
	~cltweight(void){};
};

template<int NSUBLAM, typename ScalarType>
class cltweight<false,NSUBLAM,ScalarType>{
public:
	ScalarType weight[3][NSUBLAM];
	cltweight(void){
		ScalarType z[NSUBLAM+1];
		ScalarType step = (ScalarType)2/(ScalarType)NSUBLAM;
		ScalarType z2_l,z2_lp1;
		int iSubLam;

		z[0] =-1.0;
		for (int iz=1;iz<NSUBLAM+1;iz++){
			iSubLam = iz-1;
			z[iz] = z[iz-1] + step;
			z2_l = z[iz]*z[iz];
			z2_lp1 = z[iz-1]*z[iz-1];

			weight[0][iSubLam] = 0.5*step;
			weight[1][iSubLam] = 0.5*(z2_l*z[iz]-z2_lp1*z[iz-1]);
			weight[2][iSubLam] = 0.5*(z2_l-z2_lp1);
		}
	}
	~cltweight(void){};
};

} // namespace internal
} /* namespace lampar */
#endif /* LAMINATE_CLTWEIGHT_HPP_ */
