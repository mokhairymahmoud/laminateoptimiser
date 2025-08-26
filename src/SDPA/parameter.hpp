/*
 * parameter.hpp
 *
 *  Created on: 2 mai 2015
 *      Author: Ramzi
 */

#ifndef SDPA_PARAMETER_HPP_
#define SDPA_PARAMETER_HPP_

#include <iostream>

namespace SDPA {

template<typename ScalarType = double>
class Parameter  {
public:
	static ScalarType Step_Size_Control(const ScalarType eval) {
		ScalarType step = -eval/0.95;
		step = step>1?step:1.0;
		return 1.0/step;
	}
	static ScalarType Predictor_Duality_Reduction()
	{
		return 0.0;
	}
	static ScalarType Corrector_Duality_Reduction(ScalarType oldgap, ScalarType newgap)
	{
		ScalarType step = newgap<oldgap?newgap/oldgap:1.0;
		step*=step;
		return step>0.1?step:0.1;
	}
};
} // end namespace SDPA




#endif /* SDPA_PARAMETER_HPP_ */
