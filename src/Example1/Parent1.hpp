/*
 * Example1.hpp
 *
 *  Created on: 9 mai 2015
 *      Author: Ramzi
 */

#ifndef EXAMPLE1_PARENT1_HPP_
#define EXAMPLE1_PARENT1_HPP_


#include <iostream>


template<typename ScalarType = double>
class Parent1  {


public:

	Parent1(void) {

	};

	virtual ~Parent1(){

	};

	virtual void message(void * arg, ScalarType value){

	};

};





#endif /* EXAMPLE1_PARENT1_HPP_ */
