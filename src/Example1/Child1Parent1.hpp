/*
 * Child1Parent1.hpp
 *
 *  Created on: 9 mai 2015
 *      Author: Ramzi
 */

#ifndef EXAMPLE1_CHILD1PARENT1_HPP_
#define EXAMPLE1_CHILD1PARENT1_HPP_

#include "../Example1/Parent1.hpp"
#include <iostream>


template<typename ScalarType = double>
class Child1Parent1 : public Parent1<ScalarType>{

private:
	int age_;

public:

	Child1Parent1(void):Parent1<ScalarType>()  {
		age_ =0;
	};

	~Child1Parent1(){

	};

	void message(void * arg, ScalarType value){
		int* age = (int*)arg;

		*age *= 2;
		age_ = *age;
	};

};



#endif /* EXAMPLE1_CHILD1PARENT1_HPP_ */
