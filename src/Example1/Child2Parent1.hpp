/*
 * Child2Parent.hpp
 *
 *  Created on: 9 mai 2015
 *      Author: Ramzi
 */

#ifndef EXAMPLE1_CHILD2PARENT1_HPP_
#define EXAMPLE1_CHILD2PARENT1_HPP_

#include "../Example1/Parent1.hpp"
#include <iostream>

using namespace std;
template<typename ScalarType = double>
class Child2Parent1 : public Parent1<ScalarType>{

private:
	string name_;

public:

	Child2Parent1(void) :Parent1<ScalarType>(){

	};

	~Child2Parent1(){

	};

	void message(void * arg, ScalarType value){
		string* name = (string*)arg;

		*name += " needs the book lol";
		name_ = *name;

	};

};




#endif /* EXAMPLE1_CHILD2PARENT1_HPP_ */
