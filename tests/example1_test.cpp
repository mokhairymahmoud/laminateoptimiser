/*
 * example1_test.cpp
 *
 *  Created on: 9 mai 2015
 *      Author: Ramzi
 */


#include <iostream>
#include "Example1/Parent1.hpp"
#include "Example1/Child1Parent1.hpp"
#include "Example1/Child2Parent1.hpp"


int main()
{
	Parent1<double> *p1;
	Child1Parent1<double> c1p1;
	Child2Parent1<double> c2p1;
	int a = 19;
	string name="Ramzi";
	p1 = &c1p1;
	p1->message(&a,5);
	std::cout <<"x = "<<a<<std::endl;
	p1 = &c2p1;
	p1->message(&name,5);
	std::cout <<"x = "<<name<<std::endl;
	return 0;
}




