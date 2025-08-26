/*
 * permutation_test.cpp
 *
 *  Created on: 6 déc. 2015
 *      Author: Ramzi
 */
#include <iostream>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <algorithm>
#include <vector>
typedef Eigen::SparseMatrix<double> SpMat;
using namespace Eigen;
using namespace std;

template <class T> void swapi ( T a, T b )
{
  T c;
  c = a;
  a=b;
  b=c;
}
template<typename T>
class CompareIndicesByAnotherVectorValues
{
	std::vector<T>* _values;
public:
	CompareIndicesByAnotherVectorValues(std::vector<T>* values) : _values(values) {}
public:
	bool operator() (const int& a, const int& b) const
		{ return (*_values)[a] > (*_values)[b]; }
};
template<int nVar, typename T>
class Compare
{
	Eigen::Matrix<T,nVar,1>* _values;
public:
	Compare(Eigen::Matrix<T,nVar,1>* values) : _values(values) {}
public:
	bool operator() (const int& a, const int& b) const
		{ return (_values->data())[a] > (_values->data())[b]; }
};
int myrandom (int i) { return std::rand()%i;}
int main() {
    PermutationMatrix<Dynamic,Dynamic> perm,perm1(3);
    int c, j;
    MatrixXd x = MatrixXd::Random(5,5);
   // std::sort(x.data(), x.data()+x.size());
    perm.resize(5);
    perm.setIdentity();
    int myints2[]= {61, 13 ,145, 15 };
    vector<int> v (myints2,myints2+4);
    int myints3[]= {0,1 ,2, 3 };
    vector<int> index (myints3,myints3+4);

    CompareIndicesByAnotherVectorValues<int> o(&v);
    sort(index.begin(), index.end(),o);

    Eigen::Matrix<double,5,1> V;
    V<<5.3,8.1,6.5,6.5,3;
    Compare<5,double> T(&V);
    sort(perm.indices().data(), perm.indices().data()+perm.indices().size(),T);
    SpMat y = x.sparseView();


    perm1.indices()=perm.indices().block(0,0,3,1);
  //  std::swap(perm.indices().data()[1],)perm.indices().data()[2]);
 /*   for (int i=dim-1; i>0; --i) {

    	j=myrandom(i+1);
    	c=perm.indices()[i];
    	perm.indices()[i] = perm.indices()[j];
    	perm.indices()[j] = c;

    }*/
    int xx=10, xy=20;                              // x:10 y:20
    std::swap(xx,xy);


    cout << "permutation\n" << perm.indices()<< endl << endl;
    cout << "permutation\n" << perm.toDenseMatrix();
    perm=perm.transpose();
    cout << "permutation\n" << perm.indices()<< endl << endl;
    cout << "permutation1\n" << perm1.indices()<< endl << endl;
    cout << "permutation\n" << perm.toDenseMatrix();
    cout << endl << endl;
    cout << "original x\n" << y;
    cout << endl << endl;
    cout << "permuted left x \n" << perm * y;
    cout <<endl << endl;
    cout << "permuted right x \n" << y * perm ;
    cout << endl << endl;
    cout << "permuted both x \n" << (perm * y).eval() * perm.transpose();
    cout <<endl << endl;
}



