/*
 * varsize.hpp
 *
 *  Created on: 11 avr. 2015
 *      Author: Ramzi
 */

#ifndef LAMINATE_VARSIZE_HPP_
#define LAMINATE_VARSIZE_HPP_


namespace lampar {

namespace internal {
template<bool isBalanced, bool isSymmetric>
class varsize{
public:
	static const int NVAR;
	static const int NVARTYPE;
	static const int SIZE;
};

template<>
class varsize<true, true>{
public:
	static const int NVAR = 2;
	static const int NVARTYPE = 2;
	static const int SIZE = NVAR*NVARTYPE;
};

template<>
class varsize<false, true>{
public:
	static const int NVAR = 4;
	static const int NVARTYPE = 2;
	static const int SIZE = NVAR*NVARTYPE;
};

template<>
class varsize<false, false>{
public:
	static const int NVAR = 4;
	static const int NVARTYPE = 3;
	static const int SIZE = NVAR*NVARTYPE;
};

template<>
class varsize<true, false> {
public:
	static const int NVAR = 2;
	static const int NVARTYPE = 3;
	static const int SIZE = NVAR*NVARTYPE;
};
} // namespace internal
} /* namespace lampar */

#endif /* LAMINATE_VARSIZE_HPP_ */
