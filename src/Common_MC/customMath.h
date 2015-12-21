#ifndef CUSTOMMATH_H_
#define CUSTOMMATH_H_

// STL containers and algorithms
#include <vector>
#include <set>
#include <algorithm>

typedef unsigned int uint;

// Basic math
#include <cmath>

#include "dataStructures.h"

using namespace std;

int Fact(int x, int fim);

unsigned long Binom(int a, int b);

double Det(STLMat &A, int n);

int IntLog2(int x);

template <class InputIterator>
int SumI(InputIterator dummy, InputIterator stop)
{
	int out = 0;
	for( ; dummy!=stop;++dummy)
		out+= *dummy;
	return out;
};

template<typename T1, typename T2>
T2 min(T1 X, T2 Y)
{
    return X > Y ? Y : X;
};

template<typename T1>
double cotan(T1 arg)
{
	return 1./tan(arg);
};

template<typename T1, typename T2>
T2 maxi(T1 X, T2 Y)
{
    return X < Y ? Y : X;
};

template <class InputIterator, class OutputIterator, typename T1>
void ScalarMultC(InputIterator start, InputIterator stop,
				    OutputIterator result, T1 a)
{
	for( ; start!=stop ;++start)
	{
		(*result) = a*(*start);
		++result;
	}
};

template <class InputIterator, typename T1>
void ScalarMultC(InputIterator start, InputIterator stop, T1 a)
{
	for( ; start!=stop ;++start)
	{
		(*start) = a*(*start);
	}
};

template <class InputIterator1, class InputIterator2, class OutputIterator>
void VectorMultC(InputIterator1 start1, InputIterator2 start2,
			 InputIterator1 stop, OutputIterator result)
{
	for( ; start1!=stop ;++start1)
	{
		(*result) = (*start2)*(*start1);
		++start2;
		++result;
	}
};

template <class InputIterator1, class InputIterator2, typename T>
T ScalarProdC(InputIterator1 start1, InputIterator2 start2, InputIterator1 stop)
{
	T out = 0;
	for( ; start1!=stop ;++start1)
	{
		out += (*start2)*(*start1);
		++start2;
	}

	return out;
};

double NormC(vector<double>& vector_in);

double RenormC(vector<double>& vector_in);


template <class InputIterator, class OutputIterator>
void VectorAddC(InputIterator start, InputIterator stop,
				  OutputIterator result)
{
	for( ; start!=stop ;++start)
	{
		(*result) = (*start) + (*result);
		++result;
	}
};

template <typename T>
int Sign(T n)
{
	if(n >= 0)
		return 1;
	else if(n < 0)
		return -1;
	else
		return 0;
};

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
};

template<typename T>
bool IsEqual(const T& vectorA,const T& vectorB)
{
	typename T::const_iterator A_it 	= vectorA.begin();
	typename T::const_iterator A_end 	= vectorA.end();

	typename T::const_iterator B_it	= vectorB.begin();

	return equal(A_it,A_end,B_it);
};

template<typename T1, typename T2>
bool CompareFloats(T1 a,T2 b,double diff)
{
	return (abs(a - b) < diff );
};

template<typename T>
double CalculateIPR(T VecInput[],int size)
{
	double sumSqr = 0.;
	double sumQuart = 0.;

	double out = 0.;
	for(int iii = 0; iii < size; ++iii)
	{
		sumSqr += pow(VecInput[iii],2);
		sumQuart += pow(VecInput[iii],4);
	}
	out = sumQuart / pow(sumSqr,2);
	return out;
};

bool IsEven(int v);

inline int idxConv(int NNN, int iii, int jjj)
{
	return NNN*iii + jjj;
}

inline int neighIdxConv(int level, int jjj)
{
	return 3*level*(level-1) + jjj;
}
#endif
