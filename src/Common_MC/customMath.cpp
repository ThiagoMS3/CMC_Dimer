#include "customMath.h"

/*
 *    Simple mathematical functions - definitions
 * 
 */

int IntLog2(int x)
{	
	double dummy = log(x)/log(2);
	return (int)dummy;
};

// Factorial
int Fact(int x, int fim=1)
{
	if (x == fim)
		return 1;
		return (x * Fact(x-1,fim));
}

// Binomial
unsigned long Binom(int a, int b)
{
	if (a < b)				return 0;
	else if	(b==0 || a==b)		return 1;
	else if 	(a>=0 && b>=0)		return Binom(a-1,b) + Binom(a-1,b-1);
	return 0;
}

// Calculate the determinant "out" of the matrix "A"
double Det(STLMat &A, int n)
{
	// Parameters initialization
	if(n == 1)
		return A[0][0];
	
	double out = 0.;
	
	STLMat 	minorMat(n-1,Dvector(n-1));
	
	for(int iii = 0; iii < n; ++iii)
	{
		for(int aux_i = 0; aux_i < iii; ++aux_i)
			for(int aux_j = 0; aux_j < n-1; ++aux_j)
				minorMat[aux_i][aux_j] = A[aux_i][aux_j+1];
		for(int aux_i = iii; aux_i < n-1; ++aux_i)
			for(int aux_j = 0; aux_j < n-1; ++aux_j)
				minorMat[aux_i][aux_j] = A[aux_i+1][aux_j+1];
		out += pow(-1.,iii) * A[iii][0]*Det(minorMat,n-1);
	}
	
 	DeallocateVector(minorMat);
	
	return out;
};

double NormC(vector<double>& vector_in)
{
	double out = 0.;
	
	vector<double>::iterator vec_it, vec_end;
	
	vec_it = vector_in.begin();
	vec_end = vector_in.end();
	
	for( ; vec_it!=vec_end ;++vec_it)
	{
		out += (*vec_it)*(*vec_it);
	}
	
	return out;
};

double RenormC(vector<double>& vector_in)
{
	double norm = 0.;
	
	vector<double>::iterator vec_it, vec_end;
	
	vec_it = vector_in.begin();
	vec_end = vector_in.end();
	
	for( ; vec_it!=vec_end ;++vec_it)
	{
		norm += (*vec_it)*(*vec_it);
	}
	
	norm = sqrt(norm);
	vec_it = vector_in.begin();
	
	for( ; vec_it!=vec_end ;++vec_it)
	{
		(*vec_it) = (*vec_it)/norm;
	}

	return norm;
};

bool IsEven(int v) // word value to compute the parity of
{ 
	return !(v%2);
}