#ifndef DATASTRUC_H_
#define DATASTRUC_H_

// // STL containers and algorithms
#include <vector>
#include "boost/random.hpp"

using namespace std;

// 1D partitions
typedef vector<char> part1D;
typedef vector<part1D> part1DList;

// XXZ chains
typedef vector<char> XXZChain;
typedef vector<part1D> XXZChainList;

// 2D partitions
typedef vector<char> part2D;
typedef vector<part2D> part2DList;

// typedef map<int,vector<char> > partMap;

typedef vector<double> Dvector;
typedef vector<Dvector > DvectorList;
typedef vector<Dvector > STLMat;

// Index lists 
typedef vector<int> idxList;
typedef vector<idxList> idxTable;

// Index lists 
typedef vector<int> Ivector;
typedef vector<Ivector> IvectorList;

// Descent paths
// -> Structure for the descent paths
struct Descent
{
	idxList	descIndex;	// The descent path
	int		descNumber;	// Number of descents
};

// -> Descent paths list
typedef vector<Descent> DescentList;

// Eigenlists
// -> Structure for the eigenpairs
// struct EigPair
// {
// 	PetscScalar		value;
// 	Vec			 	vec;
// };
// 
// // -> Eigenpairs list
// typedef vector<EigPair> EigList;

struct range
{
	double rbegin;
	double rend;
};

template<typename T>
void DeallocateVector(vector<T> &vec)
{
	vec.clear();
	vector<T>().swap(vec);
};
/*
template<typename T>
void DeallocateMatrix(vector<T> &vec)
{
	for(unsigned int iii = 0; iii < vec.size(); ++iii)
		DeallocateVector(vec[iii]);
	
	vec.clear();
	vector<T>().swap(vec);
};*/

struct coord
{
	int	x;
	int	y;
};

inline coord operator-(const coord &c1, const coord &c2)
{
	coord dummy = {c1.x - c2.x,c1.y - c2.y};
	return dummy;
};

inline coord operator+(const coord &c1, const coord &c2)
{
	coord dummy = {c1.x + c2.x,c1.y + c2.y};
	return dummy;
};

inline bool operator<(const coord &c1, const coord &c2)
{
	if(c1.y < c2.y || (c1.y==c2.y&&c1.x<c2.x))
	{
		return true;
	}
	else
	{
		return false;
	}
};

typedef vector<char> spinLayer;
typedef vector<spinLayer> spinStack;

typedef unsigned long ulong;

#endif
