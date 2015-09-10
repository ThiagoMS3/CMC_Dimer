#ifndef MCISING_H_
#define MCISING_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <iostream>
#include <fstream>
#include <map>

#include "dataStructures.h"
#include "customMath.h"

// #include "SysConf.h"
#include "OrderParam.h"
#include "MCParameters.h"
#include "MCInputParams.h"

#include "boost/random.hpp"
#include <time.h>

// Definitions
using namespace std;
using namespace boost;

extern lagged_fibonacci607 rng;

// *** Inline functions and templates ***
// **************************************

// --- Observables / energy ---
// ----------------------------
inline double MeanMagnetization(vector<char> const &spinConf,int L, int N)
{
	double Mag = 0;

	for(int iii = 0; iii < N; ++iii)
	{
		for(int jjj = 0; jjj < L; ++jjj)
		{
			Mag += (double)spinConf[idxConv(N,jjj,iii)];
		}
	}
	return pow(Mag/(L*N),2);
};

inline double CalculateE(double Kz,
				double Kt,
				double Vt,
				
				int TotalOldNbF,
				int N,
				int L,
				vector<char> const &spinConf)
{
	double Energy = Kz * TotalOldNbF * Vt;
	
	for(int nnn = 0; nnn < N-1; ++nnn)
		for(int iii = 0; iii < L; ++iii)
		{
			Energy -= Kt*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn)+1];
		}
		
	for(int iii = 0; iii < L; ++iii)
	{
		Energy -= Kt*spinConf[idxConv(N,iii,0) + N-1]*spinConf[idxConv(N,iii,0)];
	}
	
	return Energy;
};

// --- Get random spin / build cluster ---
// ---------------------------------------
inline int GetRandomSpin(
			    int OldNbF
			   )
{
	random::uniform_int_distribution<> dist(0, OldNbF - 1);
	return dist(rng);
};

inline bool GetRandomAccept(
			    double accept
			   )
{
	bernoulli_distribution<> bernoulliMC(accept);
	return bernoulliMC(rng);
};

inline int GetAbsIndex(int TotalOldNbF)
{
	random::uniform_int_distribution<> dist(0, TotalOldNbF - 1);
	return dist(rng);
};

inline double GetAcceptance(double Kz,
				double Kt,
				double Vt,
				
				int TotalOldNbF,
				int TotalNewNbF,
				
				int SpinToFlip,
				int nnn,
				int N,
				spinStack& spinConf,
				double& DE
			    )
{
	// > Calculate the energy
	int SpinAbove, SpinBelow;
	if(nnn == 0)
	{
		SpinBelow = spinConf[SpinToFlip][N-1];
	}
	else
	{
		SpinBelow = spinConf[SpinToFlip][nnn-1];
	}
	
	if(nnn == N-1)
	{
		SpinAbove = spinConf[SpinToFlip][0];
	}
	else
	{
		SpinAbove = spinConf[SpinToFlip][nnn+1];
	}
	
	DE = Kz*Vt*(TotalNewNbF - TotalOldNbF) + 
				2*Kt*(int)spinConf[SpinToFlip][nnn]*(SpinAbove + SpinBelow);
				
	// Calculate the acceptance
	double accept = min(1,exp(-DE)*(double)TotalOldNbF/TotalNewNbF);
	
	return accept;
};

inline double GetAcceptanceCluster(double Kz,
				double Kt,
				double q,
				double Vt,
				
				int TotalOldNbF,
				int TotalNewNbF,
			        int nbToChangeAtFirst,
				int nbToChangeAtUpper,
				int nbToChangeAtLower,
				
				int UpperSpinIsFlippable,
				int LowerSpinIsFlippable,
				
				int firstSpin,
				int upperSpin,
				int lowerSpin,
				double& DE
			    )
{
	// > Calculate the energy	
	double accept = 0;
	
	// VERIFIED
	DE = Kz*Vt*(TotalNewNbF - TotalOldNbF) + 
				2*Kt*firstSpin*(upperSpin + lowerSpin);
					
	// Extra weight given by the cluster update
	double weight = exp(-Kz*Vt*(nbToChangeAtFirst)-2*Kt*firstSpin*(upperSpin + lowerSpin));
								// Spin flip - OK
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*nbToChangeAtUpper)),-UpperSpinIsFlippable*firstSpin*upperSpin);
								// Cluster boreder m
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*nbToChangeAtLower)),-LowerSpinIsFlippable*firstSpin*lowerSpin);
								// Cluster boreder m'
	
	// Calculate the acceptance
	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);
	
	return accept;
};

inline bool TestFlippable(vector<char>	const &spinConf,
			    int layer,
			    int N,
			    vector<int> const	&neighbours,
			    int ChosenNeigh,
			    int NbOfNeights
			   )
{
	int test = 0;
	int index = 0;
	for(int iii = 0; iii < NbOfNeights;++iii)
	{
		index = idxConv(NbOfNeights,ChosenNeigh,iii);
		test += spinConf[idxConv(N,neighbours[index],layer)];
	}
	
	return (test==0);
};

// --- I/O funtions ---
//---------------------

// *** Functions declarations ***
// ******************************

// --- Build spin configuration ---
// --------------------------------
// void SetNeighbours(int L,					// Lattice dimensions
// 			    int LTot,
// 			    int NbOfNeights,			//
// 			    int nx,				//
// 			    int ny,				//
// 			    int p,				//
// 			    vector<coord> &offset,
// 			    idxTable &neighboursTable		// Output
// 			   );

// void SetNeighboursMoessner(int nx,				// lines
// 			    int ny,				// cols
// 			    vector<int> &neighboursTable		// Output
// 			   );

// void SetInitialSpinConf(int N,
// 				    int L,
// 				    int LTot,
// 				    int NbOfNeights,
// 				    idxTable &neighboursTable,
// 				    
// 				    spinStack &spinConf,
// 				    idxList &flippableSpinLists,
// 				    idxList &flippableSpinPositions,
// 				    int &TotalOldNbF
// 				   );

// void SetInitialSpinConfMoessner(int N,
// 				    int L,
// 				int nx,
// 				int ny,
// 				    int NbOfNeights,
// 				    vector<int> const &neighboursTable,
// 				    
// 				    vector<char> &spinConf,
// 				    vector<int> &flippableSpinLists,
// 				    vector<int> &flippableSpinPositions,
// 				    int &TotalOldNbF
// 				   );

// --- Observables / energy ---
// ----------------------------
void CountNf(vector<char> const	&spinConf,
	     int 	 	N,
	     int		L,
	     
	     vector<int> const	&neighbours,
	     
	     int 	 	NbOfNeights,
	     vector<double> 	&Ncount
	    );

void CalculateStagMag(spinStack& spinConf,
		      idxTable& neighboursTable,
		      int LTot, int LBorder,int L, int N,
		      int nx, int ny, int p,
		      double& stagMag);

// --- Get random spin / build cluster ---
// ---------------------------------------
void SpinUpdate(spinStack& spinConf, int nnn, int spinToFlip, int TotalNewNbF, int& TotalOldNbF,
			 idxList& flippableSpinLists,idxList& flippableSpinPositions, idxList& toAdd, idxList& toRemove,
			 int nbToAdd, int nbToRemove);

void ClusterUpdate(vector<char>& spinConf,
		   int kkk,
		   int clusterSize,
		   int clusterStart,
		   int clusterEnd,
		   int N,

		   int TotalNewNbF,
		   int& TotalOldNbF,

		   idxList& flippableSpinLists,
		   idxList& flippableSpinPositions, 
		   
		   idxTable& toAdd, 
		   idxTable& toRemove,
		   idxList& nbToAdd,
		   idxList& nbToRemove
		  );

void TestUpdate(vector<char>& spinConf, int nnn, int N, int kkk, int L, int NbOfNeights, vector<int>& neighboursTable,
		idxList& flippableSpinPositions,
		vector<int>& toRemove, int& nbToRemove, vector<int>& toAdd, int& nbToAdd);

void BuildCluster(// Chosen spin
		  int nnn,
		  int kkk,
		  char firstSpin,
		  
		  // Energy parameters
		  double Kz,
		  double Vt,
		  double baseProb,
		  
		  // Geometry
		  int L,
		  int N, 
		  int NbOfNeights,
		  int TotalOldNbF,
		  vector<int> const& neighboursTable,
		  
		  // Spin configuration
		  idxList& flippableSpinLists,
		  idxList& flippableSpinPositions,
		  vector<char>& spinConf,

		  // Updated parameters
		  int& TotalNewNbF,
		  idxTable& toAdd,
		  idxList& nbToAdd,
		  idxTable& toRemove,
		  idxList& nbToRemove,
		  
		  int& clusterStart,
		  int& clusterEnd, 
		  int& clusterSize,
		  int& UpperBorderIsF,
		  int& LowerBorderIsF
		 );

// --- Bunching / binning ---
// --------------------------
// void BunchingErrors(vector<double>& data,
// 		    int order,
// 		    int NbMeasures,
// 		    vector<double>& errData
// 		   );
// 
// void Bunch(vector<double>& dataIn,
// 	   vector<double>& dataOut,
// 	   double& errData,
// 	   int NbOfPoints
// 	  );

// void DataBunch(double meas,vector<double> & dataBin, int & binSize, int minBin, int & binCounter, int & fillBin);

// --- Autocorrelations ---
// ------------------------
// void CalculateAutocorrelation(vector<double>& autocorr, double& integAutocorr, vector<double>& data, int dataSize, int acSize);

// --- I/O funtions ---
//---------------------
void MergeIntegACOutput(char filenameBase[], string& commonFile,
			int nodes, double VtStart, int VtN, double dVt, int data, int erase);

void MergeACOutput(const char filenameBase[], string& commonFile,
		   int nodes, int VtN, int AutocorrLength, int erase);

void MergeMCOutput(const char filenameBase[], string& commonFile,
		   int nodes,int VtN,int NValue, int erase);

void MergeBunchOutput(const char filenameBase[], string& commonFile,
		      int nodes, double VtStart, int VtN, double dVt, int OrderBunch, int erase);

// void rngLoop();

#endif