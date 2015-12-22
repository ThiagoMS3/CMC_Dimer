#ifndef SYSCONF_H
#define SYSCONF_H

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <complex>

#include "dataStructures.h"
#include "customMath.h"
#include "MCInputParams.h"
#include "MCParameters.h"

#include "boost/random.hpp"

// Definitions
using namespace std;
using namespace boost;

class SysConf
{
private:
	// *******************************************************************
	// ---> Method used to print / read the class elements.
	//	   Uses "boost::serialization"
	//      !!! This library must be compiled beforehand !!!
	// *******************************************************************

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
//		ar & 	neighboursTable;
		ar & 	spinConf;
//		ar & 	flippableSpinLists;
//		ar & 	flippableSpinPositions;
		ar & 	CoreStart;
		ar & 	CoreEnd;
		ar & 	CoreSize;
		ar & 	initCondType;
		ar & 	nx;
		ar & 	ny;
		ar & 	p;
		ar & 	L;
		ar & 	LTot;
		ar & 	LBorder;
		ar & 	N;
		ar & 	NDimers;
		ar & 	NbOfNeights;
		ar & 	TotalOldNbF;
		ar & 	TotalNewNbF;
		ar & 	SimType;
		ar &	RNGSeed;
		ar &	BorderType;
		ar & 	ConfigFile;
		ar & 	Kz;
		ar & 	Kt;
		ar & 	Vt;

		ar & alpha0;
		ar & alpha3;

		ar & shift0;
		ar & shift3;

		ar & measurementsDone;
	}

	// *******************************************************************
	// Private variables
	// *******************************************************************

	// RNG
	lagged_fibonacci607 m_rng;	// RNG engine ---> has to be sent by other methods
	int 				RNGSeed;	// Seed

	// Temporary store the changes of the flippable spins table
	vector<int>     		nbToRemove; 	// Nb of spins which must be changed (removed/added to the
	vector<int>     		nbToAdd;		//	flippableSpin lists)
	vector<vector<int> >	toRemove;
	vector<vector<int> >	toAdd;
	vector<int>			dummyToAdd;
	vector<int>			dummyToRemove;

	int maxExp;
	int minExp;

	// For the diag = N0 case, we only need /Delta N0
	vector<int> 			DeltaN0;		//
	int					Total_DeltaN0;

	// Acceptance variables
	bool 		IsAccepted;
	double 		accept;
	ulong 		totalAccept;
	double		baseProb;

	// Energy
	double 		energy;
	double 		DE;

	// MC dummy variables ---> won't be sent by 'serialize', put here only to avoid memory re-allocations
	int 			absIndex;				// Absolute index of the spin inside 'flippableSpinLists'
	int 			spinToFlip;			// Position of the spin over the whole stack
	int 			nnnCluster;			// Spin's layer
	int			kkkCluster;			// Spin's position inside the layer
	char 		firstSpin;

	int 	clusterStart, clusterEnd, clusterSize;
	int	LowerBorderIsF, UpperBorderIsF;
		// Are the borders of the cluster flippable?
		// ____BoundIsF	= 0 -> false
		//	  			= 1 -> true

	// Vectors used for the staggered magnetism
	vector<int> 	spinChain;
	vector<char> 	equivPart;
	vector<int>		horizontalDimerPos;

	// Variables used to convert from a 2D index system to a 1D table
	vector<coord> offset;

	// Core geometry -> Internal hexagon
	vector<int> CoreStart;
	vector<int> CoreEnd;
	int CoreSize;

	vector<int> rowStart;
	vector<int> rowEnd;
	vector<int> rowSize;

	// Vectors keeping the calculated exponentials
	vector<double> expZ;

	vector<char> neighsSpin;
	int TypeOfFlipTest;

public:
	// *******************************************************************
	// Public variables
	// *******************************************************************

	// ---> Spin configuration
	vector<int>	neighboursTable;
			// Table keeping the neighbours
			// neighboursTable[6*iii+jjj]
			//				= jjj'th neighbour of the iii'th spin

	vector<int>	secondNeighboursTable;
	vector<vector<int> > nthNeighboursTable;
	int nbOfNeighLevels;

	vector<double> m_x;
	vector<double> m_y;

	vector<int> SqrDistancesTable;
	vector<int> SqrDistancesKey;
	vector<int> SqrDistancesPos;

	int maximumSqrDist;
	int nbOfDists;

	vector<int>	weightTable;
			// Table keeping the weight of each neighbouring relation
			// weightTable[6*iii+jjj]
			//				= jjj'th neighbour weight of the iii'th spin
	vector<char>	spinConf;
			// Table keeping the spin configurations
			// spinConf[N*iii+nnn]
			//				= spin at the layer nnn, position iii

	/*	We can choose the flippable spins in three ways :
		1) vector<vector<char> >	flippableSpinTable
			Choose an random spin @ position [iii][jjj] and see if it's flippable.
			We have to keep an auxiliary table saying if the spin at the position [iii][jjj] is flippable.

				Advantages   : - no need to reallocate memory
							- memory access is O[1]
			Disavantages : - the rejection rate will be high if the spin density is small
							- will use a lot of memory

				Better suited for high N_f/N*L values.

		2) vector<vector<int> >	flippableSpinLists
				Keep a flippable spins list and choose one spin randomly

				Advantages   : - uses less memory than 1)
							- we can manipulate the memory usage with resize()/reserve() to adapt
							the memory usage
							- no rejections
				Disavantages : - needs to insert/remove elements, which are O[N]
							- we have to do searches inside the vectors when we update then
							- removal of elements

				Better suited for low N_f/N*L values
	*/
	vector<int>	flippableSpinLists;
			// Table saving the spins that can be flipped
			// flippableSpinLists[iii] = (unordered) array of the flippable spin
			//		in the layer 'iii'
	vector<int>	flippableSpinPositions;
			// Table saving the position of an spin inside 'flippableSpinLists'
			// flippableSpinPositions[iii] = position of spin 'iii' inside
			//		'flippableSpinLists'

	// ---> Geometry
	int		nx;
	int		ny;
	int		p;

	int 		L;
	int 		LTot;
	int		LBorder;

	int		N;
	int 		NDimers;
	int		NbOfNeights;

	int 		TotalOldNbF;		// Number of flippable sites of the actual configuration
	int 		TotalNewNbF;		// Number of flippable sites of the proposed configuration

	vector<int>	SublatticeSize;		// Number of elements in the sublattices
	vector<int>	SublatticePositions;
	vector<int>	SublatticeA;		// Sublattices
	vector<int>	SublatticeB;		//
	vector<int>	SublatticeC;		//

	vector<double>	dummySubs;

	vector<int>	DimerDensityProfile;		// Number of elements in the sublattices
	vector<int>	N3_Map;		// Number of elements in the sublattices

	vector<complex<double> > Jindex;

	// ---> MC parameters
	double		Kz;
	double		Kt;

	double		q;
	double		Vt;

	double alpha0;
	double alpha3;

	double shift0;
	double shift3;

	// ---> Simulation type
	int 			BorderType;	// Boundary conditions
	string			SimType;		// Simulation type (Moessner, Partition)
	string			ConfigFile;		// Neighbour table file
	int 			initCondType;	// Inital condition (Star3, Fixed borders, Random)

	// ---> Simulation parameters
	int				measurementsDone;

	// *******************************************************************
	// Constructors and associated methods
	// *******************************************************************

	// Default constructor
	SysConf() {SetRNG(42);}

	// Set "Seed"
	SysConf(double Seed) {SetRNG(Seed);}

	// Construct from parameters
	SysConf(const input_params& input,const string& TypeIn);

	// Copy (not really used ... ImportObject takes care of it)
	SysConf(const SysConf& original);

	// Methods used in the constructors
	void SetParam(const input_params& input,const string& TypeIn);
	void Copy(const SysConf& original);

	void SetNeighbours();
	void SetFlipTables();
	void SetSecondNeightbours();
//	void Move(int& pos, int direction);
//	void SetNthNeightbours();

	// *******************************************************************
	// Initial state / SimType == "Moessner"
	// *******************************************************************
	void SetNeighboursMoessner();
	void SetInitialMoessner();

	// *******************************************************************
	// Initial state / SimType == "Manual"
	// *******************************************************************
	void SetConfigFileManual(string& configName);

	void SetNeighboursManual();
	void SetInitialManual(string& configName);

	void SetAntiWeights();
	void ReadAntiWeights();

	// *******************************************************************
	// Initial state / SimType == "Partition"
	// *******************************************************************
	void SetHexOffset();
	void SetNeighboursPart();

	void SetInitialPart();			// Random partition
    void SetInitialPartEmpty();			// Random partition
	void SetInitialPartAntiFerro();		// V/t -> -\infty partition
	void SetInitialPartStar3();		// Star3 conf - used for open boundaries
	void SetInitialPartGeneral();

	// *******************************************************************
	// RNG setup
	// *******************************************************************
	void SetRNG(double Seed);
	void ExportRNG(const char outfile[]);
	void ImportRNG(const char infile[]);
	void PrintRNG_DEBUG(int nnn);

	// *******************************************************************
	// MC methods - Common
	// *******************************************************************
	bool TestFlippable(int spin, int layer);	// Test if a spin can be flipped
	void ClusterUpdate();					// If the cluster update was accepted, update it

	bool GetRandomAccept(double inAccept);		// Random acceptance
	int  GetAbsIndex(int InTotalOldNbF);		// Random absolute index

	void TestUpdate(int layer);	// Find out the changes of the spin configuration if a spin is flipped

	void RaiseMeasures(int UpdateInterval);

	// *******************************************************************
	// MC methods - Diagonal = N3
	// *******************************************************************
	void MC_Prepare(MCParameters& inParams);	// >>>> Prepare MC simulation
	void CalculateE();			// Calculate the energy

	void ClusterIter();			// Do a MC run :
	void BuildCluster();		// > Build the cluster

	void GrowCluster(int& workPoint, int& borderStatus, int step, int border, int borderRedirect);
	void GetAcceptanceCluster();	// > Get its acceptance

	// *******************************************************************
	// MC methods - Diagonal = N0
	// *******************************************************************
	void MC_Prepare_Zero(MCParameters& inParams);	// >>>> Prepare MC simulation
	void CalculateE_Zero();			// Calculate the energy

	void ClusterIter_Zero();			// Do a MC run :
	void BuildCluster_Zero();		// > Build the cluster

	void GrowCluster_Zero(int& workPoint, int& borderStatus, int step, int border, int borderRedirect);
	void GetAcceptanceCluster_Zero();	// > Get its acceptance

	void Update_DeltaN0(int layer);		// Find out \Delta N0 for a given layer

	// *******************************************************************
	// MC methods - Diagonal = a_0 * N_0 + a_3 * N_3
	// *******************************************************************
	void MC_Prepare_Mixed(MCParameters& inParams);	// >>>> Prepare MC simulation
	void CalculateE_Mixed();			// Calculate the energy

	void ClusterIter_Mixed();			// Do a MC run :
	void BuildCluster_Mixed();		// > Build the cluster

	void GrowCluster_Mixed(int& workPoint, int& borderStatus, int step, int border, int borderRedirect);
	void GetAcceptanceCluster_Mixed();	// > Get its acceptance

	// *******************************************************************
	// Order Parameters / getters
	// *******************************************************************

	// ---> Basic MC
	double GetEnergy();						// Energy
	double GetClusterSize();				// Cluster size
	double GetMagnetization();				// RMS magnetization (2D)
	double GetCoreMagnetization();			// RMS magnetization of the core hexagon
	ulong  GetTotalAccept();				// Percentage of accepted updates

	// ---> 1D XXZ spin chain
	void   BuildEquivalentSpinChain();		// Build the equivalent XXZ spin chain
	double CalculateStagMag();				// Calculate the staggered mag. for the spin chain
	double CalculateChainEnergy();			// Calculate the staggered mag. for the spin chain
	void   CalculateChainCorrelations(vector<double>& CorrSz, vector<double>& CorrStag);

	// ---> Mean partition equivalence
	void   CalculateMeanPartition();
	void   GetHeights(int nbOfParts, int lll);
	void   SetParts(int nbOfParts, int nxInitCol, int nyInitRow);

    void   GetMeanPart(vector<double> &output);

	// ---> Nf
	void   GetNf(vector<double> &Ncount);	// Total number of sites with {0, 1, 2, 3} dimers
	int    GetN0();							// Total number of sites with 0 dimers
	int    HasZeroDimers(int spin, int layer);
	int    GetLocalField(int lll, int nnn);

	void   SetSublatticeManual();
	void   SetSublatticeMoessner();
//	void   GetSublatticePeriodic(vector<double>& output, vector<double> & NSiteCount);
	void   GetSublattice(vector<double>& output, vector<double> & NSiteCount);

	void   GetSiteNf(	vector<double> & NSiteCount,
								vector<double> & Ncount,
								vector<double> & Ndimer,

								vector<double> & MeanSublattice,

								vector<complex<double> >& complexPhaseVector,
								vector<double>& realPhaseVector,
								complex<double>& complexPhase,
								double& realPhase);

	void   GetSiteNf(	vector<double> & NSiteCount,
						vector<double> & Ncount,
						vector<double> & Ndimer,

						vector<double> & MeanSublattice);

	void   GetSiteNf(	vector<double> & NSiteCount,
						vector<double> & Ncount,
						vector<double> & Ndimer);

    double GetPlaquetteTest();

	void   GetDimerDimerCorrelation(vector<double> & corr);
	void   GetSzSzCorrelation(vector<double> & corr);
	void   GetSxSxCorrelation(vector<double>& corr);

	// >>>> Site - site
	void   GetN3Correlation(double &corr, vector<double> &Star3Network);
	void   GetN3N3SpacialCorrelation(vector<double>& N3SpatialCorrelation, vector<double>& N3LocalDensities);

	// ---> energies
	double GetQEnergy_N3();
	double GetQEnergy_N0();
	double GetQEnergy_Mixed();

	double GetQNewEnergy_N3();
	double GetQNewEnergy_N0();
	double GetQNewEnergy_Mixed();

	double GetPotEnergy_N3();
	double GetPotEnergy_N0();
	double GetPotEnergy_Mixed();

	double GetKinEnergy();
	double GetNewKinEnergy();

	void   GetKinEnergyDensity(vector<double> & KinOut);
	void   GetPotV3EnergyDensity(vector<double> & PotOut);
	void   GetPotV0EnergyDensity(vector<double> & PotOut);

	void   GetAllQEnergy_Mixed(vector<double> & dummyVec);
	void   GetAllNewQEnergy_Mixed(vector<double> & dummyVec);

	// *******************************************************************
	// Conf randomizers
	// *******************************************************************
	void SetRandomColFlips();
	void RandomSpinFlips(int flips);		// Single spins
	void RandomColFlip();		// Whole columns

	// *******************************************************************
	// Output/debug
	// *******************************************************************
	void PrintConfigStatus();
	void PrintDebug(ostream& out);

	// *******************************************************************
	// Correlation operators
	// *******************************************************************
	void SetPositions();
	int MinimalSqrDistance(int pointA, int pointB);
	void SetDistances();
};

#endif
