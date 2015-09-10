#ifndef MCINPUTPARAMS_H_
#define MCINPUTPARAMS_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "dataStructures.h"
#include "customMath.h"

// Definitions
using namespace std;
using namespace boost;

// --- I/O funtions ---
//---------------------
template<typename T>
void Jump(T& filestream, int numberOfLines)
{
	string dummy;
	for(int iii = 0; iii < numberOfLines; ++iii)
		getline(filestream,dummy);
};

template<typename T>
string ToString(const T &arg)
{
	stringstream out;
	out << arg;
	return out.str();
};

// --- Input struct ---
// --------------------
class input_params
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & nx;
		ar & ny;
		ar & p;		// Dimensions
		ar & NbOfNeights;		// Number of neighbours of each site

		// Hamiltonian parameters
		ar & VtType;
		ar & VtN;
		ar & VtValues;
		ar & VtStart;
		ar & VtEnd;

		ar & alpha0;
		ar & alpha3;

		ar & shift0;
		ar & shift3;

		// Classical MC parameters
		ar & N; 			// Number of stacks
		ar & Dbeta;			// Temperature step
		ar & NbMaxUpdate;		// Number of updates before calculating an observable

		ar & NbMeasures;		// Number of measures that we want;
		ar & MeasuresToTake;

		ar & AutocorrLength;
		ar & minBinNb;

		// -> Computer config
		ar & compType;
		ar & numProc;

		// -> Output files
		ar & outFolderName;
		ar & execName;
		ar & runType;

		ar & initCondType;
		ar & RNGSeed;
	}

public:
	// -> Partition / lattice geometry
	int 		nx, ny, p;		// Dimensions
	int 		NbOfNeights;		// Number of neighbours of each site

	// Hamiltonian parameters
	int		VtType;
	int  		VtN;
	vector<double>	VtValues;
	double 		VtStart;
	double 		VtEnd;

	double alpha0;
	double alpha3;

	double shift0;
	double shift3;

	// Classical MC parameters
	int 		N; 			// Number of stacks
	double		Dbeta;			// Temperature step
	ulong		NbMaxUpdate;		// Number of updates before calculating an observable

	int		NbMeasures;		// Number of measures that we want;
	int		MeasuresToTake;

	int		AutocorrLength;
	int		minBinNb;

	// -> Computer config
	int 		compType;
	int		numProc;
	int 		numNodes;

	// -> Output files
	string		outFolderName;
	string		execName;
	string		runType;

	int		initCondType;
	int 		RNGSeed;

	int 		BorderType;

	input_params() {}

	void SetParamMoessner(ifstream& inputF)
	{
		runType = "Moessner";
		p = -1;
		initCondType = 0;

		Jump(inputF,4);

		double dVt;

		// Get partition problem geometry
		inputF >> nx; Jump(inputF,1);
		inputF >> ny; Jump(inputF,1);
		inputF >> NbOfNeights; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtType;Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtN; Jump(inputF,1);
		VtValues.resize(VtN);

		Jump(inputF,1);
		if(VtType==0)
		{
			inputF >> VtStart >> VtEnd;Jump(inputF,1);
			dVt = (VtEnd - VtStart)/(VtN - 1);
			for(int iii = 0; iii < VtN; ++iii)
			{
				VtValues[iii] = VtStart + iii * dVt;
			}
		}
		else
		{
			for(int iii = 0; iii < VtN; ++iii)
			{
				inputF >> VtValues[iii];
			}
			VtStart = VtValues[0];
			VtEnd = VtValues[VtN - 1];
			Jump(inputF,1);
		}
		Jump(inputF,2);

		inputF >> alpha0; Jump(inputF,1);
		inputF >> alpha3; Jump(inputF,1);

		Jump(inputF,2);

		inputF >> shift0; Jump(inputF,1);
		inputF >> shift3; Jump(inputF,1);

		Jump(inputF,2);

		// Get MC parameters
		inputF >> N; Jump(inputF,1);
		inputF >> Dbeta; Jump(inputF,1);
		inputF >> NbMaxUpdate; Jump(inputF,2);

		inputF >> NbMeasures; Jump(inputF,1);
		inputF >> MeasuresToTake; Jump(inputF,1);
		inputF >> AutocorrLength; Jump(inputF,1);
		inputF >> minBinNb; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> compType; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> numProc; Jump(inputF,1);
		inputF >> BorderType; Jump(inputF,1);
		Jump(inputF,3);
		inputF >> RNGSeed; Jump(inputF,1);
		Jump(inputF,2);

		getline(inputF,outFolderName);
		Jump(inputF,2);

		getline(inputF,execName);
	}

	void SetParamManual(ifstream& inputF)
	{
		runType = "Manual";
		p = -1;
		initCondType = 0;

		Jump(inputF,4);

		double dVt;

		// Get partition problem geometry
		inputF >> nx; Jump(inputF,1);
		inputF >> ny; Jump(inputF,1);
		inputF >> NbOfNeights; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtType;Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtN; Jump(inputF,1);
		VtValues.resize(VtN);

		Jump(inputF,1);
		if(VtType==0)
		{
			inputF >> VtStart >> VtEnd;Jump(inputF,1);
			dVt = (VtEnd - VtStart)/(VtN - 1);
			for(int iii = 0; iii < VtN; ++iii)
			{
				VtValues[iii] = VtStart + iii * dVt;
			}
		}
		else
		{
			for(int iii = 0; iii < VtN; ++iii)
			{
				inputF >> VtValues[iii];
			}
			Jump(inputF,1);
		}
		Jump(inputF,2);

		inputF >> alpha0; Jump(inputF,1);
		inputF >> alpha3; Jump(inputF,1);

		Jump(inputF,2);

		inputF >> shift0; Jump(inputF,1);
		inputF >> shift3; Jump(inputF,1);

		Jump(inputF,2);

		// Get MC parameters
		inputF >> N; Jump(inputF,1);
		inputF >> Dbeta; Jump(inputF,1);
		inputF >> NbMaxUpdate; Jump(inputF,2);

		inputF >> NbMeasures; Jump(inputF,1);
		inputF >> MeasuresToTake; Jump(inputF,1);
		inputF >> AutocorrLength; Jump(inputF,1);
		inputF >> minBinNb; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> compType; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> numProc; Jump(inputF,1);
		inputF >> BorderType; Jump(inputF,1);
		Jump(inputF,3);
		inputF >> RNGSeed; Jump(inputF,1);
		Jump(inputF,2);

		getline(inputF,outFolderName);
		Jump(inputF,2);

		getline(inputF,execName);
		PrintDebug(cout);
		cout.flush();
	}

	void SetParamPart(ifstream& inputF)
	{
		runType = "Part";
		Jump(inputF,4);

		double dVt;

		// Get partition problem geometry
		inputF >> nx; Jump(inputF,1);
		inputF >> ny; Jump(inputF,1);
		inputF >> p; Jump(inputF,1);
		inputF >> NbOfNeights; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtType;Jump(inputF,1);
		Jump(inputF,2);

		inputF >> VtN; Jump(inputF,1);
		VtValues.resize(VtN);

		Jump(inputF,1);
		if(VtType==0)
		{
			inputF >> VtStart >> VtEnd;Jump(inputF,1);
			dVt = (VtEnd - VtStart)/(VtN - 1);
			for(int iii = 0; iii < VtN; ++iii)
			{
				VtValues[iii] = VtStart + iii * dVt;
			}
		}
		else
		{
			for(int iii = 0; iii < VtN; ++iii)
			{
				inputF >> VtValues[iii];
			}
			Jump(inputF,1);
		}
		Jump(inputF,2);

		inputF >> alpha0; Jump(inputF,1);
		inputF >> alpha3; Jump(inputF,1);

		Jump(inputF,2);

		inputF >> shift0; Jump(inputF,1);
		inputF >> shift3; Jump(inputF,1);

		Jump(inputF,2);

		// Get MC parameters
		inputF >> N; Jump(inputF,1);
		inputF >> Dbeta; Jump(inputF,1);
		inputF >> NbMaxUpdate; Jump(inputF,2);

		inputF >> NbMeasures; Jump(inputF,1);
		inputF >> MeasuresToTake; Jump(inputF,1);
		inputF >> AutocorrLength; Jump(inputF,1);
		inputF >> minBinNb; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> compType; Jump(inputF,1);
		Jump(inputF,2);

		inputF >> numProc; Jump(inputF,1);
		inputF >> initCondType; Jump(inputF,1);
		Jump(inputF,3);
		inputF >> BorderType; Jump(inputF,1);
		Jump(inputF,3);
		inputF >> RNGSeed; Jump(inputF,1);
		Jump(inputF,2);

		getline(inputF,outFolderName);
		Jump(inputF,2);

		getline(inputF,execName);
	}

	void PrintDebug(ostream& out)
	{
		out << " ### DEBUG Print input_params ###" << endl;
		out << " type                   = " << runType << endl;
		out << " nx                     = " << nx << endl;
		out << " ny                     = " << ny << endl;
		out << " p                      = " << p << endl;
		out << " NbOfNeights            = " << NbOfNeights << endl << endl;

		out << " VtType                 = " << VtType << endl;
		out << " VtN                    = " << VtN << endl;
		out << " Vt's                   = " << VtN << endl;

		for(int iii = 0; iii < VtN; ++iii)
		{
			out << " " << VtValues[iii];
		}
		out << endl << endl;

		// Get MC parameters
		out << " N                      = " << N << endl;
		out << " Dbeta                  = " << Dbeta << endl;
		out << " NbMaxUpdate            = " << NbMaxUpdate << endl << endl;

		out << " NbMeasures             = " << NbMeasures << endl;
		out << " MeasuresToTake         = " << MeasuresToTake << endl;
		out << " AutocorrLength         = " << AutocorrLength << endl << endl;

		out << " CompType               = " << compType << endl;
		out << " numProc                = " << numProc<< endl;

		out << " outFolderName          = " << outFolderName << endl;
		out << " ExecName               = " << execName << endl;
	}
};

#endif
