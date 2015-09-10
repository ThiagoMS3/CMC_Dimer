#ifndef ORDERPARAM_H_
#define ORDERPARAM_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#include "dataStructures.h"
#include "customMath.h"
#include "MCInputParams.h"


// Definitions
using namespace std;
using namespace boost;

// --- Big Huge Class for order parameters ---
// -------------------------------------------
class OrderParam
{
private:
// Method used to print / read the class elements. Uses "boost::serialization"
// !!! This library must be compiled beforehand !!!
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & 	m_Type;

		ar & 	m_data;
		ar & 	m_bin;
		ar & 	m_error;
		ar & 	m_errorBinSize;

		ar & 	m_minBinNb;
		ar & 	m_totalPoints;
		ar & 	m_binCount;
		ar & 	m_binSize;
		ar & 	m_fillBin;

		ar & 	m_AC;
		ar & 	m_IntegAC;
		ar & 	m_ACInterval;

		ar & 	m_mean;
		ar & 	m_corr;
		ar & 	m_corrSize;

		ar & 	m_nbOfMeas;

		ar & 	m_folder;
		ar & 	m_outMean;
		ar & 	m_outErr;
		ar & 	m_outAC;
		ar & 	m_outACInteg;

		ar & 	m_Vt;
	}

	// Type
	string m_Type;

	// Binning
	vector<double> m_data;		// Raw data used to do binning
	vector<double> m_bin;		// Vector used do for each binning step
	vector<double> m_error;		// Error values
	vector<double> m_errorBinSize;	// Bin sizes

	int m_minBinNb;		// Minimum number of bins
	int m_totalPoints;	// Total number of error points
	int m_binCount;		// number of bins filled
	int m_binSize;		// Final size of the bins
	int m_fillBin;		// Filled elements inside the actual bin

	// Autocorrelation
	vector<double> m_AC;
	double m_IntegAC;
	int m_ACInterval;

	// Mean values
	double m_mean;		// Mean value
	int m_nbOfMeas;		// Number of measurements

	// Correlation
	vector<double> m_corr;
	int m_corrSize;

	// Data files
	string m_folder;
	string m_outMean;
	string m_outErr;
	string m_outAC;
	string m_outACInteg;

	ofstream output;

	double m_Vt;

public:
	// ---> Constructors
	OrderParam() {}

	OrderParam(input_params& input);

	// ---> Setters
// 	void SetParam(vector<double>& data, int minBinNb, int binCount, int binSize, int fillBin, double mean, int nbOfMeas, int ACInterval);

	void SetFolder(const char folder[]);
	void SetType(const char dataType[]);
	void SetCorrLength(int size);
	int GetCorrLength();
	void CreateFiles(const char folder[],const char dataType[]);
	void CreateFiles();

	void SetVt(double Vt);

	// ---> Procedures
	// >>>> Add data to the mean/binning
	void AddData(double measurement);

	// >>>> Do binning
	void CalculateBinning(int& totalPoints);
	void CalculateAC();
	void CalculateError(int& totalPoints);

	void CalculateRMSBinning(int& totalPoints);
	void CalculateRMSAC();
	void CalculateRMSError(int& totalPoints);

	// >>>> Reset the mean value for a new series of measurements
	void RestartMean(int extraMeas);
	void RestartMean();
	void RestartCorr(int extraMeas);
	void RestartCorr();
	void RaiseMeasures(int extraMeas);

	void AddCorr(vector<double>& inp_corr);
	void GetCorr(vector<double>& outCorr);
	void GetCorr(vector<double>& outCorr,int length);

	// >>>> Get errors
	void GetError(vector<double>& error, vector<double>& sizes);

	// >>>> Get mean values (normal and RMS)
	double GetMean();
	double ConvertMean();
	int GetNbOfMeas();

	// ---> I/O
// 	void PrintMean();
// 	void PrintRMSMean();
	void PrintError();
	void PrintAC();
	void PrintBin();

	int GetNumberOfPoints()
	{
		return m_totalPoints;
	}

	int GetIACLength()
	{
		return m_ACInterval;
	}

	double GetIntegAC()
	{
		return m_IntegAC;
	}

// 	void PrintIntegAC();

	void CalculateRMSMean();
	void ResetRMSMean();

	double GetVt()
	{
		return m_Vt;
	}

	void PrintDebug(ostream& out);
};

#endif
