#ifndef TWODORDERPARAM_H_
#define TWODORDERPARAM_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "dataStructures.h"
#include "customMath.h"
#include "MCInputParams.h"
#include "SysConf.h"

// Definitions
using namespace std;
using namespace boost;

// --- Big Huge Class for order parameters ---
// -------------------------------------------
class TwoDOrderParam
{
private:
// Method used to print / read the class elements. Uses "boost::serialization"
// !!! This library must be compiled beforehand !!!
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & 	m_Type;
		ar & 	m_mean;
		ar & 	m_nbOfMeas;

		ar & 	m_row_number;
		ar & 	m_col_number;

		ar & 	m_folder;
		ar & 	m_outMean;

		ar & 	m_Vt;
		ar & 	m_parent;
//		ar &    m_x;
//		ar &    m_y;
	}

	// Type
	string m_Type;

	// Data
	vector<double> m_mean;		// Raw data
	int m_nbOfMeas;		// Number of measurements
	int m_row_number;
	int m_col_number;

	vector<double> m_x;
	vector<double> m_y;

	// Data files
	string m_folder;
	string m_outMean;

	ofstream output;

	double m_Vt;
	string m_parent;

public:
	// ---> Constructors
	TwoDOrderParam() {m_mean.resize(2,0);};
// 	TwoDOrderParam(int minBinNb, int nbOfMeas, int ACInterval,string filename,const char dataType[]);

	// ---> Setters
// 	void SetParam(vector<double>& data, int minBinNb, int binCount, int binSize, int fillBin, double mean, int nbOfMeas, int ACInterval);

	void Initialize(int row_number, int col_number);

	void SetFolder(const char folder[]);
	void SetType(const char dataType[]);
	void CreateFiles(const char folder[],const char dataType[]);

	void SetVt(double Vt);
	void GetParent(string& filename);

	// ---> Procedures
	// >>>> Add data to the mean/binning
	void AddData(vector<double> & measurement);

	// >>>> Get mean values (normal and RMS)
	void GetMean(vector<double> & output);
	int  GetNbOfMeas();
	void RestartMean(int extraMeas);
	void RestartMean();
	void ConvertMean();
	void RaiseMeasures(int extraMeas);

	int GetDataSize()
	{
		return m_row_number*m_col_number;
	}

	double GetVt()
	{
		return m_Vt;
	}

	void Export2D_Data(string& outputString);
	void Export2D_Matrix(string& outputString, int nxColMax, int nyRowMax);
	void ReadPositions(string& posTable);
	void SetPositions(input_params& input);

	void GetSuperlatticeNf(vector<double>& output,int nx, int ny);
	void RectManualGetSuperlatticeNf(vector<double>& output,int nx, int ny, vector<int>& neighTable);
	void RhombusManualGetSuperlatticeNf(vector<double>& output,vector<int>& SubSize, int nx, int ny, vector<int>& neighTable);
};

#endif
