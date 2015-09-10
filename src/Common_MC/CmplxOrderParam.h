#ifndef CMPLXORDERPARAM_H_
#define CMPLXORDERPARAM_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/complex.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

#include "dataStructures.h"
#include "customMath.h"
#include "MCInputParams.h"


// Definitions
using namespace std;
using namespace boost;

// --- Big Huge Class for order parameters ---
// -------------------------------------------
class CmplxOrderParam
{
private:
// Method used to print / read the class elements. Uses "boost::serialization"
// !!! This library must be compiled beforehand !!!
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & 	m_Type;

		ar & 	m_hist;

		ar & 	m_histArgSize;
		ar & 	m_histMagSize;
		ar & 	m_histMaxMag;

		ar & 	m_histDeltaMag;
		ar & 	m_histDeltaArg;

		ar & 	m_mean;
		ar & 	m_nbOfMeas;

		ar & 	m_folder;
		ar & 	m_outMean;

		ar & 	m_Vt;
	}

	// Type
	string m_Type;

	// Histogram


	// Mean values
	complex<double>  			m_mean;			// Mean value
	int 						m_nbOfMeas;		// Number of measurements

	// Data files
	string 						m_folder;
	string 						m_outMean;
	ofstream 					output;

	double 						m_Vt;

public:
	// ---> Constructors
	CmplxOrderParam() {}

	int 						m_histArgSize;
	int 						m_histMagSize;
	double						m_histMaxMag;

	double						m_histDeltaMag;
	double						m_histDeltaArg;

	vector<int> 				m_hist;			// Base of the histogram

	CmplxOrderParam(input_params& input,int histArgSize, int histMagSize, double histMaxMag);

	// ---> Setters
	void 			SetFolder(const char folder[]);
	void 			SetType(const char dataType[]);
	void 			CreateFiles(const char folder[],const char dataType[]);

	void 			SetVt(double Vt);

	// ---> Procedures
	// >>>> Add data to the mean/binning
	void 			AddData(complex<double>  measurement);

	// >>>> Reset the mean value for a new series of measurements
	void 			RestartMean(int extraMeas);
	void 			RestartMean();
	void 			RaiseMeasures(int extraMeas);

	// >>>> Get mean values (normal and RMS)
	complex<double> GetMean();
	complex<double> ConvertMean();
	int 			GetNbOfMeas();

	void			GetHist(vector<int> & histOut);

	double GetVt()
	{
		return m_Vt;
	}
};

#endif
