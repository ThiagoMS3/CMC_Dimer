#ifndef MCPARAMETERS_H_
#define MCPARAMETERS_H_

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
// #include "SysConf.h"

// Definitions
using namespace std;
using namespace boost;

class MCParameters
{
// Method used to print / read the class elements. Uses "boost::serialization"
// !!! This library must be compiled beforehand !!!
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & beta;
		ar & Kz;
		ar & Kt;
		ar & q;
		ar & NbMaxUpdate;
		ar & Vt;
		ar & folderName;
		ar & NbMeasures;
		ar & MeasuresToTake;

		ar & alpha0;
		ar & alpha3;

		ar & shift0;
		ar & shift3;
	}

public:
	// ---> Members
	double	beta;
	double	Kz;
	double	Kt;
	double	q;

	double alpha0;
	double alpha3;

	double shift0;
	double shift3;

	int 	NbMaxUpdate;

	double	Vt;
	string	folderName;

	int 	NbMeasures;
	int 	MeasuresToTake;

	// ---> Constructors
	MCParameters() {}


	MCParameters(const input_params& input,int L, const char outputFolder[])
	{
		SetParam(input,L,outputFolder);
	}

	MCParameters(MCParameters& toCopy)
	{
		Copy(toCopy);
	}

	// ---> Methods
	void SetParam(const input_params& input,int L, const char outputFolder[])
	{
		beta	= input.Dbeta * input.N;
		Kz		= input.Dbeta;
		Kt		= -log(tanh((input.Dbeta)))/2.;
		q 		= exp(-2*Kt);

		if(input.NbMaxUpdate==0)
			NbMaxUpdate = L*input.N;
		else
			NbMaxUpdate = input.NbMaxUpdate;

		NbMeasures = input.NbMeasures;
		MeasuresToTake = input.MeasuresToTake;
		folderName	= outputFolder;

		alpha0 = input.alpha0;
		alpha3 = input.alpha3;

		shift0 = input.shift0;
		shift3 = input.shift3;
	}

	void UpdateData(const input_params& input)
	{
		NbMeasures = input.NbMeasures;
		MeasuresToTake = input.MeasuresToTake;

		beta	= input.Dbeta * input.N;
		Kz		= input.Dbeta;
		Kt		= -log(tanh((input.Dbeta)))/2.;
		q 		= exp(-2*Kt);
	}

	void Copy(MCParameters& toCopy)
	{
		beta = toCopy.beta;
		Kz = toCopy.Kz;
		Kt = toCopy.Kt;
		q = toCopy.q;
		NbMaxUpdate = toCopy.NbMaxUpdate;
		Vt = toCopy.Vt;
		folderName = toCopy.folderName;
		NbMeasures = toCopy.NbMeasures;
		MeasuresToTake = toCopy.MeasuresToTake;

		alpha0 = toCopy.alpha0;
		alpha3 = toCopy.alpha3;

		shift0 = toCopy.shift0;
		shift3 = toCopy.shift3;
	}

	void SetVt(double VtIn)
	{
		Vt = VtIn;
	}

	void PrintDebug(ostream& out)
	{
		out << " ### DEBUG Print MCParameters ###" << endl;
		out << " beta                   = " << beta << endl;
		out << " Kz                     = " << Kz << endl;
		out << " Kt                     = " << Kt << endl;
		out << " q                      = " << q << endl;
		out << " NbMaxUpdate            = " << NbMaxUpdate << endl;
		out << " V/t                    = " << Vt << endl;
		out << " NbMeasures             = " << NbMeasures << endl;
		out << " MeasuresToTake         = " << MeasuresToTake << endl;
	}
};

#endif
