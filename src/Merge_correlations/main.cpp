#include "main.h"
/*
 *    Classical Monte Carlo simulation for the QDM - Local updates version
 *    >>> Based on Thomas' notes (14 May 2013)
 *
 */

void MergeTwoDCorrMC(const char dummyCorrType[], const char dummyMeanType[], vector<string> &VtFolder,string& WorkFolder, input_params& inputPars)
{
	int L = inputPars.nx*inputPars.ny;

	string dataCorrType = dummyCorrType;
	string dataMeanType = dummyMeanType;

	string ObjCorrFile = "/Obj_" + dataCorrType +".odat";
	string ObjMeanFile = "/Obj_" + dataMeanType +".odat";

	string CovarianceFileName = WorkFolder + "/Final_" + dataCorrType + "_Covariance.dat";
	string CorrAFileName = WorkFolder + "/Final_" + dataCorrType + "_Corr_SubA.dat";
	string CorrBFileName = WorkFolder + "/Final_" + dataCorrType + "_Corr_SubB.dat";
	string CorrCFileName = WorkFolder + "/Final_" + dataCorrType + "_Corr_SubC.dat";

	string SqrDistancesFile = WorkFolder + "/All_SqrDistances.dat";
	string SqrDistancesKeyFile = WorkFolder + "/SqrDistances_Key.dat";

	vector<vector<int> >	SqrDistances(L,vector<int>(L,-1));
	vector<int>				SqrDistancesKey;

	ifstream 	SqrDistancesStream(SqrDistancesFile.c_str());
	ifstream 	SqrDistancesKeyStream(SqrDistancesKeyFile.c_str());

	int intDummy;

	while(!SqrDistancesKeyStream.eof())
	{
		SqrDistancesKeyStream >> intDummy;
		SqrDistancesKey.push_back(intDummy);
	}

	for(int iii = 0; iii < L; ++iii)
	{
		for(int jjj = 0; jjj < L; ++jjj)
		{
			SqrDistancesStream >> SqrDistances[iii][jjj];
		}
	}

	SqrDistancesStream.close();
	SqrDistancesKeyStream.close();

	string fileName;

	cout << " >>>> Will open files " << ObjCorrFile << " and " << ObjMeanFile << endl;

	TwoDOrderParam CorrObj;
	TwoDOrderParam MeanObj;

	fileName = VtFolder[0] + ObjMeanFile;
	ImportObject(CorrObj,fileName.c_str());

	fileName = VtFolder[0] + ObjCorrFile;
	ImportObject(MeanObj,fileName.c_str());

	int CorrLen = CorrObj.GetDataSize();
	int MeanLen = MeanObj.GetDataSize();

	int nbOfDists = CorrLen/MeanLen;

	vector<vector<double> > CorrOut(inputPars.VtN,vector<double>(CorrLen,0));
	vector<vector<double> > MeanOut(inputPars.VtN,vector<double>(MeanLen,0));

	CorrObj.GetMean(MeanOut[0]);
	MeanObj.GetMean(MeanOut[0]);

	for(int iii = 1; iii < inputPars.VtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjMeanFile;
		ImportObject(CorrObj,fileName.c_str());

		fileName = VtFolder[iii] + ObjCorrFile;
		ImportObject(MeanObj,fileName.c_str());

		CorrObj.GetMean(MeanOut[iii]);
		MeanObj.GetMean(MeanOut[iii]);
	}

	// Now, prepare the data
	vector<vector<double> > CovarianceOut(inputPars.VtN,vector<double>(nbOfDists,0));
	vector<vector<double> > CorrSubAOut(inputPars.VtN,vector<double>(nbOfDists,0));
	vector<vector<double> > CorrSubBOut(inputPars.VtN,vector<double>(nbOfDists,0));
	vector<vector<double> > CorrSubCOut(inputPars.VtN,vector<double>(nbOfDists,0));

	double dummyDouble = -1;
	map<int,double>			 DummyDistMap;
	for(int iii = 0; iii < nbOfDists; ++iii)
	{
		DummyDistMap.insert(pair<int,double>(SqrDistancesKey[iii],0));
	}

	vector<map<int,double> > MeanN3PerDist(inputPars.VtN,DummyDistMap);
	int SqrDist = 0;

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{

		// Calculate < n_3_d >

		/* For each value of V/t*/
		for(int jjj = 0; jjj < L; ++jjj)
		{
			/* For each site at dist 0 */
			for(int kkk = 0; kkk < L; ++kkk)
			{
				/* For each site at dist d */
				SqrDist = SqrDistances[jjj][kkk];
				dummyDouble = MeanOut[iii][kkk];
				MeanN3PerDist[iii][SqrDist] = MeanN3PerDist[iii][SqrDist] + dummyDouble;
			}
		}

		for(int jjj = 0; jjj < nbOfDists; ++jjj)
		{
			MeanN3PerDist[iii][SqrDistancesKey[jjj]]=MeanN3PerDist[iii][SqrDistancesKey[jjj]]/MeanN3PerDist[iii].count(SqrDistancesKey[jjj]);
		}

		// Calculate <

		/* For each value of V/t*/
		for(int jjj = 0; jjj < L; ++jjj)
		{
			/* For each site at dist 0 */
			for(int kkk = 0; kkk < L; ++kkk)
			{
				/* For each site at dist d */
				SqrDist = SqrDistances[jjj][kkk];
				dummyDouble = CorrOut[iii][kkk];
				MeanN3PerDist[iii][SqrDist] = MeanN3PerDist[iii][SqrDist] + dummyDouble;
			}
		}
	}

//	ofstream output;

//	cout << " ---> Merging means" << endl;
//
//	double V = 0;
//	double U = 0;
//
//	double dummyAlpha = sqrt(inputPars.alpha0*inputPars.alpha0 + inputPars.alpha3*inputPars.alpha3);
//
//	double alpha0 = inputPars.alpha0/dummyAlpha;
//	double alpha3 = inputPars.alpha3/dummyAlpha;
//
//	double shift0 = inputPars.shift0;
//	double shift3 = inputPars.shift3;
//
//	output.open(MeanFileName.c_str(),ios::trunc);
//	output.precision(15);
//	for(int iii = 0; iii < inputPars.VtN; ++iii)
//	{
//		V = alpha3*inputPars.VtValues[iii] + shift3;
//		U = alpha0*inputPars.VtValues[iii] + shift0;
//		output << V << "\t" << U << "\t";
//
//		for(int jjj = 0; jjj < len; ++jjj)
//		{
//				output << " " << MeanOut[iii][jjj];
//		}
//		output << endl;
//	}
//	output.close();

//	cout << " ---> Merging IntegAC" << endl;
//
//	output.open(IACFileName.c_str(),ios::trunc);
//	output.precision(15);
//	for(int iii = 0; iii < inputPars.VtN; ++iii)
//	{
//		output << inputPars.VtValues[iii] << " " << IntegACOut[iii] << endl;
//	}
//	output.close();

//	fileName = "AC_" + dataType + ".dat";
//	MergeAutoCorrOutput(fileName.c_str(),VtFolder,WorkFolder,inputPars);
//
//	fileName = "err_" + dataType +".dat";
//	MergeBunchOutput(fileName,VtFolder,WorkFolder,inputPars,MeanObj.GetNumberOfPoints());
};

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// >>>> Initialize
	// ---> This program expects 2 arguments : the configuration file, and the MC parameters file
	if(argc<3)
	{
		cout << "Too few arguments : need a input file and a working folder" << endl;
		return 0;
	}

	if(argc>3)
	{
		cout << "Too many arguments : need a input file and a working folder" << endl;
		return 0;
	}

	// ---> Read config files
	string inputFile 		= argv[1];
	string WorkFolder 		= argv[2];

	input_params inputPars;

	ImportObject(inputPars,inputFile.c_str());

	// >>>> Set up file names
	// ---> Folders
	vector<string> 		  vtFolder(inputPars.VtN,WorkFolder + "/Vt_");
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		vtFolder[iii] = vtFolder[iii] + ToString(inputPars.VtValues[iii]);
		mkdir(vtFolder[iii].c_str(),0777);
	}

	MergeTwoDCorrMC("SpatialCorrN3N3","MeanLocalN3",vtFolder,WorkFolder,inputPars);

	return 0;
}
