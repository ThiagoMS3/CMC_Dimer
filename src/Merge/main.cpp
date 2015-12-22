#include "main.h"
/*
 *    Classical Monte Carlo simulation for the QDM - Local updates version
 *    >>> Based on Thomas' notes (14 May 2013)
 *
 */

void MergeMCOutput(const char filename[], vector<string> &folders,string& WorkFolder)
{
	string dummy = filename;
	string inputFilename;
	string outputFilename = WorkFolder + "/Final_" + dummy;

	ifstream inputs;
	ofstream output;


// 	Set up V/t
	double Vt = 0;
	double value = 0;

	output.open(outputFilename.c_str(),ios::trunc);
	output.precision(15);

	for(uint rank = 0; rank < folders.size(); ++rank)
	{
		inputFilename = folders[rank] + "/" + dummy;
		inputs.open(inputFilename.c_str());

		inputs >> Vt;
		inputs >> value;

		output << Vt << " " << value << endl;
		output.flush();
		inputs.close();
	}
	output.close();
};

void MergeMCMultiple(const char filename[], vector<string> &folders, string& WorkFolder, int len)
{
	string dummy = filename;
	string inputFilename;
	string outputFilename = WorkFolder + "/Final_" + dummy + ".dat";

	ifstream inputs;
	ofstream output;


// 	Set up V/t
	double Vt = 0;
	double value = 0;

	output.open(outputFilename.c_str(),ios::trunc);
	output.precision(15);

	for(uint rank = 0; rank < folders.size(); ++rank)
	{
		inputFilename = folders[rank] + "/" + dummy + ".odat";
		inputs.open(inputFilename.c_str());

		inputs >> Vt;
		output << Vt;

		for(int iii = 0; iii < len; ++iii)
		{
			inputs >> value;

			output << " " << value;
		}

		output << endl;
		output.flush();
		inputs.close();
	}
	output.close();
};

void MergeAutoCorrOutput(const char filename[], vector<string> &folders,string& WorkFolder, input_params& inputPars)
{
	string dummy = filename;
	string inputFilename;
	string outputFilename = WorkFolder + "/Final_" + dummy;

	ifstream inputs;
	ofstream output;

	vector<vector<double> > dataBuffer(inputPars.AutocorrLength,vector<double>(inputPars.VtN,0));

	for(uint rank = 0; rank < folders.size(); ++rank)
	{
		inputFilename = folders[rank] + "/" + dummy;
		inputs.open(inputFilename.c_str());

		for(int jjj = 0; jjj < inputPars.AutocorrLength; ++jjj)
		{
			inputs >> dataBuffer[jjj][rank];
		}
		inputs.close();
	}

	output.open(outputFilename.c_str(),ios::trunc);
	output.precision(15);
	for(int iii = 0; iii < inputPars.AutocorrLength; ++iii)
	{
		output << iii;
		for(int jjj =0; jjj < inputPars.VtN; ++jjj)
		{
			 output << " " << dataBuffer[iii][jjj];
		}
		output << endl;
	}
	output.close();


};

void MergeBunchOutput(string& filename, vector<string> &folders,string& WorkFolder, input_params& inputPars,int NumberOfBinPoints)
{
	string dummy = filename;
	string inputFilename;
	string outputFilename = WorkFolder + "/Final_" + dummy;

	ifstream inputs;
	ofstream output;

	vector<vector<double> > dataBuffer(NumberOfBinPoints,vector<double>(inputPars.VtN + 1,0));

	for(int rank = 0; rank < inputPars.VtN; ++rank)
	{
		inputFilename = folders[rank] + "/" + dummy;
		inputs.open(inputFilename.c_str());

		for(int jjj = 0; jjj < NumberOfBinPoints; ++jjj)
		{
			inputs >> dataBuffer[jjj][inputPars.VtN];
			inputs >> dataBuffer[jjj][rank];
		}
		inputs.close();
	}

	output.open(outputFilename.c_str(),ios::trunc);
	output.precision(15);
	for(int iii = 0; iii < NumberOfBinPoints; ++iii)
	{
		output << dataBuffer[iii][inputPars.VtN];
		for(int jjj =0; jjj < inputPars.VtN; ++jjj)
		{
			 output << " " << dataBuffer[iii][jjj];
		}
		output << endl;
	}
	output.close();
};

void MergeMC(const char dummyType[], vector<string> &VtFolder,string& WorkFolder, input_params& inputPars)
{
	string dataType = dummyType;
	string ObjFile = "/Obj_Mean" + dataType +".odat";
	string MeanFileName = WorkFolder + "/Final_Mean" + dataType + ".dat";
	string IACFileName = WorkFolder + "/Final_IntegAC_" + dataType + ".dat";
	string fileName;
	cout << " >>>> Will open file " << ObjFile << endl;
	OrderParam MeanObj;

	vector<double> MeanOut(inputPars.VtN,0);
	vector<double> IntegACOut(inputPars.VtN,0);

	ofstream output;

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanOut[iii] = MeanObj.GetMean();
		IntegACOut[iii] = MeanObj.GetIntegAC();
	}

	cout << " ---> Merging means" << endl;

	output.open(MeanFileName.c_str(),ios::trunc);
	output.precision(15);

	double V = 0;
	double U = 0;

	double dummyAlpha = sqrt(inputPars.alpha0*inputPars.alpha0 + inputPars.alpha3*inputPars.alpha3);

	double alpha0 = inputPars.alpha0/dummyAlpha;
	double alpha3 = inputPars.alpha3/dummyAlpha;

	double shift0 = inputPars.shift0;
	double shift3 = inputPars.shift3;

	output << "# V   U   Value" << endl;
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		V = alpha3*inputPars.VtValues[iii] + shift3;
		U = alpha0*inputPars.VtValues[iii] + shift0;
		output << V << "\t" << U << "\t" << MeanOut[iii] << endl;
	}
	output.close();

	cout << " ---> Merging IntegAC" << endl;

	output.open(IACFileName.c_str(),ios::trunc);
	output.precision(15);
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		V = alpha3*inputPars.VtValues[iii] + shift3;
		U = alpha0*inputPars.VtValues[iii] + shift0;
		output << V << "\t" << U << "\t" << IntegACOut[iii] << endl;
	}
	output.close();

	fileName = "AC_" + dataType + ".dat";
	MergeAutoCorrOutput(fileName.c_str(),VtFolder,WorkFolder,inputPars);

	fileName = "err_" + dataType +".dat";
	MergeBunchOutput(fileName,VtFolder,WorkFolder,inputPars,MeanObj.GetNumberOfPoints());
};

void MergeCmplxMC(const char dummyType[], vector<string> &VtFolder,string& WorkFolder, input_params& inputPars)
{
	string dataType = dummyType;
	string ObjFile = "/Obj_Mean" + dataType +".odat";
	string MeanFileName = WorkFolder + "/Final_Mean" + dataType + ".dat";
	string fileName;
	cout << " >>>> Will open file " << ObjFile << endl;
	CmplxOrderParam MeanObj;

	vector<complex<double> > MeanOut(inputPars.VtN,0);

	ofstream output;

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanOut[iii] = MeanObj.GetMean();
	}

	cout << " ---> Merging means" << endl;

	output.open(MeanFileName.c_str(),ios::trunc);
	output.precision(15);

	double V = 0;
	double U = 0;

	double dummyAlpha = sqrt(inputPars.alpha0*inputPars.alpha0 + inputPars.alpha3*inputPars.alpha3);

	double alpha0 = inputPars.alpha0/dummyAlpha;
	double alpha3 = inputPars.alpha3/dummyAlpha;

	double shift0 = inputPars.shift0;
	double shift3 = inputPars.shift3;

	output << "# V   U   Value" << endl;
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		V = alpha3*inputPars.VtValues[iii] + shift3;
		U = alpha0*inputPars.VtValues[iii] + shift0;
		output << V << "\t" << U << "\t" << MeanOut[iii].real() << "\t" << MeanOut[iii].imag() << endl;
	}
	output.close();

	string histOutput;
	int magSize = -1;
	int argSize = -1;
	int idx = -1;
	vector<int> vecHist;

	double deltaMag = -1;
	double deltaArg = -1;

	double magnitude = -1;
	double argument = -1;

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		histOutput = VtFolder[iii] + "_cmplx_hist.dat";
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());
		magSize = MeanObj.m_histMagSize;
		argSize = MeanObj.m_histArgSize;
		deltaMag = MeanObj.m_histDeltaMag;
		deltaArg = MeanObj.m_histDeltaArg;

		MeanObj.GetHist(vecHist);

		output.open(histOutput.c_str(),ios::trunc);
		output.precision(15);

		for(int kkk = 0; kkk < magSize; ++kkk)
		{
			for(int jjj = 0; jjj < argSize; ++jjj)
			{
				magnitude = kkk*deltaMag+deltaMag/2;
				argument  = jjj*deltaArg - M_PI+deltaArg/2;
				idx = idxConv(argSize,kkk,jjj);
				output << magnitude << " " << argument << " " << vecHist[idx] << endl;
			}
			argument = -M_PI+deltaArg/2;
			idx = idxConv(argSize,kkk,0);
			output << magnitude << " " << argument << " " << vecHist[idx] << endl;
			output << endl;
		}
		output.close();

	}

};

void MergeTwoDMC(const char dummyType[], vector<string> &VtFolder,string& WorkFolder, input_params& inputPars)
{
	string dataType = dummyType;
	string ObjFile = "/Obj_" + dataType +".odat";
	string MeanFileName = WorkFolder + "/Final_" + dataType + ".dat";
//	string IACFileName = WorkFolder + "/Final_IntegAC_" + dataType + ".dat";
	string fileName;
	cout << " >>>> Will open file " << ObjFile << endl;

	TwoDOrderParam MeanObj;

	fileName = VtFolder[0] + ObjFile;
	ImportObject(MeanObj,fileName.c_str());

	int len = MeanObj.GetDataSize();
	vector<vector<double> > MeanOut(inputPars.VtN,vector<double>(len,0));

	MeanObj.GetMean(MeanOut[0]);
//	vector<double> IntegACOut(inputPars.VtN,0);

	ofstream output;

	for(int iii = 1; iii < inputPars.VtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanObj.GetMean(MeanOut[iii]);
//		IntegACOut[iii] = MeanObj.GetIntegAC();
	}

	cout << " ---> Merging means" << endl;

	double V = 0;
	double U = 0;

	double dummyAlpha = sqrt(inputPars.alpha0*inputPars.alpha0 + inputPars.alpha3*inputPars.alpha3);

	double alpha0 = inputPars.alpha0/dummyAlpha;
	double alpha3 = inputPars.alpha3/dummyAlpha;

	double shift0 = inputPars.shift0;
	double shift3 = inputPars.shift3;

	output.open(MeanFileName.c_str(),ios::trunc);
	output.precision(15);
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		V = alpha3*inputPars.VtValues[iii] + shift3;
		U = alpha0*inputPars.VtValues[iii] + shift0;
		output << V << "\t" << U << "\t";

		for(int jjj = 0; jjj < len; ++jjj)
		{
				output << " " << MeanOut[iii][jjj];
		}
		output << endl;
	}
	output.close();

};

void MergeCorrMC(const char dummyType[], vector<string> &VtFolder,string& WorkFolder, input_params& inputPars,int CorrSize)
{
	string dataType = dummyType;
	string CorrObjFile = "/Obj_Corr" + dataType +".odat";

	string fileName;
	cout << " >>>> Will open file " << CorrObjFile << endl;
	OrderParam 		CorrObj;
	TwoDOrderParam 	MeanObj;

	vector<vector<double> > CorrOut(inputPars.VtN,vector<double>(CorrSize,0));

	vector<double> ConstMean(inputPars.VtN,0);

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		fileName = VtFolder[iii] + CorrObjFile;
		ImportObject(CorrObj,fileName.c_str());
		CorrObj.GetCorr(CorrOut[iii],CorrSize);
	}

	cout << " ---> Merging correlations" << endl;

	string MeanFileName = WorkFolder + "/Final_Corr" + dataType + ".dat";
	ofstream output(MeanFileName.c_str(),ios::trunc);
	output.precision(15);
	for(int jjj = 0; jjj < CorrSize; ++jjj)
	{
		output << jjj;
		for(int iii = 0; iii < inputPars.VtN; ++iii)
		{
			output << " " << CorrOut[iii][jjj];
		}
		output << endl;
	}
	output.close();
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

	MergeMC("Mag",vtFolder,WorkFolder,inputPars);
	MergeMC("Energy",vtFolder,WorkFolder,inputPars);
	MergeMC("QEnergy",vtFolder,WorkFolder,inputPars);
	MergeMC("QNewEnergy",vtFolder,WorkFolder,inputPars);
	MergeMC("PotEnergy",vtFolder,WorkFolder,inputPars);
	MergeMC("KinEnergy",vtFolder,WorkFolder,inputPars);
	MergeMC("NewKinEnergy",vtFolder,WorkFolder,inputPars);
	MergeMC("N0",vtFolder,WorkFolder,inputPars);
	MergeMC("N1",vtFolder,WorkFolder,inputPars);
	MergeMC("N2",vtFolder,WorkFolder,inputPars);
	MergeMC("N3",vtFolder,WorkFolder,inputPars);

	if(inputPars.runType.compare("Part")!=0)
	{
		MergeMC("SubA",vtFolder,WorkFolder,inputPars);
		MergeMC("SubB",vtFolder,WorkFolder,inputPars);
		MergeMC("SubC",vtFolder,WorkFolder,inputPars);
	}

	if(inputPars.runType.compare("Moessner")==0)
	{
		MergeTwoDMC("SpatialCorrN3N3",vtFolder,WorkFolder,inputPars);
		MergeTwoDMC("MeanLocalN3",vtFolder,WorkFolder,inputPars);

#ifdef COMPLEX_PARAMS
		MergeMC("SymmetryParameter",vtFolder,WorkFolder,inputPars);
		MergeMC("Radius",vtFolder,WorkFolder,inputPars);
		MergeCmplxMC("ComplexPhase",vtFolder,WorkFolder,inputPars);

		MergeMC("SymmetryParameterPerLayer",vtFolder,WorkFolder,inputPars);
		MergeMC("RadiusPerLayer",vtFolder,WorkFolder,inputPars);
		MergeCmplxMC("ComplexPhasePerLayer",vtFolder,WorkFolder,inputPars);
#endif
	}

	MergeCorrMC("SzSz",vtFolder,WorkFolder,inputPars,inputPars.N);
	MergeCorrMC("SxSx",vtFolder,WorkFolder,inputPars,inputPars.N);
	MergeCorrMC("Dimer",vtFolder,WorkFolder,inputPars,inputPars.N);

	if(inputPars.ny==1)
	{
		MergeMC("StagMag",vtFolder,WorkFolder,inputPars);
		MergeMCMultiple("MeanPart",vtFolder,WorkFolder,inputPars.nx);
		MergeMC("ChainEnergy",vtFolder,WorkFolder,inputPars);
	}

	if(inputPars.ny!=1&&inputPars.runType.compare("Part")==0&&inputPars.initCondType==0)
	{
		MergeMC("CoreMag",vtFolder,WorkFolder,inputPars);
	}

	MergeMCOutput("Accept.dat",vtFolder,WorkFolder);
	MergeMCOutput("ClusterNumber.dat",vtFolder,WorkFolder);
	MergeMCOutput("ClusterSize.dat",vtFolder,WorkFolder);

	return 0;
}
