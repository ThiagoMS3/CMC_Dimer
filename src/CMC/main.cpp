#include "main.h"

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// >>>> Initialize
	// ---> This program expects 2 arguments : the configuration file, and the MC parameters file
	if(argc<3)
	{
		cout << "Too few arguments : need a configuration file and a MC parameters file" << endl;
		return 0;
	}

	if(argc>3)
	{
		cout << "Too many arguments : need a configuration file and a MC parameters file" << endl;
		return 0;
	}

	// ---> File/ folder names
	string confFile		= argv[1];	// Spin cofiguration file
	string WorkFolder		= argv[2];	// Working folder

	string outConf 		= WorkFolder + "/Obj_conf.odat";	// Object files
	string MCparametersFile	= WorkFolder + "/Obj_MCparams.odat";	//
	string RNGengineFile	= WorkFolder + "/Obj_Engine.odat";	//

	string OrderParamFile;	// (dummy) Order parameters filename

	// ---> Import objects
	// Configuration
	SysConf 	conf;
	cout << " -> Importing " << confFile << endl;
	ImportObject(conf,confFile.c_str());

	// RNG
	cout << " -> Importing " << RNGengineFile << endl;
	conf.ImportRNG(RNGengineFile.c_str());

	// MC parameters
	MCParameters 	parameters;
	cout << " -> Importing " << MCparametersFile << endl;
	ImportObject(parameters,MCparametersFile.c_str());

	// Order parameters
	OrderParam	MeanMag;
	OrderParam	MeanEnergy;
	OrderParam	MeanN0;
	OrderParam	MeanN1;
	OrderParam	MeanN2;
	OrderParam	MeanN3;
	OrderParam	MeanQEnergy;
	OrderParam	MeanPotEnergy;
	OrderParam	MeanKinEnergy;
	OrderParam	CorrN3;

	OrderParam	MeanStagMag;
	OrderParam	MeanChainEnergy;
	OrderParam	MeanCoreMag;

	OrderParam	CorrMag;
	OrderParam	CorrStagMag;

	// 2D Order parameters

	TwoDOrderParam	MeanLocalDimer;
	TwoDOrderParam	MeanLocalNf;
	TwoDOrderParam	MeanStar3;
	TwoDOrderParam	MeanPartition;

	/*
	 *    > Name of the files
	 *    > Import it
	 *    > Update number of measurements
	 */

	OrderParamFile = WorkFolder + "/Obj_MeanMag.odat";
	ImportObject(MeanMag,OrderParamFile.c_str());
	MeanMag.ResetRMSMean();
	MeanMag.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanEnergy.odat";
	ImportObject(MeanEnergy,OrderParamFile.c_str());
	MeanEnergy.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanN0.odat";
	ImportObject(MeanN0,OrderParamFile.c_str());
	MeanN0.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanN1.odat";
	ImportObject(MeanN1,OrderParamFile.c_str());
	MeanN1.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanN2.odat";
	ImportObject(MeanN2,OrderParamFile.c_str());
	MeanN2.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanN3.odat";
	ImportObject(MeanN3,OrderParamFile.c_str());
	MeanN3.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanQEnergy.odat";
	ImportObject(MeanQEnergy,OrderParamFile.c_str());
	MeanQEnergy.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanPotEnergy.odat";
	ImportObject(MeanPotEnergy,OrderParamFile.c_str());
	MeanPotEnergy.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanKinEnergy.odat";
	ImportObject(MeanKinEnergy,OrderParamFile.c_str());
	MeanKinEnergy.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanCorrN3.odat";
	ImportObject(CorrN3,OrderParamFile.c_str());
	CorrN3.RestartMean(parameters.MeasuresToTake);

	if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
	{
		if(conf.ny==1)
		{
			OrderParamFile = WorkFolder + "/Obj_MeanStagMag.odat";
			ImportObject(MeanStagMag,OrderParamFile.c_str());
			MeanStagMag.RestartMean(parameters.MeasuresToTake);

			OrderParamFile = WorkFolder + "/Obj_MeanChainEnergy.odat";
			ImportObject(MeanChainEnergy,OrderParamFile.c_str());
			MeanChainEnergy.RestartMean(parameters.MeasuresToTake);

			OrderParamFile = WorkFolder + "/Obj_CorrMag.odat";
			ImportObject(CorrMag,OrderParamFile.c_str());
			CorrMag.RestartCorr(parameters.MeasuresToTake);

			OrderParamFile = WorkFolder + "/Obj_CorrStagMag.odat";
			ImportObject(CorrStagMag,OrderParamFile.c_str());
			CorrStagMag.RestartCorr(parameters.MeasuresToTake);

			OrderParamFile = WorkFolder + "/Obj_MeanPart.odat";
			ImportObject(MeanPartition,OrderParamFile.c_str());
			MeanPartition.RestartMean(parameters.MeasuresToTake);
		}
		else
		{
			OrderParamFile = WorkFolder + "/Obj_MeanCoreMag.odat";
			ImportObject(MeanCoreMag,OrderParamFile.c_str());
			MeanCoreMag.RestartMean(parameters.MeasuresToTake);
		}
	}

	OrderParamFile = WorkFolder + "/Obj_MeanLocalNf.odat";
	ImportObject(MeanLocalNf,OrderParamFile.c_str());
	MeanLocalNf.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_MeanLocalDimer.odat";
	ImportObject(MeanLocalDimer,OrderParamFile.c_str());
	MeanLocalDimer.RestartMean(parameters.MeasuresToTake);

	OrderParamFile = WorkFolder + "/Obj_Star3.odat";
	ImportObject(MeanStar3,OrderParamFile.c_str());
	MeanStar3.RestartMean(parameters.MeasuresToTake);

	// >>>> MC
	// ---> Preparations

	// Extra file names
	string 	outAccept 		= WorkFolder + "/Accept.dat";		// Accepance rates
	string 	outMeanClusterSize 	= WorkFolder + "/ClusterSize.dat";	// Mean cluster size
	string 	outMeanClusterNumber 	= WorkFolder + "/ClusterNumber.dat";	// Mean cluster number

	// Counters
	int 	NbUpdates 		= 0;		// Number of updates done
	int 	MeasureIndex 		= 0;

	int 	FinalNumberOfClusters	= 0;
	double 	FinalMeanClusterSize 	= 0;

	// Some auxiliary variables
	double 		dummyParam = 0;

	vector<double>	dummyCorrMag(conf.nx + conf.p,0);
	vector<double>	dummyCorrStagMag(conf.nx + conf.p,0);
	vector<double> dummyNf(4,-1);
	vector<double>	dummyMeanLocalNf(conf.L);
	vector<double> dummyMeanLocalDimer(conf.L*conf.NbOfNeights);
	vector<double> dummyMeanStar3(conf.L,0);
	vector<double>	dummyPart(conf.nx,0);

	int 		OrderBunch;

	// ---> Run (your fools) !
	conf.MC_Prepare(parameters);

	while(MeasureIndex < parameters.NbMeasures)
	{
		conf.ClusterIter();
		NbUpdates 			+= conf.GetClusterSize();
		FinalMeanClusterSize	+= conf.GetClusterSize();
		FinalNumberOfClusters	+= 1;

		if(NbUpdates>=parameters.NbMaxUpdate)
		{
			if(MeasureIndex>=parameters.NbMeasures-parameters.MeasuresToTake)
			{
				// Do measurements
				dummyParam 		= conf.GetMagnetization();
				MeanMag.AddData(dummyParam);

				dummyParam		= conf.GetEnergy()/(conf.N);
				MeanEnergy.AddData(dummyParam);

				conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer);
				MeanN0.AddData(dummyNf[0]);
				MeanN1.AddData(dummyNf[1]);
				MeanN2.AddData(dummyNf[2]);
				MeanN3.AddData(dummyNf[3]);

				dummyParam		= conf.GetQEnergy_N3();
				MeanQEnergy.AddData(dummyParam);

				dummyParam		= conf.GetPotEnergy_N3();
				MeanPotEnergy.AddData(dummyParam);

				dummyParam		= conf.GetKinEnergy();
				MeanKinEnergy.AddData(dummyParam);

				MeanLocalNf.AddData(dummyMeanLocalNf);
				MeanLocalDimer.AddData(dummyMeanLocalDimer);

				conf.GetN3Correlation(dummyParam,dummyMeanStar3);
				dummyParam = dummyParam/dummyNf[3];

				CorrN3.AddData(dummyParam);
				MeanStar3.AddData(dummyMeanStar3);

				if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
				{
					if(conf.ny==1)
					{
						conf.BuildEquivalentSpinChain();

						dummyParam = conf.CalculateStagMag();
						MeanStagMag.AddData(dummyParam);

						conf.GetMeanPart(dummyPart);
						MeanPartition.AddData(dummyPart);

						dummyParam = conf.CalculateChainEnergy();
						MeanChainEnergy.AddData(dummyParam);

						conf.CalculateChainCorrelations(dummyCorrMag,dummyCorrStagMag);
						CorrMag.AddCorr(dummyCorrMag);
						CorrStagMag.AddCorr(dummyCorrStagMag);
					}
					else
					{
						dummyParam = conf.GetCoreMagnetization();
						MeanCoreMag.AddData(dummyParam);
					}
				}
			}

			++MeasureIndex;
			NbUpdates = 0;
		}
	}

	// >>>> Post-QMC calculations
	// > Binning / bunching and autocorrelations
//	MeanMag.CalculateRMSMean();
	MeanMag.CalculateError(OrderBunch);
	MeanMag.PrintError();

	MeanEnergy.CalculateError(OrderBunch);
	MeanEnergy.PrintError();

	MeanN0.CalculateError(OrderBunch);
	MeanN0.PrintError();

	MeanN1.CalculateError(OrderBunch);
	MeanN1.PrintError();

	MeanN2.CalculateError(OrderBunch);
	MeanN2.PrintError();

	MeanN3.CalculateError(OrderBunch);
	MeanN3.PrintError();

	MeanQEnergy.CalculateError(OrderBunch);
	MeanQEnergy.PrintError();

	MeanPotEnergy.CalculateError(OrderBunch);
	MeanPotEnergy.PrintError();

	MeanKinEnergy.CalculateError(OrderBunch);
	MeanKinEnergy.PrintError();

	CorrN3.CalculateError(OrderBunch);
	CorrN3.PrintError();

	if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
	{
		if(conf.ny==1)
		{
			MeanStagMag.CalculateError(OrderBunch);
			MeanStagMag.PrintError();

			MeanChainEnergy.CalculateError(OrderBunch);
			MeanChainEnergy.PrintError();
		}
		else
		{
			MeanCoreMag.CalculateRMSMean();
			MeanCoreMag.CalculateError(OrderBunch);
			MeanCoreMag.PrintError();
		}
	}

	ofstream output;

	output.open(outAccept.c_str(), ios::out | ios::app);
	output << conf.Vt << " " << (double)conf.GetTotalAccept()/FinalNumberOfClusters << endl;
	output.close();

	output.open(outMeanClusterNumber.c_str(), ios::out | ios::app);
	output << conf.Vt << " " << (double)FinalNumberOfClusters << endl;
	output.close();

	output.open(outMeanClusterSize.c_str(), ios::out | ios::app);
	output << conf.Vt << " " << (double)FinalMeanClusterSize/FinalNumberOfClusters << endl;
	output.close();

	// >>>> Export data
	OrderParamFile = WorkFolder + "/Obj_MeanMag.odat";
	ExportObject(MeanMag,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanEnergy.odat";
	ExportObject(MeanEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN0.odat";
	ExportObject(MeanN0,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN1.odat";
	ExportObject(MeanN1,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN2.odat";
	ExportObject(MeanN2,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN3.odat";
	ExportObject(MeanN3,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanQEnergy.odat";
	ExportObject(MeanQEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanPotEnergy.odat";
	ExportObject(MeanPotEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanKinEnergy.odat";
	ExportObject(MeanKinEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanCorrN3.odat";
	ExportObject(CorrN3,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanLocalDimer.odat";
	ExportObject(MeanLocalDimer,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanLocalNf.odat";
	ExportObject(MeanLocalNf,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_Star3.odat";
	ExportObject(MeanStar3,OrderParamFile.c_str());

	if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
	{
		if(conf.ny==1)
		{
			OrderParamFile = WorkFolder + "/Obj_MeanStagMag.odat";
			ExportObject(MeanStagMag,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanPart.odat";
			ExportObject(MeanPartition,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanChainEnergy.odat";
			ExportObject(MeanChainEnergy,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_CorrMag.odat";
			ExportObject(CorrMag,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_CorrStagMag.odat";
			ExportObject(CorrStagMag,OrderParamFile.c_str());
		}
		else
		{
			OrderParamFile = WorkFolder + "/Obj_MeanCoreMag.odat";
			ExportObject(MeanCoreMag,OrderParamFile.c_str());
		}
	}

	ExportObject(conf,outConf.c_str());
	conf.ExportRNG(RNGengineFile.c_str());

	return 0;
}
