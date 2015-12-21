#include "main.h"

void ConvertTime(int time_input, int& time_interval_sec, int& time_interval_min, int& time_interval_hours, int& time_interval_days)
{
	int time_interval = time_input;
	time_interval_sec = time_interval % 60;
	time_interval = time_interval/60; // Minutes

	time_interval_min = time_interval % 60;
	time_interval = time_interval/60; // Hours

	time_interval_hours = time_interval % 24;
	time_interval = time_interval/24; // Days

	time_interval_days = time_interval;
}

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
	string WorkFolder	= argv[2];	// Working folder

	string outConf 			= WorkFolder + "/Obj_conf.odat";	// Object files
	string MCparametersFile	= WorkFolder + "/Obj_MCparams.odat";	//
	string RNGengineFile	= WorkFolder + "/Obj_Engine.odat";	//

	string OrderParamFile;	// (dummy) Order parameters filename

	// ---> Import objects
	// Configuration
	SysConf 	conf;
	ImportObject(conf,confFile.c_str());

	// RNG
	conf.ImportRNG(RNGengineFile.c_str());

	// MC parameters
	MCParameters 	parameters;
	ImportObject(parameters,MCparametersFile.c_str());

	// Order parameters
	OrderParam	MeanMag;
	OrderParam	MeanEnergy;

	OrderParam	MeanN0;
	OrderParam	MeanN1;
	OrderParam	MeanN2;
	OrderParam	MeanN3;

#ifdef COMPLEX_PARAMS
	CmplxOrderParam	MeanComplexPhase;
	OrderParam	MeanSymmetryParameter;
	OrderParam	MeanRadius;

	CmplxOrderParam	MeanComplexPhasePerLayer;
	OrderParam	MeanSymmetryParameterPerLayer;
	OrderParam	MeanRadiusPerLayer;
#endif


	OrderParam	MeanQEnergy;
	OrderParam	MeanQKinEnergy;
	OrderParam	MeanQPotEnergy;

	OrderParam      MeanQNewEnergy;
	OrderParam      MeanQNewKinEnergy;
	
	OrderParam	CorrN3;

	OrderParam	SzSzCorrelation;
	TwoDOrderParam	SzSzCorrelationLocalMean;

	OrderParam	DimerCorrelation;
	TwoDOrderParam	DimerCorrelationLocalMean;

	OrderParam	MeanStagMag;
	OrderParam	MeanCoreMag;

	OrderParam	CorrMag;
	OrderParam	CorrStagMag;

	// 2D Order parameters
	TwoDOrderParam	MeanLocalDimer;
	TwoDOrderParam	MeanPart;
	TwoDOrderParam	MeanLocalNf;
	TwoDOrderParam  MeanKinEnergyDensity;
	TwoDOrderParam	SpatialCorrN3N3;
	TwoDOrderParam	MeanLocalN3;

	OrderParam		MeanSubA;
	OrderParam		MeanSubB;
	OrderParam		MeanSubC;

	/*
	 *    > Name of the files
	 *    > Import it
	 *    > Update number of measurements
	 */

	if(parameters.MeasuresToTake!=0)
	{
		OrderParamFile = WorkFolder + "/Obj_MeanMag.odat";
		ImportObject(MeanMag,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanEnergy.odat";
		ImportObject(MeanEnergy,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanN0.odat";
		ImportObject(MeanN0,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanN1.odat";
		ImportObject(MeanN1,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanN2.odat";
		ImportObject(MeanN2,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanN3.odat";
		ImportObject(MeanN3,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanQEnergy.odat";
		ImportObject(MeanQEnergy,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanKinEnergy.odat";
		ImportObject(MeanQKinEnergy,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanQNewEnergy.odat";
		ImportObject(MeanQNewEnergy,OrderParamFile.c_str());
		
		OrderParamFile = WorkFolder + "/Obj_MeanNewKinEnergy.odat";
		ImportObject(MeanQNewKinEnergy,OrderParamFile.c_str());
		
		OrderParamFile = WorkFolder + "/Obj_MeanPotEnergy.odat";
		ImportObject(MeanQPotEnergy,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanKinEnergyDensity.odat";
		ImportObject(MeanKinEnergyDensity,OrderParamFile.c_str());

		if(conf.SimType.compare("Part")==0&&(conf.initCondType == 0 || conf.initCondType == 1))
		{
			OrderParamFile = WorkFolder + "/Obj_MeanPart.odat";
			ImportObject(MeanPart,OrderParamFile.c_str());

			if(conf.ny==1)
			{
				OrderParamFile = WorkFolder + "/Obj_MeanStagMag.odat";
				ImportObject(MeanStagMag,OrderParamFile.c_str());

				OrderParamFile = WorkFolder + "/Obj_CorrMag.odat";
				ImportObject(CorrMag,OrderParamFile.c_str());

				OrderParamFile = WorkFolder + "/Obj_CorrStagMag.odat";
				ImportObject(CorrStagMag,OrderParamFile.c_str());
			}
			else
			{
				OrderParamFile = WorkFolder + "/Obj_MeanCoreMag.odat";
				ImportObject(MeanCoreMag,OrderParamFile.c_str());
			}
		}

		if(conf.SimType.compare("Part")!=0)
		{
			OrderParamFile = WorkFolder + "/Obj_MeanSubA.odat";
			ImportObject(MeanSubA,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanSubB.odat";
			ImportObject(MeanSubB,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanSubC.odat";
			ImportObject(MeanSubC,OrderParamFile.c_str());

			// Complex!
#ifdef COMPLEX_PARAMS
			OrderParamFile = WorkFolder + "/Obj_MeanComplexPhase.odat";
			ImportObject(MeanComplexPhase,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanSymmetryParameter.odat";
			ImportObject(MeanSymmetryParameter,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanRadius.odat";
			ImportObject(MeanRadius,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanComplexPhasePerLayer.odat";
			ImportObject(MeanComplexPhasePerLayer,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanSymmetryParameterPerLayer.odat";
			ImportObject(MeanSymmetryParameterPerLayer,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanRadiusPerLayer.odat";
			ImportObject(MeanRadiusPerLayer,OrderParamFile.c_str());
#endif
		}

		if(conf.SimType.compare("Moessner")==0)
		{
			OrderParamFile = WorkFolder + "/Obj_SpatialCorrN3N3.odat";
			ImportObject(SpatialCorrN3N3,OrderParamFile.c_str());

			OrderParamFile = WorkFolder + "/Obj_MeanLocalN3.odat";
			ImportObject(MeanLocalN3,OrderParamFile.c_str());
		}

		OrderParamFile = WorkFolder + "/Obj_MeanLocalNf.odat";
		ImportObject(MeanLocalNf,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_MeanLocalDimer.odat";
		ImportObject(MeanLocalDimer,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_CorrSzSz.odat";
		ImportObject(SzSzCorrelation,OrderParamFile.c_str());

		OrderParamFile = WorkFolder + "/Obj_CorrDimer.odat";
		ImportObject(DimerCorrelation,OrderParamFile.c_str());
	}
	// >>>> MC
	// ---> Preparations

	// Extra file names
	string 	outAccept 		= WorkFolder + "/Accept.dat";		// Accepance rates
	string 	outMeanClusterSize 	= WorkFolder + "/ClusterSize.dat";	// Mean cluster size
	string 	outMeanClusterNumber 	= WorkFolder + "/ClusterNumber.dat";	// Mean cluster number

	// Counters
	long 	NbUpdates 			= 0;		// Number of updates done
	int 	MeasureIndex 		= 0;
	int dummyIACLength = MeanMag.GetIACLength();
	long	UpdateInterval		= -1;
	if(parameters.MeasuresToTake>dummyIACLength*10)
	{
		UpdateInterval = parameters.MeasuresToTake/10;
	}
	else
	{
		UpdateInterval = parameters.MeasuresToTake;
	}

	int		TimesUpdated 		= 0;

	double 	FinalNumberOfClusters	= 0;
	double 	FinalMeanClusterSize 	= 0;

	// ---> Prepare
#if BORDER_TYPE == 4
// ANTI-PERIODIC BORDERS - SET WEIGHTS
	conf.ReadAntiWeights();
#endif
	conf.SetNeighbours();
	conf.SetSecondNeightbours();

	conf.SetPositions();
	if(conf.SimType.compare("Moessner")==0)
	{
		conf.SetDistances();
	}


	conf.SetFlipTables();
	conf.MC_Prepare_Mixed(parameters);

	// Some auxiliary variables
	double 		dummyParam = 0;

	vector<double>	dummyCorrMag(conf.nx + conf.p,0);
	vector<double>	dummyCorrStagMag(conf.nx + conf.p,0);
	vector<double>	dummyNf(4,-1);
	vector<double>	dummySublattice(3,-1);
	vector<double>	dummyMeanLocalNf(conf.L,0);
	vector<double>	dummyMeanKinEnergy(conf.L,0);
	vector<double> 	dummyMeanLocalDimer(conf.L*conf.NbOfNeights,0);
	vector<double> 	dummyMeanStar3(conf.L,0);
	vector<double> 	dummyCorr(conf.N,0);
	vector<double> 	dummyEnergies(2,0);
	vector<double>  dummyNewEnergies(2,0);
	vector<double> 	dummySiteMean(conf.L*conf.NbOfNeights,0);
	vector<double>  dummyPart(conf.nx*conf.ny,0);

	vector<double>  dummySpatialCorrelation(conf.L*conf.nbOfDists,0);
	vector<double>  dummyLocalDensities(conf.L,0);

#ifdef COMPLEX_PARAMS
	complex<double> dummyPhase = 0.;
	double			dummyRadius = 0.;

	vector<complex<double> > dummyPhasePerLayer(conf.N,0.);
	vector<double> dummyRadiusPerLayer(conf.N,0.);
#endif

	int 		OrderBunch;

	// ---> Set output and timestamp
	ofstream	output;
	int pid = getpid();

	char HostName[64];
	gethostname(HostName,64);

	string 		timestampFilename = WorkFolder + "/../status/status_" + HostName + "_" + ToString(conf.Vt);
	ofstream	timestampOutput(timestampFilename.c_str(), ios::out | ios::trunc);

	cout << timestampFilename << endl;

	time_t 		time_start = time(NULL);
	time_t 		time_end = 0;
	int 		time_interval = 0;
	int 		time_interval_sec = 0;
	int 		time_interval_min = 0;
	int			time_interval_hours = 0;
	int 		time_interval_days = 0;

	cout << "-->  V/t = " << conf.Vt << " @ " << HostName <<  ", PID " << pid << ", started at " << ctime(&time_start);

	timestampOutput << " * Simulation started at " << ctime(&time_start);
	timestampOutput.close();

	if(conf.SimType.compare("Manual")==0)
	{
		conf.SetSublatticeManual();
	}
	else if(conf.SimType.compare("Moessner")==0)
	{
		conf.SetSublatticeMoessner();
	}

	// ---> Run (you fools) !
	while(MeasureIndex < parameters.NbMeasures)
	{
		conf.ClusterIter_Mixed();
		NbUpdates 				+= conf.GetClusterSize();
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

        		        conf.GetAllQEnergy_Mixed(dummyEnergies);
				MeanQEnergy.AddData(dummyEnergies[0]+dummyEnergies[1]);
        		        MeanQKinEnergy.AddData(dummyEnergies[0]);
                		MeanQPotEnergy.AddData(dummyEnergies[1]);

				conf.GetAllNewQEnergy_Mixed(dummyNewEnergies);
				MeanQNewEnergy.AddData(dummyNewEnergies[0]+dummyNewEnergies[1]);				
				MeanQNewKinEnergy.AddData(dummyNewEnergies[0]);

                		conf.GetKinEnergyDensity(dummyMeanKinEnergy);
				MeanKinEnergyDensity.AddData(dummyMeanKinEnergy);

				if(conf.SimType.compare("Part")==0)
				{
					conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer);
				}
				else
				{

#ifdef COMPLEX_PARAMS
					conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer,dummySublattice,
									dummyPhasePerLayer,dummyRadiusPerLayer,dummyPhase,dummyRadius);
#else
					conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer,dummySublattice);
#endif

					MeanSubA.AddData(dummySublattice[0]);
					MeanSubB.AddData(dummySublattice[1]);
					MeanSubC.AddData(dummySublattice[2]);

#ifdef COMPLEX_PARAMS
					MeanSymmetryParameter.AddData(dummyRadius);
					MeanRadius.AddData(abs(dummyPhase));
					MeanComplexPhase.AddData(dummyPhase);

					for(int nnn = 0; nnn < conf.N; ++nnn)
					{
						MeanSymmetryParameterPerLayer.AddData(dummyRadiusPerLayer[nnn]);
						MeanRadiusPerLayer.AddData(abs(dummyPhasePerLayer[nnn]));
						MeanComplexPhasePerLayer.AddData(dummyPhasePerLayer[nnn]);
					}
#endif
				}

				MeanN0.AddData(dummyNf[0]);
				MeanN1.AddData(dummyNf[1]);
				MeanN2.AddData(dummyNf[2]);
				MeanN3.AddData(dummyNf[3]);

				MeanLocalNf.AddData(dummyMeanLocalNf);
				MeanLocalDimer.AddData(dummyMeanLocalDimer);

				conf.GetSzSzCorrelation(dummyCorr);
				SzSzCorrelation.AddCorr(dummyCorr);

				conf.GetDimerDimerCorrelation(dummyCorr);
				DimerCorrelation.AddCorr(dummyCorr);

				if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
				{
					if(conf.ny==1)
					{
						conf.BuildEquivalentSpinChain();

						dummyParam = conf.CalculateStagMag();
						MeanStagMag.AddData(dummyParam);

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

				if(conf.SimType.compare("Part")==0&&conf.initCondType < 2)
				{
					conf.CalculateMeanPartition();
					conf.GetMeanPart(dummyPart);
					MeanPart.AddData(dummyPart);
				}
			}

			++MeasureIndex;
			NbUpdates = 0;
		}

		// Temporary output
		if((MeasureIndex>=(TimesUpdated+1)*UpdateInterval+parameters.NbMeasures-parameters.MeasuresToTake)&&parameters.MeasuresToTake!=0)
		{
			++TimesUpdated;
			conf.RaiseMeasures(UpdateInterval);

			time_end = time(NULL);

			time_interval = difftime(time_end,time_start);

			ConvertTime(time_interval,time_interval_sec,time_interval_min,time_interval_hours,time_interval_days);

			timestampOutput.open(timestampFilename.c_str(), ios::out | ios::app);
			timestampOutput << " *         Update no. " << TimesUpdated << " : " << ctime(&time_end);
			timestampOutput << "                          " << time_interval_days << "d " << time_interval_hours << "h " << time_interval_min << "m " << time_interval_sec << "s" << endl;
			timestampOutput << "                          " << MeasureIndex << " measures for this run" << endl;
			timestampOutput << "                          " << conf.measurementsDone << " measures done in total" << endl;
			timestampOutput.close();

			// >>>> Post-QMC calculations
			// > Binning / bunching and autocorrelations
			output.open(outAccept.c_str(), ios::out | ios::trunc);
			output << conf.Vt << " " << conf.GetTotalAccept()/FinalNumberOfClusters << endl;
			output.close();

			output.open(outMeanClusterNumber.c_str(), ios::out | ios::trunc);
			output << conf.Vt << " " << FinalNumberOfClusters << endl;
			output.close();

			output.open(outMeanClusterSize.c_str(), ios::out | ios::trunc);
			output << conf.Vt << " " << FinalMeanClusterSize/FinalNumberOfClusters << endl;
			output.close();

			if(TimesUpdated*UpdateInterval > dummyIACLength)
			{
				MeanMag.CalculateError(OrderBunch);
				MeanMag.PrintError();

				MeanEnergy.CalculateError(OrderBunch);
				MeanEnergy.PrintError();

				MeanQEnergy.CalculateError(OrderBunch);
				MeanQEnergy.PrintError();
				
				MeanQKinEnergy.CalculateError(OrderBunch);
				MeanQKinEnergy.PrintError();

				MeanQPotEnergy.CalculateError(OrderBunch);
				MeanQPotEnergy.PrintError();

				MeanQNewEnergy.CalculateError(OrderBunch);
				MeanQNewEnergy.PrintError();

				MeanQNewKinEnergy.CalculateError(OrderBunch);
				MeanQNewKinEnergy.PrintError();

				MeanN0.CalculateError(OrderBunch);
				MeanN0.PrintError();

				MeanN1.CalculateError(OrderBunch);
				MeanN1.PrintError();

				MeanN2.CalculateError(OrderBunch);
				MeanN2.PrintError();

				MeanN3.CalculateError(OrderBunch);
				MeanN3.PrintError();

				if(conf.SimType.compare("Part")!=0)
				{
					MeanSubA.CalculateError(OrderBunch);
					MeanSubA.PrintError();
					MeanSubB.CalculateError(OrderBunch);
					MeanSubB.PrintError();
					MeanSubC.CalculateError(OrderBunch);
					MeanSubC.PrintError();
				}

				if(conf.SimType.compare("Part")==0&&conf.initCondType == 0)
				{
					if(conf.ny==1)
					{
						MeanStagMag.CalculateError(OrderBunch);
						MeanStagMag.PrintError();
					}
					else
					{
						MeanCoreMag.CalculateRMSMean();
						MeanCoreMag.CalculateError(OrderBunch);
						MeanCoreMag.PrintError();
					}
				}
			}

			// >>>> Export data
			MeanMag.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanMag.odat";
			ExportObject(MeanMag,OrderParamFile.c_str());

			MeanEnergy.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanEnergy.odat";
			ExportObject(MeanEnergy,OrderParamFile.c_str());

			MeanQEnergy.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanQEnergy.odat";
			ExportObject(MeanQEnergy,OrderParamFile.c_str());

                        MeanQNewEnergy.RaiseMeasures(UpdateInterval);
                        OrderParamFile = WorkFolder + "/Obj_MeanQNewEnergy.odat";
                        ExportObject(MeanQNewEnergy,OrderParamFile.c_str());

			MeanQPotEnergy.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanPotEnergy.odat";
			ExportObject(MeanQPotEnergy,OrderParamFile.c_str());

			MeanQKinEnergy.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanKinEnergy.odat";
			ExportObject(MeanQKinEnergy,OrderParamFile.c_str());

			MeanQNewKinEnergy.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanNewKinEnergy.odat";
			ExportObject(MeanQNewKinEnergy,OrderParamFile.c_str());

			MeanKinEnergyDensity.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanKinEnergyDensity.odat";
			ExportObject(MeanKinEnergyDensity,OrderParamFile.c_str());

			MeanN0.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanN0.odat";
			ExportObject(MeanN0,OrderParamFile.c_str());

			MeanN1.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanN1.odat";
			ExportObject(MeanN1,OrderParamFile.c_str());

			MeanN2.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanN2.odat";
			ExportObject(MeanN2,OrderParamFile.c_str());

			MeanN3.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanN3.odat";
			ExportObject(MeanN3,OrderParamFile.c_str());

			if(conf.SimType.compare("Moessner")==0)
			{
				SpatialCorrN3N3.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder +  "/Obj_SpatialCorrN3N3.odat";
				ExportObject(SpatialCorrN3N3,OrderParamFile.c_str());

				MeanLocalN3.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder +  "/Obj_MeanLocalN3.odat";
				ExportObject(MeanLocalN3,OrderParamFile.c_str());

#ifdef COMPLEX_PARAMS
				MeanSymmetryParameter.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanSymmetryParameter.odat";
				ExportObject(MeanSymmetryParameter,OrderParamFile.c_str());

				MeanRadius.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanRadius.odat";
				ExportObject(MeanRadius,OrderParamFile.c_str());

				MeanComplexPhase.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder +  "/Obj_MeanComplexPhase.odat";
				ExportObject(MeanComplexPhase,OrderParamFile.c_str());

				MeanSymmetryParameterPerLayer.RaiseMeasures(conf.N*UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanSymmetryParameterPerLayer.odat";
				ExportObject(MeanSymmetryParameterPerLayer,OrderParamFile.c_str());

				MeanRadiusPerLayer.RaiseMeasures(conf.N*UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanRadiusPerLayer.odat";
				ExportObject(MeanRadiusPerLayer,OrderParamFile.c_str());

				MeanComplexPhasePerLayer.RaiseMeasures(conf.N*UpdateInterval);
				OrderParamFile = WorkFolder +  "/Obj_MeanComplexPhasePerLayer.odat";
				ExportObject(MeanComplexPhasePerLayer,OrderParamFile.c_str());
#endif
			}

			if(conf.SimType.compare("Part")!=0)
			{
				MeanSubA.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanSubA.odat";
				ExportObject(MeanSubA,OrderParamFile.c_str());
				MeanSubB.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanSubB.odat";
				ExportObject(MeanSubB,OrderParamFile.c_str());
				MeanSubC.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanSubC.odat";
				ExportObject(MeanSubC,OrderParamFile.c_str());
			}

			MeanLocalDimer.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanLocalDimer.odat";
			ExportObject(MeanLocalDimer,OrderParamFile.c_str());

			MeanLocalNf.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_MeanLocalNf.odat";
			ExportObject(MeanLocalNf,OrderParamFile.c_str());

			SzSzCorrelation.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_CorrSzSz.odat";
			ExportObject(SzSzCorrelation,OrderParamFile.c_str());

			DimerCorrelation.RaiseMeasures(UpdateInterval);
			OrderParamFile = WorkFolder + "/Obj_CorrDimer.odat";
			ExportObject(DimerCorrelation,OrderParamFile.c_str());

			if(conf.SimType.compare("Part")==0&&(conf.initCondType == 0 || conf.initCondType == 1))
			{
				if(conf.ny==1)
				{
					MeanStagMag.RaiseMeasures(UpdateInterval);
					OrderParamFile = WorkFolder + "/Obj_MeanStagMag.odat";
					ExportObject(MeanStagMag,OrderParamFile.c_str());

					CorrMag.RaiseMeasures(UpdateInterval);
					OrderParamFile = WorkFolder + "/Obj_CorrMag.odat";
					ExportObject(CorrMag,OrderParamFile.c_str());

					CorrStagMag.RaiseMeasures(UpdateInterval);
					OrderParamFile = WorkFolder + "/Obj_CorrStagMag.odat";
					ExportObject(CorrStagMag,OrderParamFile.c_str());
				}
				else
				{
					MeanCoreMag.RaiseMeasures(UpdateInterval);
					OrderParamFile = WorkFolder + "/Obj_MeanCoreMag.odat";
					ExportObject(MeanCoreMag,OrderParamFile.c_str());
				}

				MeanPart.RaiseMeasures(UpdateInterval);
				OrderParamFile = WorkFolder + "/Obj_MeanPart.odat";
				ExportObject(MeanPart,OrderParamFile.c_str());
			}

			ExportObject(conf,outConf.c_str());
			conf.ExportRNG(RNGengineFile.c_str());
		}
	}

	if(parameters.MeasuresToTake==0)
	{
		timestampOutput.open(timestampFilename.c_str(), ios::out | ios::app);
		timestampOutput << " * No measurements were saved " << endl;
		timestampOutput.close();

		ExportObject(conf,outConf.c_str());
		conf.ExportRNG(RNGengineFile.c_str());
	}

	time_end = time(NULL);

	time_interval = difftime(time_end,time_start);

	ConvertTime(time_interval,time_interval_sec,time_interval_min,time_interval_hours,time_interval_days);

	cout << "     Total duration : " << time_interval_days << "d " << time_interval_hours << "h " << time_interval_min << "m " << time_interval_sec << "s" << endl;
	return 0;
}
