#include "main.h"
/*
 *    Classical Monte Carlo simulation for the QDM - Local updates version
 *    >>> Based on Thomas' notes (14 May 2013)
 *
 */

std::string get_working_path()
{
   char temp[256];
   return ( getcwd(temp, 256) ? std::string( temp ) : std::string("") );
}

void CpFile(const char inFile[], const char outFile[])
{
	ofstream  	dst(outFile);
	ifstream  	src(inFile);

	dst << src.rdbuf();

	dst.close();
	src.close();
}

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// ---> Set parameters
	// >>>> Input
	input_params input;

	string outFolder;
	string param_file;
	string rng_filename;

	int rng_steps;
	int rng_seed;

	cout << get_working_path() << endl;

	// Get file
	if(argc > 2)
	{
		outFolder = argv[1];

		istringstream ss(argv[2]);
		if (!(ss >> rng_steps))
		{
			cout << "Invalid number " << argv[2] << '\n';
			return 0;
		}
	}
	else
	{
		cout << "Not enough parameters - need at least a filename and a number of random steps" << endl;
		return 0;
	}

	string filename 	= outFolder + "Obj_conf.odat";
	string backup_name 	= outFolder + "OLD_Obj_conf.odat";

	SysConf conf;

	ImportObject(conf,filename.c_str());

	if(argc > 3)
	{
		istringstream ss(argv[3]);
		if (!(ss >> rng_seed))
		{
			cout << "Invalid number " << argv[3] << '\n';
			return 0;
		}
		conf.SetRNG(rng_seed);
	}

	CpFile(filename.c_str(),backup_name.c_str());

	// Set observables
	string OrderParamFile;

	OrderParam	MeanN0;
	OrderParam	MeanN1;
	OrderParam	MeanN2;
	OrderParam	MeanN3;

	TwoDOrderParam	MeanLocalDimer;
	TwoDOrderParam	MeanLocalNf;

	vector<double> dummyNf(4,-1);

	vector<double>	dummyMeanLocalNf(conf.L);
	vector<double> dummyMeanLocalDimer(conf.L*conf.NbOfNeights);


//	OrderParamFile = outFolder + "Obj_MeanN0.odat";
//	ImportObject(MeanN0,OrderParamFile.c_str());
//	MeanN0.RestartMean(rng_steps);
//
//	OrderParamFile = outFolder + "Obj_MeanN1.odat";
//	ImportObject(MeanN1,OrderParamFile.c_str());
//	MeanN1.RestartMean(rng_steps);
//
//	OrderParamFile = outFolder + "Obj_MeanN2.odat";
//	ImportObject(MeanN2,OrderParamFile.c_str());
//	MeanN2.RestartMean(rng_steps);
//
//	OrderParamFile = outFolder + "Obj_MeanN3.odat";
//	ImportObject(MeanN3,OrderParamFile.c_str());
//	MeanN3.RestartMean(rng_steps);

	OrderParamFile = outFolder + "Obj_MeanLocalNf.odat";
	ImportObject(MeanLocalNf,OrderParamFile.c_str());
	MeanLocalNf.RestartMean(rng_steps);

	OrderParamFile = outFolder + "Obj_MeanLocalDimer.odat";
	ImportObject(MeanLocalDimer,OrderParamFile.c_str());
	MeanLocalDimer.RestartMean(rng_steps);


	// Randomize !!!
//	int counter = 0;

	conf.SetNeighbours();
	conf.SetFlipTables();
	conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer);

	cout << " Initial n_i's : " << dummyNf[0] << " " << dummyNf[1] << " " << dummyNf[2] << " " << dummyNf[3] << endl;

	conf.RandomSpinFlips(rng_steps);


//	conf.SetRandomColFlips();
//	for(int iii = 0; iii < rng_steps; ++iii)
//	{
//		conf.RandomColFlip();
//	}

//	conf.initCondType = 1;
	cout << " Ended flips" << endl;
	conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer);
	MeanLocalNf.AddData(dummyMeanLocalNf);
	MeanLocalDimer.AddData(dummyMeanLocalDimer);
	cout << " Final n_i's : " << dummyNf[0] << " " << dummyNf[1] << " " << dummyNf[2] << " " << dummyNf[3] << endl;
	ExportObject(conf,filename.c_str());

//	int OrderBunch = 0;
//
//	MeanN0.CalculateError(OrderBunch);
//	MeanN0.PrintError();
//
//	MeanN1.CalculateError(OrderBunch);
//	MeanN1.PrintError();
//
//	MeanN2.CalculateError(OrderBunch);
//	MeanN2.PrintError();
//
//	MeanN3.CalculateError(OrderBunch);
//	MeanN3.PrintError();


//	OrderParamFile = outFolder + "Obj_MeanN0.odat";
//	ExportObject(MeanN0,OrderParamFile.c_str());
//
//	OrderParamFile = outFolder + "Obj_MeanN1.odat";
//	ExportObject(MeanN1,OrderParamFile.c_str());
//
//	OrderParamFile = outFolder + "Obj_MeanN2.odat";
//	ExportObject(MeanN2,OrderParamFile.c_str());
//
//	OrderParamFile = outFolder + "Obj_MeanN3.odat";
//	ExportObject(MeanN3,OrderParamFile.c_str());
//
	OrderParamFile = outFolder + "Obj_RNG_MeanLocalDimer.odat";
	ExportObject(MeanLocalDimer,OrderParamFile.c_str());

	OrderParamFile = outFolder + "Obj_RNG_MeanLocalNf.odat";
	ExportObject(MeanLocalNf,OrderParamFile.c_str());
//
//	cout << "Got " << counter << " n_0's" << endl;

	return 0;
}
