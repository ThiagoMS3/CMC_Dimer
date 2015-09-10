#include "main.h"
/*
 *    Classical Monte Carlo simulation for the QDM - Local updates version
 *    >>> Based on Thomas' notes (14 May 2013)
 *
 */

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
	string inputString;
	string filename;
	string foldername;
	double Vt = 0;

	// Get file
	if(argc > 1)
	{
		inputString = argv[1];
	}
	else
	{
		cout << "Wrong parameters - need at least a filename" << endl;
		return 0;
	}

	if(argc > 2)
	{
		istringstream ss(argv[2]);
		if (!(ss >> Vt))
		{
			cout << "Invalid number " << argv[2] << '\n';
			return 0;
		}
	}


	filename 		= "./" + inputString + "/Obj_conf.odat";
	foldername 		= "./" + inputString;
	SysConf conf;

	ImportObject(conf,filename.c_str());

	conf.SetNeighbours();
	conf.SetFlipTables();

	vector<double> dummyNf(4,-1);
	vector<double>	dummyMeanLocalNf(conf.L);
	vector<double>	dummyMeanStar3(conf.L);
	vector<double> dummyMeanLocalDimer(conf.L*conf.NbOfNeights);

	TwoDOrderParam dummyLocalNf;
	TwoDOrderParam dummyLocalDimer;
	TwoDOrderParam dummyStar3;

	cout << " ----------------------------------------------------------------------------- " << endl << endl;
	cout << " ---> Config file : " << filename << endl << endl;
	cout << " ----------------------------------------------------------------------------- " << endl << endl;
	conf.PrintConfigStatus();
	cout << " ----------------------------------------------------------------------------- " << endl;
	cout << " ---> Exporting local Nf's : " <<endl;

	dummyLocalNf.Initialize(conf.L,1);
	dummyLocalNf.SetFolder(foldername.c_str());
	dummyLocalNf.SetType("MeanLocalNf");
	dummyLocalNf.RestartMean(1);

	dummyLocalDimer.Initialize(conf.L,conf.NbOfNeights);
	dummyLocalDimer.SetFolder(foldername.c_str());
	dummyLocalDimer.SetType("MeanLocalDimer");
	dummyLocalDimer.RestartMean(1);

	dummyStar3.Initialize(conf.L,1);
	dummyStar3.SetFolder(foldername.c_str());
	dummyStar3.SetType("Star3Lattice");
	dummyStar3.RestartMean(1);

	conf.GetSiteNf(dummyMeanLocalNf,dummyNf,dummyMeanLocalDimer);

	dummyLocalNf.AddData(dummyMeanLocalNf);
	ExportObject(dummyLocalNf,"Obj_MeanLocalNf.odat");

	dummyLocalDimer.AddData(dummyMeanLocalDimer);
	ExportObject(dummyLocalDimer,"Obj_MeanLocalDimer.odat");

	dummyStar3.AddData(dummyMeanStar3);
	ExportObject(dummyStar3,"Obj_Star3.odat");

	cout << " ----------------------------------------------------------------------------- " << endl;

	return 0;
}
