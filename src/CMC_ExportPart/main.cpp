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
	string outputString = "out_Part.dat";

	string filename;
	string confName;
	string posTable;

	// Get file
	if(argc > 1)
	{
		filename = argv[1];
	}
	else
	{
		cout << "Wrong parameters - need at least a filename" << endl;
		return 0;
	}

	if(argc > 2)
	{
		outputString = argv[2];
	}

	TwoDOrderParam param;
	ImportObject(param,filename.c_str());
	param.GetParent(confName);

	SysConf conf;
	ImportObject(conf,confName.c_str());

	param.Export2D_Matrix(outputString,conf.nx,conf.ny);

	return 0;
}
