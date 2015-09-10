#include "main.h"

/*
 *    Classical Monte Carlo simulation for the QDM - Local updates version
 *    >>> Based on Thomas' notes (14 May 2013)
 *
 */

//void MergeSuperlattice(vector<string> &VtFolder,string& WorkFolder, int totalVtN, int startVtN, input_params inputPars)
//{
//	string ObjFile = "/Obj_MeanLocalNf.odat";
//	string MeanFileName = WorkFolder + "/Final_SuperlatticeNf.dat";
//
//	string fileName;
//
//	TwoDOrderParam MeanObj;
//
//	vector<vector<double> > MeanOut(totalVtN,vector<double>(3,0));
//
//	ofstream output;
//
//	for(int iii = 0; iii < totalVtN; ++iii)
//	{
//		fileName = VtFolder[iii] + ObjFile;
//		ImportObject(MeanObj,fileName.c_str());
//
//		MeanObj.NewGetSuperlatticeNf(MeanOut[iii],inputPars.nx,inputPars.ny);
//	}
//
//	output.open(MeanFileName.c_str(),ios::trunc);
//	for(int iii = 0; iii < totalVtN; ++iii)
//	{
//		output << inputPars.VtValues[iii + startVtN] << " " << MeanOut[iii][0] << " " << MeanOut[iii][1] << " " << MeanOut[iii][2]<< endl;
//	}
//	output.close();
//};


void MergeSuperlattice(vector<string> &VtFolder,string& WorkFolder, int totalVtN, int startVtN, input_params inputPars)
{
	string ObjFile = "/Obj_MeanLocalNf.odat";
	string MeanFileName = WorkFolder + "/Final_SuperlatticeNf.dat";

	string fileName;

	TwoDOrderParam MeanObj;

	vector<vector<double> > MeanOut(totalVtN,vector<double>(3,0));

	ofstream output;

	cout << MeanFileName << endl;
	for(int iii = 0; iii < totalVtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanObj.GetSuperlatticeNf(MeanOut[iii],inputPars.nx,inputPars.ny);
	}

	output.open(MeanFileName.c_str(),ios::trunc);
	for(int iii = 0; iii < totalVtN; ++iii)
	{
		output << inputPars.VtValues[iii + startVtN] << " " << MeanOut[iii][0] << " " << MeanOut[iii][1] << " " << MeanOut[iii][2]<< endl;
	}
	output.close();
};

void RectManualMergeSuperlattice(vector<string> &VtFolder,string& WorkFolder, int totalVtN, int startVtN, input_params inputPars, vector<int>& neighTable)
{
	string ObjFile = "/Obj_MeanLocalNf.odat";
	string MeanFileName = WorkFolder + "/Final_SuperlatticeNf.dat";

	string fileName;

	TwoDOrderParam MeanObj;

	vector<vector<double> > MeanOut(totalVtN,vector<double>(3,0));

	ofstream output;

	for(int iii = 0; iii < totalVtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanObj.RectManualGetSuperlatticeNf(MeanOut[iii],inputPars.nx,inputPars.ny,neighTable);
	}

	output.open(MeanFileName.c_str(),ios::trunc);
	for(int iii = 0; iii < totalVtN; ++iii)
	{
		output << inputPars.VtValues[iii + startVtN] << " " << MeanOut[iii][0] << " " << MeanOut[iii][1] << " " << MeanOut[iii][2]<< endl;
	}
	output.close();
};

void RhombusManualMergeSuperlattice(vector<string> &VtFolder,string& WorkFolder, int totalVtN, int startVtN, input_params inputPars, vector<int>& neighTable)
{
	string ObjFile = "/Obj_MeanLocalNf.odat";
	string MeanFileName = WorkFolder + "/Final_SuperlatticeNf.dat";

	string fileName;

	TwoDOrderParam MeanObj;

	vector<vector<double> > MeanOut(totalVtN,vector<double>(3,0));
	vector<int> SubSize(3,0);

	ofstream output;

	for(int iii = 0; iii < totalVtN; ++iii)
	{
		fileName = VtFolder[iii] + ObjFile;
		ImportObject(MeanObj,fileName.c_str());

		MeanObj.RhombusManualGetSuperlatticeNf(MeanOut[iii],SubSize,inputPars.nx,inputPars.ny,neighTable);
	}

	cout << "Sub-lattice sizes:" << endl;
	cout << "   A: " << SubSize[0] << "/" << inputPars.nx*inputPars.ny << endl;
	cout << "   B: " << SubSize[1] << "/" << inputPars.nx*inputPars.ny << endl;
	cout << "   C: " << SubSize[2] << "/" << inputPars.nx*inputPars.ny << endl;

	output.open(MeanFileName.c_str(),ios::trunc);
	for(int iii = 0; iii < totalVtN; ++iii)
	{
		output << inputPars.VtValues[iii + startVtN] << " " << MeanOut[iii][0] << " " << MeanOut[iii][1] << " " << MeanOut[iii][2]<< endl;
	}
	output.close();
};

void GetNeighbours(string& neighFile, int nx, int ny, vector<int>& neighTable)
{
	ifstream input(neighFile.c_str());

	int L = nx*ny;
	for(int iii = 0; iii < L; ++iii)
	{
		for(int jjj = 0; jjj < 6; ++jjj)
		{
			input >> neighTable[idxConv(6,iii,jjj)];
		}
		Jump(input,1);
	}

	input.close();
};

void HelpMsg()
{
	cout << "Usage" << endl;
	cout << "    GetSuperlattice [options] [input parameters file] [output folder]" << endl;
	cout << "Options" << endl;
	cout << "    -h    Displays this help message" << endl << endl;
	cout << "    -T [type]" << endl;
	cout << "          Lattice type. Default is 'Periodic'" << endl;
	cout << "             Periodic  -> will suppose periodic boundary conditions" << endl;
	cout << "             Rectangle -> will suppose open boundary conditions on a rectangular lattice" << endl;
	cout << "             Rhombus   -> will suppose open boundary conditions on a rhombus lattice" << endl;
	cout << "    -N [file]" << endl;
	cout << "          Will read the neighbour relations from [file]" << endl;
};

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	int count;

	// Flags
	int typeFlag 	= 0;   // 0 = Periodic, 1 = Rect., 2 = Rhombus
	int neighFlag 	= 0;

	string neighFile;
	string typeString;

	while ((count = getopt(argc,argv,"hT:N:"))!=-1)
	{
		switch(count)
		{
		case 'h':
			HelpMsg();
			return 0;
			break;
		case 'T':
			typeString = optarg;
			if(!typeString.compare("Periodic"))
			{
				typeFlag = 0; // 0 = Periodic
			}
			else if(!typeString.compare("Rectangle"))
			{
				typeFlag = 1; // 1 = Rect.
			}
			else if(!typeString.compare("Rhombus"))
			{
				typeFlag = 2; // 2 = Rhombus
			}
			break;
		case 'N':
			neighFlag = 1;
			neighFile = optarg;

			break;
		}
	}

	string inputFile;
	string WorkFolder;

	if(argc!=optind + 2)
	{
		cout << "Error: Not enough arguments!" << endl;
		HelpMsg();
		return 0;
	}
	else
	{
		inputFile = argv[optind];
		WorkFolder = argv[optind+1];
	}

	// ---> Read config files
	input_params inputPars;
	ImportObject(inputPars,inputFile.c_str());

	vector<int> neighTable(6*inputPars.nx*inputPars.ny,-1);

	if(neighFlag)
	{
		GetNeighbours(neighFile,inputPars.nx,inputPars.ny,neighTable);
	}

	cout << "got it" << endl;
	// >>>> Set up file names
	// ---> Folders
	vector<string> 		  vtFolder(inputPars.VtN,WorkFolder + "/Vt_");
	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		vtFolder[iii] = vtFolder[iii] + ToString(inputPars.VtValues[iii]);
	}

	switch(typeFlag)
	{
	case 0:
		MergeSuperlattice(vtFolder,WorkFolder,inputPars.VtN,0,inputPars);
		break;
	case 1:
		RectManualMergeSuperlattice(vtFolder,WorkFolder,inputPars.VtN,0,inputPars,neighTable);
		break;
	case 2:
		RhombusManualMergeSuperlattice(vtFolder,WorkFolder,inputPars.VtN,0,inputPars,neighTable);
		break;
	}

	return 0;
}
