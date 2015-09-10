#include "main.h"
/*
 *		Program used to generate the 2D data files : local dimer density
 *      - Local dimer density for each plaquette
 *		- Local dimer density for each edge
 *		- Local kinetic energy density
 *		- Local potential energy density - V0 and V3
 */

void CpFile(const char inFile[], const char outFile[])
{
	ofstream  	dst(outFile);
	ifstream  	src(inFile);

	dst << src.rdbuf();

	dst.close();
	src.close();
};

// Test if a file exists
inline bool FileExists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
};

// Export a data file
void ExportTwoDData(bool FileTypeExists, string & fileType, const string & inputFolder, const string & outputFolder, input_params & inputPars, TwoDOrderParam & param, string & posTable)
{
	if(FileTypeExists)
	{
		// Set output/input files
		string outputFile = outputFolder + "/" + fileType + ".dat";
		string inputFile = inputFolder + "/Obj_" + fileType + ".odat";

		ImportObject(param,inputFile.c_str());

		// Set up lattice positions
		if(!inputPars.runType.compare("Manual"))
		{
			param.ReadPositions(posTable);
		}
		else
		{
			param.SetPositions(inputPars);
		}
			param.Export2D_Data(outputFile);
	}
};

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// ---> Initialize the program
	string foldername;

	// Get folder
	if(argc > 1)
	{
		foldername = argv[1];
	}
	else
	{
		cout << "Wrong input parameters" << endl;
		cout << "---> Need at least a folder name as an argument!!!" << endl;
		return 0;
	}

	// Get input parameters object
	input_params inputPars;
	string inputFilename = foldername + "/Obj_InputParams.odat";
	ImportObject(inputPars,inputFilename.c_str());

	// Get the positions table
	string posTable;

	/* If the simulation was done defining manually the lattice, a second argument,
	   with the name of the file containing the lattice's coordinates, must be given */
	if(!inputPars.runType.compare("Manual"))
	{
		if(argc > 2)
		{
			posTable = argv[2];
		}
		else
		{
			cout << "Manually defined lattice not found" << endl;
			cout << "---> Need the file with the lattice coordinates as a second argument!!!" << endl;
			return 0;
		}
	}

	// ---> Build the folders and test for the file existence

	// Test which 2D parameters exist
	string testFolder = foldername + "/Vt_" + ToString(inputPars.VtValues[0]);
	string testFile;

	/* Local dimer density on each plaquette */
	testFile =  testFolder + "/Obj_MeanLocalNf.odat";
	bool NfDensityBool = FileExists(testFile);
	if(NfDensityBool)
	{
		cout << " - Local dimer density per plaquette : OK " << endl;
	}
	else
	{
		cout << " - Local dimer density per plaquette : NOT FOUND " << endl;
	}

	/* Local dimer density on each edge */
	testFile =  testFolder + "/Obj_MeanLocalDimer.odat";
	bool DimerDensityBool = FileExists(testFile);
	if(DimerDensityBool)
	{
		cout << " - Local dimer density per edge      : OK " << endl;
	}
	else
	{
		cout << " - Local dimer density per edge      : NOT FOUND " << endl;
	}

	/* Kinetic energy density */
	testFile =  testFolder + "/Obj_MeanKinEnergyDensity.odat";
	bool KinDensityBool = FileExists(testFile);
	if(KinDensityBool)
	{
		cout << " - Kinetic energy density            : OK " << endl;
	}
	else
	{
		cout << " - Kinetic energy density            : NOT FOUND " << endl;
	}

	/* Potential energy density - V3 */
	testFile =  testFolder + "/Obj_MeanPotV3EnergyDensity.odat";
	bool PotV3DensityBool = FileExists(testFile);
	if(PotV3DensityBool)
	{
		cout << " - Potential V3 energy density       : OK " << endl;
	}
	else
	{
		cout << " - Potential V3 energy density       : NOT FOUND " << endl;
	}

	/* Potential energy density - V0 */
	testFile =  testFolder + "/Obj_MeanPotV0EnergyDensity.odat";
	bool PotV0DensityBool = FileExists(testFile);
	if(PotV0DensityBool)
	{
		cout << " - Potential V0 energy density       : OK " << endl;
	}
	else
	{
		cout << " - Potential V0 energy density       : NOT FOUND " << endl;
	}

	// Set up output folders and export data
	string commonOutputFolder = foldername + "/2D_Data";
	mkdir(commonOutputFolder.c_str(),0777);

	string outputFolder;
	string inputFolder;
	string fileType;

	TwoDOrderParam param;

	for(int iii = 0; iii < inputPars.VtN; ++iii)
	{
		outputFolder = commonOutputFolder + "/Vt_" + ToString(inputPars.VtValues[iii]);
		inputFolder  = foldername + "/Vt_" + ToString(inputPars.VtValues[iii]);

		mkdir(outputFolder.c_str(),0777);

		/* Local dimer density on each plaquette */
		fileType = "MeanLocalNf";
		ExportTwoDData(NfDensityBool,fileType,inputFolder,outputFolder,inputPars,param,posTable);

		/* Local dimer density on each edge */
		fileType = "MeanLocalDimer";
		ExportTwoDData(DimerDensityBool,fileType,inputFolder,outputFolder,inputPars,param,posTable);

		/* Kinetic energy density */
		fileType = "MeanKinEnergyDensity";
		ExportTwoDData(KinDensityBool,fileType,inputFolder,outputFolder,inputPars,param,posTable);

		/* Potential energy density - V3 */
		fileType = "MeanPotV3EnergyDensity";
		ExportTwoDData(PotV3DensityBool,fileType,inputFolder,outputFolder,inputPars,param,posTable);

		/* Potential energy density - V0 */
		fileType = "MeanPotV0EnergyDensity";
		ExportTwoDData(PotV0DensityBool,fileType,inputFolder,outputFolder,inputPars,param,posTable);
	}

	return 0;
}
