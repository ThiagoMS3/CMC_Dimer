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
};

bool TestFolder(string& folderName)
{
	struct stat sb;
	return stat(folderName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
};

void UsageHelp()
{
	cout << "--> Update Partition simulations program" << endl << endl;

	cout << "    Usage: ./Update_Part_Lab [folder] [input file]" << endl << endl;

	cout << "  > Obligatory arguments" << endl;
	cout << "       [folder]:      path to the simulation folder to be updated" << endl << endl;

	cout << "  > Optional arguments" << endl;
	cout << "       [input file]:  path to the file containing the new parameters" << endl;
	cout << "                      Default value = \"params_Part.in\"" << endl << endl;
};

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// ---> Set parameters
	// >>>> Input
	input_params input;
	int pid;
	string filename;		// Input file
	string outFolder;		// Folder of the simulation to be updated

	// Get file and folder
	if(argc == 1)
	{
		cout << "--------------------------------------------------------------------------------" << endl << endl;
		UsageHelp();
		cout << "--------------------------------------------------------------------------------" << endl;
		return 0;
	}

	if(argc == 2)
	{
		outFolder = argv[1];
		filename = "params_Part.in";
	}

	if(argc == 3)
	{
		outFolder = argv[1];
		filename = argv[2];
	}

	if(argc > 3)
	{
		cout << "--------------------------------------------------------------------------------" << endl << endl;

		cout << "ERROR: more than two arguments given, simulation WON'T be updated!" << endl << endl;
		UsageHelp();
		cout << "--------------------------------------------------------------------------------" << endl;
		return 0;
	}

	if(!TestFolder(outFolder))
	{
		cout << " ERROR: folder \"" << outFolder << "\" not found!" << endl;
		return 0;
	}

	// Read input
	ifstream inputF(filename.c_str());

	if(inputF.good())
	{
		pid = getpid();
		input.SetParamPart(inputF);
		inputF.close();
	}
	else
	{
		cout << " ERROR: input file \"" << filename << "\" not found!" << endl;
		inputF.close();
		return 0;
	}

	// ---> Geometry parameters
	string SimType = "Part";

	// ---> Prepare the output file/folder names
	// >>>> V/t folders
	vector<string> 		  vtFolder(input.VtN,outFolder + "Vt_");
	for(int iii = 0; iii < input.VtN; ++iii)
	{
		vtFolder[iii] = vtFolder[iii] + ToString(input.VtValues[iii]);
		mkdir(vtFolder[iii].c_str(),0777);
	}

	// >>>> Filenames
	string dummyString;
	string dummyStringBis;
	string outBase 		= "./" + outFolder + "/";
	string scriptFilename	= "./run_Part_" + ToString(pid) + ".sh";
	string commandsFilename	= "comm_Part_" + ToString(pid) + ".sh";
	string sysConfFilename;
	string 	param_out 	= outBase + filename;
	CpFile(filename.c_str(),param_out.c_str());

	// ---> MC parameters
//	int L = input.nx*input.ny;
//	MCParameters params(input,L,outFolder.c_str());
	MCParameters params;

	// ---> Set up files
	// >>>> For the geometry

	dummyString = "./" + outFolder + "/Obj_InputParams.odat";
	ExportObject(input,dummyString.c_str());

	for(int iii=0; iii < input.VtN; ++iii)
	{
		// >>>> For the MC parameters
		dummyString 	= "./" + vtFolder[iii] + "/Obj_MCparams.odat";

		ImportObject(params,dummyString.c_str());

		params.UpdateData(input);

		ExportObject(params,dummyString.c_str());
	}


	if(input.compType==0)
	{
		ofstream commandsFile(commandsFilename.c_str(),ios::trunc);
		int procCounter = 0;
		for(int iii = 0; iii < input.VtN; ++iii)
		{
			commandsFile 	<< "nice -19 ./" << input.execName << " " << vtFolder[iii]
							<< "/Obj_conf.odat" << " " << vtFolder[iii]
							<< " > ${PWD}/" << vtFolder[iii] << "/out_log.txt &" <<  endl;
			++procCounter;
			if(procCounter==input.numProc)
			{
				commandsFile << "wait" << endl;
				procCounter = 0;
			}
		}
		commandsFile.close();
	}
	else if(input.compType==1||input.compType==2)
	{
		ofstream commandsFile(commandsFilename.c_str(),ios::trunc);
		for(int iii = 0; iii < input.VtN; ++iii)
		{
			commandsFile 	<< "nice -19 ./" << input.execName << " " << vtFolder[iii]
							<< "/Obj_conf.odat" << " " << vtFolder[iii] <<  endl;
		}
		commandsFile.close();
	}
	else if(input.compType==3)
	{
		ofstream commandsFile(commandsFilename.c_str(),ios::trunc);
		for(int iii = 0; iii < input.VtN; ++iii)
		{
			commandsFile 	<< "srun --exclusive -N1 -n1 nice -19 ./" << input.execName << " " << vtFolder[iii]
							<< "/Obj_conf.odat" << " " << vtFolder[iii]  <<  endl;
		}
		commandsFile.close();
	}

	if(input.compType==0)
	{
		// Direct
		cout << "---> Output type        : direct - Lab (local)" << endl;
		cout << "     Execution script   : " << scriptFilename << endl;
		cout << "     Commands list      : " << commandsFilename << endl;

		ofstream scriptFile(scriptFilename.c_str(),ios::trunc);
		scriptFile << ". " << commandsFilename << endl;
		scriptFile << "./Merge " << outFolder << "/Obj_InputParams.odat" << " " << outFolder <<  endl;
		scriptFile.close();
	}
	if(input.compType==1)
	{
		// Parallel - simple
		cout << "---> Output type        : parallel - Lab (local)" << endl;
		cout << "     Execution script   : " << scriptFilename << endl;
		cout << "     Commands list      : " << commandsFilename << endl;

		ofstream scriptFile(scriptFilename.c_str(),ios::trunc);
		scriptFile << "parallel -j " << input.numProc << " < " << commandsFilename << endl;
		scriptFile << "./Merge " << outFolder << "/Obj_InputParams.odat" << " " << outFolder <<  endl;
		scriptFile.close();
	}
	if(input.compType==2)
	{
		// Torque - parallel
		string mergeScript 	= "merge_script_"  + ToString(pid) + ".sh";

		cout << "---> Output type        : TORQUE - Totoro, MESU" << endl;
		cout << "     Execution script   : none" << endl;
		cout << "     Commands list      : " << commandsFilename << endl;
		cout << "     Merge script       : " << mergeScript << endl;

		ofstream mergeFile(mergeScript.c_str(),ios::trunc);
		mergeFile << "./Merge " << outFolder << "/Obj_InputParams.odat" << " " << outFolder <<  endl;
		mergeFile.close();
	}
	if(input.compType==3)
	{
		// SLURM - parallel
		string SLURMbase 	= "parallel_base.sh";
		string mergeScript 	= "merge_script_"  + ToString(pid) + ".sh";

		cout << "---> Output type        : SLURM - Lab (cluster)" << endl;
		cout << "     Execution script   : " << scriptFilename << endl;
		cout << "     Commands list      : " << commandsFilename << endl;
		cout << "     Base file          : " << SLURMbase << endl;
		cout << "     Merge script       : " << mergeScript << endl << endl;
		cout << "--------------------------------------------------------------------------------" << endl << endl;
		cout << "---> Open '" << scriptFilename << "' to edit :" << endl;
		cout << "        Job name        : #SBATCH --job-name=[name]" << endl;
		cout << "        Number of CPU's : #SBATCH --ntask=[number]" <<  endl << endl;
		cout << "--------------------------------------------------------------------------------" << endl << endl;
		cout << "---> Submit job with       : 'sbatch " << scriptFilename << "'" << endl ;
		cout << "                             (from any computer, no need to use 'screen')" << endl << endl;
		cout << "     View job queue        : 'squeue'" << endl << endl;
		cout << "     View your jobs        : 'squeue -u [your username]'" << endl << endl;
		cout << "     To kill a job         : > get JOB ID no. from 'squeue -u [your username]'" << endl;
		cout << "                             > then 'scancel [job ID]'" << endl << endl;
		cout << "     View cluster status   : 'sinfo'" << endl;
		cout << "                             > 'idle'       computers are 100% free" << endl;
		cout << "                             > 'mix'        computers are partially free" << endl;
		cout << "                             > 'alloc'      computers are 100% full" << endl;
		cout << "                             > 'drain/down' computers are offline" << endl << endl;
		cout << "     The 'status' sub-folder contains files with the progress of the simulation" << endl << endl;
		cout << "--------------------------------------------------------------------------------" << endl << endl;

		CpFile(SLURMbase.c_str(),scriptFilename.c_str());
		ofstream scriptFile(scriptFilename.c_str(),ios::app);
		ofstream mergeFile(mergeScript.c_str(),ios::trunc);

		scriptFile << "$parallel < " << commandsFilename << endl;
		mergeFile << "./Merge " << outFolder << "/Obj_InputParams.odat" << " " << outFolder <<  endl;

		scriptFile.close();
		mergeFile.close();
	}

	return 0;
}
