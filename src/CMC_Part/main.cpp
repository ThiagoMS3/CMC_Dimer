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

void UsageHelp()
{
	cout << "--> Set Partition simulations program" << endl << endl;

	cout << "    Usage: ./Set_Part_Lab [input file]" << endl << endl;

	cout << "  > Optional arguments" << endl;
	cout << "       [input file]:  path to the file containing the parameters" << endl;
	cout << "                      Default value = \"params_Part.in\"" << endl << endl;
}

// --- Main code ---
// -----------------
int main(int argc, char **argv)
{
	// ---> Set parameters
	// >>>> Input
	input_params input;
	int pid;
	string filename;		// Input file

	// Get file and folder
	if(argc == 1)
	{
		filename = "params_Part.in";
	}

	if(argc == 2)
	{
		filename = argv[1];
	}

	if(argc > 2)
	{
		cout << "--------------------------------------------------------------------------------" << endl << endl;

		cout << "ERROR: more than one argument given, simulation WON'T be set up!" << endl << endl;
		UsageHelp();
		cout << "--------------------------------------------------------------------------------" << endl;
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

	// >>>> Geometry parameters
	string SimType = "Part";
	SysConf conf(input,SimType);
	conf.SetNeighboursPart();
	if(conf.initCondType == 0)
	{
		conf.SetInitialPartAntiFerro();
	}
	else if(conf.initCondType == 1)
	{
		conf.SetInitialPartEmpty();
	}
	else if(conf.initCondType == 2)
	{
		conf.SetInitialPartStar3();
	}
	else
	{
		conf.SetInitialPart();
	}

	// >>>> Set up order parameters
	OrderParam dummyOrder(input);
	TwoDOrderParam dummy2DOrder;

	// ---> Prepare the output file/folder names
	// >>>> Output folder
	string outFolder	= input.outFolderName + "_" + ToString(input.nx) + "_" + ToString(input.ny) + "_" + ToString(input.p) + "_" + ToString(pid);
	mkdir(outFolder.c_str(),0777);

	// >>>> V/t folders
	vector<string> 		  vtFolder(input.VtN,outFolder + "/Vt_");
	for(int iii = 0; iii < input.VtN; ++iii)
	{
		vtFolder[iii] = vtFolder[iii] + ToString(input.VtValues[iii]);
		mkdir(vtFolder[iii].c_str(),0777);
	}

	string statusFolder = outFolder + "/status";
	mkdir(statusFolder.c_str(),0777);

	// >>>> Filenames
	string dummyString;
	string dummyStringBis;
	string outBase 		= "./" + outFolder + "/";
	string sysConfFilename	= outBase + "Obj_conf.odat";
	string 	param_out 	= outBase + filename;
	CpFile(filename.c_str(),param_out.c_str());

	// ---> MC parameters
	MCParameters params(input,conf.L,outFolder.c_str());

	// ---> Set up files
	// >>>> For the geometry
	ExportObject(conf,sysConfFilename.c_str());

	dummyString = "./" + outFolder + "/Obj_InputParams.odat";
	ExportObject(input,dummyString.c_str());

	for(int iii = 0; iii < input.VtN; ++iii)
	{
		// >>>> For the MC parameters
		dummyString 	= "./" + vtFolder[iii] + "/Obj_MCparams.odat";
		params.SetVt(input.VtValues[iii]);
		ExportObject(params,dummyString.c_str());

		// >>>> For the order parameters
		dummyString 	= "./" + vtFolder[iii];
		dummyOrder.SetVt(input.VtValues[iii]);
		dummyOrder.SetFolder(dummyString.c_str());

		dummy2DOrder.SetVt(input.VtValues[iii]);
		dummy2DOrder.SetFolder(dummyString.c_str());


		// Mean mag
		dummyOrder.SetType("Mag");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanMag.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		// Mean energy
		dummyOrder.SetType("Energy");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanEnergy.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("QEnergy");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanQEnergy.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("PotEnergy");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanPotEnergy.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("KinEnergy");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanKinEnergy.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummy2DOrder.Initialize(conf.L,1);
		dummy2DOrder.SetType("MeanKinEnergyDensity");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanKinEnergyDensity.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

                dummyOrder.SetType("QNewEnergy");
                dummyOrder.CreateFiles();
                dummyStringBis  = "./" + vtFolder[iii] + "/Obj_MeanQNewEnergy.odat";
                ExportObject(dummyOrder,dummyStringBis.c_str());
                
                dummyOrder.SetType("NewKinEnergy");
                dummyOrder.CreateFiles();
                dummyStringBis  = "./" + vtFolder[iii] + "/Obj_MeanNewKinEnergy.odat";
                ExportObject(dummyOrder,dummyStringBis.c_str());

		// Mean N_j
		dummyOrder.SetType("N0");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanN0.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("N1");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanN1.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("N2");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanN2.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("N3");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanN3.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("PlaquetteIndex");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanPlaquetteIndex.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("Star3Index");
		dummyOrder.CreateFiles();
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanStar3Index.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		// Mean mag stag
		if(conf.initCondType == 0)
		{
			if(input.ny==1)
			{
				dummyOrder.SetType("MeanStagMag");
				dummyOrder.CreateFiles();
				dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanStagMag.odat";
				ExportObject(dummyOrder,dummyStringBis.c_str());

				dummyOrder.SetType("MeanChainEnergy");
				dummyOrder.CreateFiles();
				dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanChainEnergy.odat";
				ExportObject(dummyOrder,dummyStringBis.c_str());

				dummyOrder.SetCorrLength(input.nx + input.p);
				dummyOrder.SetType("CorrMag");
				dummyOrder.CreateFiles();
				dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrMag.odat";
				ExportObject(dummyOrder,dummyStringBis.c_str());

				dummyOrder.SetType("CorrStagMag");
				dummyOrder.CreateFiles();
				dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrStagMag.odat";
				ExportObject(dummyOrder,dummyStringBis.c_str());
			}
			else
			{
				dummyOrder.SetType("MeanCoreMag");
				dummyOrder.CreateFiles();
				dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanCoreMag.odat";
				ExportObject(dummyOrder,dummyStringBis.c_str());
			}
		}

		// For the 2D order parameters
		dummy2DOrder.Initialize(conf.nx*conf.nx,1);
		dummy2DOrder.SetType("MeanPart");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanPart.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		dummy2DOrder.Initialize(conf.L,1);
		dummy2DOrder.SetType("MeanLocalNf");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanLocalNf.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		dummy2DOrder.Initialize(conf.L,conf.NbOfNeights);
		dummy2DOrder.SetType("MeanLocalDimer");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_MeanLocalDimer.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		dummy2DOrder.Initialize(conf.L,1);
		dummy2DOrder.SetType("Star3Lattice");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_Star3.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		dummyOrder.SetCorrLength(conf.N);
		dummyOrder.SetType("CorrSzSz");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrSzSz.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummyOrder.SetType("CorrSxSx");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrSxSx.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummy2DOrder.SetType("CorrMeanSz");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrMeanSzSz.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		dummyOrder.SetType("CorrDimer");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrDimer.odat";
		ExportObject(dummyOrder,dummyStringBis.c_str());

		dummy2DOrder.Initialize(conf.L,conf.NbOfNeights);
		dummy2DOrder.SetType("CorrMeanDimer");
		dummyStringBis	= "./" + vtFolder[iii] + "/Obj_CorrMeanDimer.odat";
		ExportObject(dummy2DOrder,dummyStringBis.c_str());

		// >>>> For the RNG engine
		dummyString = "./" + vtFolder[iii] + "/Obj_Engine.odat";
		conf.ExportRNG(dummyString.c_str());
	}

	string scriptFilename	= "./run_Part_" + ToString(pid) + ".sh";
	string commandsFilename	= "comm_Part_" + ToString(pid) + ".sh";

	if(input.compType==0)
	{
		ofstream commandsFile(commandsFilename.c_str(),ios::trunc);
		int procCounter = 0;
		for(int iii = 0; iii < input.VtN; ++iii)
		{
			commandsFile 	<< "nice -19 ./" << input.execName << " " << outFolder
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
			commandsFile 	<< "nice -19 ./" << input.execName << " " << outFolder
							<< "/Obj_conf.odat" << " " << vtFolder[iii] << endl;
		}
		commandsFile.close();
	}
	else if(input.compType==3)
	{
		ofstream commandsFile(commandsFilename.c_str(),ios::trunc);
		for(int iii = 0; iii < input.VtN; ++iii)
		{
			commandsFile 	<< "srun --exclusive -N1 -n1 nice -19 ./" << input.execName << " " << outFolder
							<< "/Obj_conf.odat" << " " << vtFolder[iii] <<  endl;
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
