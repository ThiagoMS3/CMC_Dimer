#include "main.h"

int main(int argc, char **argv)
{
	string confFile		= argv[1];
	string WorkFolder	= argv[2];

	SysConf 	conf;
	ImportObject(conf,confFile.c_str());

	OrderParam	MeanMag;
	OrderParam	MeanEnergy;

	OrderParam	MeanN0;
	OrderParam	MeanN1;
	OrderParam	MeanN2;
	OrderParam	MeanN3;

	OrderParam	MeanQEnergy;
    OrderParam	MeanQKinEnergy;
    OrderParam	MeanQPotEnergy;

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
	TwoDOrderParam	MeanSublattice;

	/*
	 *    > Name of the files
	 *    > Import it
	 *    > Update number of measurements
	 */

	string OrderParamFile = WorkFolder + "/Obj_MeanMag.odat";
	ImportObject(MeanMag,OrderParamFile.c_str());
	MeanMag.ConvertMean();
	ExportObject(MeanMag,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanEnergy.odat";
	ImportObject(MeanEnergy,OrderParamFile.c_str());
	MeanEnergy.ConvertMean();
	ExportObject(MeanEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN0.odat";
	ImportObject(MeanN0,OrderParamFile.c_str());
	MeanN0.ConvertMean();
	ExportObject(MeanN0,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN1.odat";
	ImportObject(MeanN1,OrderParamFile.c_str());
	MeanN1.ConvertMean();
	ExportObject(MeanN1,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN2.odat";
	ImportObject(MeanN2,OrderParamFile.c_str());
	MeanN2.ConvertMean();
	ExportObject(MeanN2,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanN3.odat";
	ImportObject(MeanN3,OrderParamFile.c_str());
	MeanN3.ConvertMean();
	ExportObject(MeanN3,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanQEnergy.odat";
	ImportObject(MeanQEnergy,OrderParamFile.c_str());
	MeanQEnergy.ConvertMean();
	ExportObject(MeanQEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanKinEnergy.odat";
	ImportObject(MeanQKinEnergy,OrderParamFile.c_str());
	MeanQKinEnergy.ConvertMean();
	ExportObject(MeanQKinEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanPotEnergy.odat";
	ImportObject(MeanQPotEnergy,OrderParamFile.c_str());
	MeanQPotEnergy.ConvertMean();
	ExportObject(MeanQPotEnergy,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanSublattice.odat";
	ImportObject(MeanSublattice,OrderParamFile.c_str());
	MeanSublattice.ConvertMean();
	ExportObject(MeanSublattice,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanLocalNf.odat";
	ImportObject(MeanLocalNf,OrderParamFile.c_str());
	MeanLocalNf.ConvertMean();
	ExportObject(MeanLocalNf,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_MeanLocalDimer.odat";
	ImportObject(MeanLocalDimer,OrderParamFile.c_str());
	MeanLocalDimer.ConvertMean();
	ExportObject(MeanLocalDimer,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_CorrSzSz.odat";
	ImportObject(SzSzCorrelation,OrderParamFile.c_str());
	SzSzCorrelation.ConvertMean();
	ExportObject(SzSzCorrelation,OrderParamFile.c_str());

	OrderParamFile = WorkFolder + "/Obj_CorrDimer.odat";
	ImportObject(DimerCorrelation,OrderParamFile.c_str());
	DimerCorrelation.ConvertMean();
	ExportObject(DimerCorrelation,OrderParamFile.c_str());

	vector<double> dummyvec(3,0);
	int meas = MeanSublattice.GetNbOfMeas();
	MeanSublattice.GetMean(dummyvec);

	cout << " mean sub = " << dummyvec[0]*meas << " " << dummyvec[1]*meas << " " << dummyvec[2]*meas << " " << " , " << meas << endl;
	return 0;
}
