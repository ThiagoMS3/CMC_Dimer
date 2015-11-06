#include "SysConf.h"

// *******************************************************************
// General cluster CMC methods
// ---> The descriptions are in the *.h file
// *******************************************************************

// *******************************************************************
// Constructors and associated methods
// *******************************************************************
SysConf::SysConf(const input_params& input,const string& TypeIn)
{
	SetParam(input,TypeIn);
	SetRNG(input.RNGSeed);
}

SysConf::SysConf(const SysConf& original)
{
	Copy(original);
	SetRNG(original.RNGSeed);
}

void SysConf::SetParam(const input_params& input,const string& TypeIn)
{

	NbOfNeights	= input.NbOfNeights;
	N 		= input.N;

	SimType 	= TypeIn;
	if(SimType.compare("Moessner")==0)
	{
		nx		= input.nx;
		ny		= input.ny;
		p		= -1;
		L		= input.nx*input.ny;
		NDimers		= input.nx*input.ny;
		LBorder		= 0;
		LTot		= L+1;
		initCondType	= 0; // Fixed
		ConfigFile = "";
	}

	else if(SimType.compare("Part")==0)
	{
		nx		= input.nx;
		ny		= input.ny;
		p		= input.p;
		LTot	= input.nx*input.ny + input.nx*input.p + input.ny*input.p + input.nx + input.ny + input.p +1;
		LBorder = 2*(input.nx + input.p + input.ny);
		L 	= LTot - LBorder;
		NDimers	= input.nx*input.ny + input.nx*input.p + input.ny*input.p;
		initCondType = input.initCondType;

		equivPart.resize(nx*ny*N,0);
		horizontalDimerPos.resize(ny*N,0);
		if(ny==1)
		{
			spinChain.resize((nx+p)*N,1);
		}
		ConfigFile = "";
	}

	else if(SimType.compare("Manual")==0)
	{
		nx		= input.nx;
		ny		= input.ny;
		p		= -1;
		L		= input.nx*input.ny;
		NDimers		= input.nx*input.ny;
		LBorder     = 3;
		LTot        = L+3;
		initCondType = 0; // Fixed
	}

	BorderType	= input.BorderType;

	if(BorderType==0||BorderType==2)
	{
		neighsSpin.resize(NbOfNeights,0);
	}

	neighboursTable.resize(LTot*NbOfNeights,-1);
	flippableSpinLists.resize(N*LTot,-1);
	flippableSpinPositions.resize(N*LTot,-1);

	spinConf.resize(LTot*N,0);

	CoreStart.resize(1,-1);
	CoreEnd.resize(1,-1);

	CoreSize = 0;

	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	toRemove.resize(N,vector<int>(NbOfNeights/2,-1));
	toAdd.resize(N,vector<int>(NbOfNeights/2,-1));

	Kz		= input.Dbeta;
	Kt		= -log(tanh((input.Dbeta)))/2.;
	Vt		= 0;

	alpha0  = input.alpha0;
	alpha3  = input.alpha3;

	shift0  = input.shift0;
	shift3  = input.shift3;

	measurementsDone = 0;
}

void SysConf::Copy(const SysConf& original)
{
	nx		= original.nx;
	ny		= original.ny;
	p		= original.p;
	L			= original.L;
	LTot			= original.LTot;
	LBorder			= original.LBorder;
	NDimers			= original.NDimers;

	NbOfNeights		= original.NbOfNeights;
	N			= original.N;

	TotalOldNbF		= original.TotalOldNbF;
	TotalNewNbF		= original.TotalNewNbF;

	neighboursTable 	= original.neighboursTable;
	spinConf		= original.spinConf;

	flippableSpinLists	= original.flippableSpinLists;
	flippableSpinPositions	= original.flippableSpinPositions;

	SimType			= original.SimType;

	CoreStart		= original.CoreStart;
	CoreEnd			= original.CoreEnd;

	CoreSize		= original.CoreSize;

	RNGSeed			= original.RNGSeed;

    if(SimType.compare("Part")==0)
	{
		equivPart.resize(nx*ny*N,0);
		horizontalDimerPos.resize(ny*N,0);
		if(ny==1)
		{
			spinChain.resize((nx+p)*N,1);
		}
	}
}

void SysConf::SetConfigFileManual(string& configName)
{
	ConfigFile = configName;
}

void SysConf::SetNeighbours()
{
	neighboursTable.resize(LTot*NbOfNeights,-1);
	if(SimType.compare("Part")==0)
	{
		SetNeighboursPart();
	}
	else if(SimType.compare("Moessner")==0)
	{
		SetNeighboursMoessner();
	}
	else if(SimType.compare("Manual")==0)
	{
		SetNeighboursManual();
	}
}

void SysConf::SetFlipTables()
{
	bool IsFlippable = false;
	int flippableCandidate = 0;
	int dummyNbF = 0;

	flippableSpinLists.resize(N*LTot,-1);
	flippableSpinPositions.resize(N*LTot,-1);

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int iii = 0; iii < L; ++iii)
		{
			flippableCandidate = idxConv(N,iii,nnn);
			IsFlippable = TestFlippable(iii,nnn);
			if(IsFlippable)
			{
				flippableSpinLists[dummyNbF] = flippableCandidate;
				flippableSpinPositions[flippableCandidate] = dummyNbF;
				++dummyNbF;
			}
		}
	}
}

void SysConf::SetSecondNeightbours()
{
	int firstPos = -1;
	vector<int> secondNeighs(6,0);
	vector<int> firstNeighs(6,0);
	secondNeighboursTable.resize(L*NbOfNeights*2,-1);
	int index = 0;
	bool IsZeroNeigh = false;

	for(int lll = 0; lll < L; ++lll)
	{
		index = 0;
		for(int kkk = 0; kkk < NbOfNeights; ++kkk)
		{
			firstNeighs[kkk] = neighboursTable[idxConv(NbOfNeights,lll,kkk)];
		}

		for(int kkk = 0; kkk < NbOfNeights; ++kkk)
		{
			firstPos = firstNeighs[kkk];

			for(int iii = 0; iii < NbOfNeights; ++iii)
			{
				secondNeighs[iii] = neighboursTable[idxConv(NbOfNeights,firstPos,iii)];
			}

			for(int iii = 0; iii < NbOfNeights; ++iii)
			{
				IsZeroNeigh = false;

				if(secondNeighs[iii]==lll)
				{
					// then we shouldn't add this term
					IsZeroNeigh = true;
				}

				if(!IsZeroNeigh)
				{
					for(int jjj = 0; jjj < NbOfNeights; ++jjj)
					{
						if(secondNeighs[iii]==firstNeighs[jjj])
						{
							// then we shouldn't add this term
							IsZeroNeigh = true;
							break;
						}
					}
				}

				if(!IsZeroNeigh)
				{
					for(int jjj = 0; jjj < index; ++jjj)
					{
						if(secondNeighs[iii]==secondNeighboursTable[idxConv(NbOfNeights*2,lll,jjj)])
						{
							// then we shouldn't add this term
							IsZeroNeigh = true;
							break;
						}
					}
				}

				if(!IsZeroNeigh)
				{
					secondNeighboursTable[idxConv(NbOfNeights*2,lll,index)]=secondNeighs[iii];
					++index;
				}
			}

		}
	}
}

//void SysConf::Move(int& pos, int direction, int distance)
//{
//	for(int iii = 0; iii < distance; ++iii)
//	{
//		pos = neighboursTable[idxConv(NbOfNeights,pos,direction)];
//	}
//}
//
//void SysConf::SetNthNeightbours()
//{
//	nbOfNeighLevels = min(nx,ny)/2;
//	int idxVecSize = nbOfNeighLevels*(nbOfNeighLevels-1)*3;
//
//	int kkk = -1;
//	int initialPos = -1;
//
//	for(int iii = 0; iii < L; ++iii)
//	{
//		nthNeighboursTable[iii].resize(idxVecSize,-1);
//
//		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
//		{
//			nthNeighboursTable[iii][jjj] = neighboursTable[idxConv(NbOfNeights,iii,jjj)];
//		}
//
//
//		for(int nnn = 2; nnn < nbOfNeighLevels; ++nnn)
//		{
//			for(int jjj = 0; jjj < nnn; ++jjj)
//			{
//
//			}
//		}
//	}
//
//	int firstPos = -1;
//	vector<int> secondNeighs(6,0);
//	vector<int> firstNeighs(6,0);
//	secondNeighboursTable.resize(L*NbOfNeights*2,-1);
//	int index = 0;
//	bool IsZeroNeigh = false;
//
//	for(int lll = 0; lll < L; ++lll)
//	{
//		index = 0;
//		for(int kkk = 0; kkk < NbOfNeights; ++kkk)
//		{
//			firstNeighs[kkk] = neighboursTable[idxConv(NbOfNeights,lll,kkk)];
//		}
//
//		for(int kkk = 0; kkk < NbOfNeights; ++kkk)
//		{
//			firstPos = firstNeighs[kkk];
//
//			for(int iii = 0; iii < NbOfNeights; ++iii)
//			{
//				secondNeighs[iii] = neighboursTable[idxConv(NbOfNeights,firstPos,iii)];
//			}
//
//			for(int iii = 0; iii < NbOfNeights; ++iii)
//			{
//				IsZeroNeigh = false;
//
//				if(secondNeighs[iii]==lll)
//				{
//					// then we shouldn't add this term
//					IsZeroNeigh = true;
//				}
//
//				if(!IsZeroNeigh)
//				{
//					for(int jjj = 0; jjj < NbOfNeights; ++jjj)
//					{
//						if(secondNeighs[iii]==firstNeighs[jjj])
//						{
//							// then we shouldn't add this term
//							IsZeroNeigh = true;
//							break;
//						}
//					}
//				}
//
//				if(!IsZeroNeigh)
//				{
//					for(int jjj = 0; jjj < index; ++jjj)
//					{
//						if(secondNeighs[iii]==secondNeighboursTable[idxConv(NbOfNeights*2,lll,jjj)])
//						{
//							// then we shouldn't add this term
//							IsZeroNeigh = true;
//							break;
//						}
//					}
//				}
//
//				if(!IsZeroNeigh)
//				{
//					secondNeighboursTable[idxConv(NbOfNeights*2,lll,index)]=secondNeighs[iii];
//					++index;
//				}
//			}
//
//		}
//	}
//}

// *******************************************************************
// Initial state / SimType == "Moessner"
// *******************************************************************
void SysConf::SetNeighboursMoessner()
{
	// ---> Set up the neighbour relations of the spins.
	//	   For the periodic boundary conditions, the setup of the spins is a lot easier ...
	SetHexOffset();

	// Parameters
	int index = 0;

	// Set bulk
	for(int iii = 1; iii < nx-1; ++iii)
	{
		for(int jjj = 1; jjj < ny-1; ++jjj)
		{
			index = nx*jjj + iii;
			neighboursTable[idxConv(6,index,0)] = index - 1;
			neighboursTable[idxConv(6,index,1)] = index + nx - 1;
			neighboursTable[idxConv(6,index,2)] = index + nx;
			neighboursTable[idxConv(6,index,3)] = index + 1;
			neighboursTable[idxConv(6,index,4)] = index - nx + 1;
			neighboursTable[idxConv(6,index,5)] = index - nx;
		}
	}

	// Set left border --> change 4 & 5
	for(int iii = 1; iii < nx-1; ++iii)
	{
		index = iii;
		neighboursTable[idxConv(6,index,0)] = index - 1;
		neighboursTable[idxConv(6,index,1)] = index + nx - 1;
		neighboursTable[idxConv(6,index,2)] = index + nx;
		neighboursTable[idxConv(6,index,3)] = index + 1;
		if(BorderType==1)
		{
			neighboursTable[idxConv(6,index,4)] = index + (ny-1)*nx + 1;	// !!!
			neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;		// !!!
		}
		else
		{
			neighboursTable[idxConv(6,index,4)] = L;	// !!!
			neighboursTable[idxConv(6,index,5)] = L;		// !!!
		}
	}

	// Set right border --> change 1 & 2
	for(int iii = 1; iii < nx-1; ++iii)
	{
		index = ((ny-1)*nx + iii);
		neighboursTable[idxConv(6,index,0)] = index - 1;
		if(BorderType==1)
		{
			neighboursTable[idxConv(6,index,1)] = iii -1;			// !!!
			neighboursTable[idxConv(6,index,2)] = iii;			// !!!
		}
		else
		{
			neighboursTable[idxConv(6,index,1)] = L;	// !!!
			neighboursTable[idxConv(6,index,2)] = L;		// !!!
		}
		neighboursTable[idxConv(6,index,3)] = index + 1;
		neighboursTable[idxConv(6,index,4)] = index - nx + 1;
		neighboursTable[idxConv(6,index,5)] = index - nx;
	}

	// Set upper border --> change 0 & 1
	for(int jjj = 1; jjj < ny-1; ++jjj)
	{
		index = jjj*nx;
		if(BorderType==1)
		{
			neighboursTable[idxConv(6,index,0)] = index + nx - 1;		// !!!
			neighboursTable[idxConv(6,index,1)] = index + 2*nx - 1;		// !!!
		}
		else
		{
			neighboursTable[idxConv(6,index,0)] = L;	// !!!
			neighboursTable[idxConv(6,index,1)] = L;		// !!!
		}
		neighboursTable[idxConv(6,index,2)] = index + nx;
		neighboursTable[idxConv(6,index,3)] = index + 1;
		neighboursTable[idxConv(6,index,4)] = index - nx + 1;
		neighboursTable[idxConv(6,index,5)] = index - nx;
	}

	// Set lower border --> change 3 & 4
	for(int jjj = 1; jjj < ny-1; ++jjj)
	{
		index = ((jjj + 1)*nx - 1);
		neighboursTable[idxConv(6,index,0)] = index - 1;
		neighboursTable[idxConv(6,index,1)] = index + nx - 1;
		neighboursTable[idxConv(6,index,2)] = index + nx;
		if(BorderType==1)
		{
			neighboursTable[idxConv(6,index,3)] = index - nx + 1;		// !!!
			neighboursTable[idxConv(6,index,4)] = index - 2*nx + 1;		// !!!
		}
		else
		{
			neighboursTable[idxConv(6,index,3)] = L;	// !!!
			neighboursTable[idxConv(6,index,4)] = L;		// !!!
		}
		neighboursTable[idxConv(6,index,5)] = index - nx;
	}

	// > Finally, set corners
	// Upper left --> change 0, 1, 4 & 5
	index = 0;
	if(BorderType==1)
	{
		neighboursTable[idxConv(6,index,0)] = index + nx - 1;			// !!!
		neighboursTable[idxConv(6,index,1)] = index + 2*nx - 1;			// !!!
	}
	else
	{
		neighboursTable[idxConv(6,index,0)] = L;	// !!!
		neighboursTable[idxConv(6,index,1)] = L;		// !!!
	}
	neighboursTable[idxConv(6,index,2)] = index + nx;
	neighboursTable[idxConv(6,index,3)] = index + 1;
	if(BorderType==1)
	{
		neighboursTable[idxConv(6,index,4)] = index + (ny-1)*nx + 1;		// !!!
		neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;			// !!!
	}
	else
	{
		neighboursTable[idxConv(6,index,4)] = L;	// !!!
		neighboursTable[idxConv(6,index,5)] = L;		// !!!
	}


	// Lower left --> change 3, 4 & 5
	index = (nx - 1);
	neighboursTable[idxConv(6,index,0)] = index - 1;
	neighboursTable[idxConv(6,index,1)] = index + nx - 1;
	neighboursTable[idxConv(6,index,2)] = index + nx;
	if(BorderType==1)
	{
		neighboursTable[idxConv(6,index,3)] = index - nx + 1;			// !!!
		neighboursTable[idxConv(6,index,4)] = index + (ny-2)*nx + 1;		// !!! !!!
		neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;			// !!!
	}
	else
	{
		neighboursTable[idxConv(6,index,3)] = L;		// !!!
		neighboursTable[idxConv(6,index,4)] = L;	// !!!
		neighboursTable[idxConv(6,index,5)] = L;		// !!!
	}


	// Upper right --> change 0, 1 & 2
	index = (ny-1)*nx;
	if(BorderType==1)
	{
		neighboursTable[idxConv(6,index,0)] = index + nx - 1;			// !!!
		neighboursTable[idxConv(6,index,1)] = index - (ny-2)*nx - 1;		// !!! !!!
		neighboursTable[idxConv(6,index,2)] = index - (ny-1)*nx;			// !!!
	}
	else
	{
		neighboursTable[idxConv(6,index,0)] = L;		// !!!
		neighboursTable[idxConv(6,index,1)] = L;	// !!!
		neighboursTable[idxConv(6,index,2)] = L;		// !!!
	}
	neighboursTable[idxConv(6,index,3)] = index + 1;
	neighboursTable[idxConv(6,index,4)] = index - nx + 1;
	neighboursTable[idxConv(6,index,5)] = index - nx;

	// Lower right --> change 1, 2, 3 & 4
	index = (nx*ny-1);
	neighboursTable[idxConv(6,index,0)] = index - 1;
	if(BorderType==1)
	{
		neighboursTable[idxConv(6,index,1)] = nx - 2;
		neighboursTable[idxConv(6,index,2)] = nx - 1;
		neighboursTable[idxConv(6,index,3)] = (ny-1)*nx;
		neighboursTable[idxConv(6,index,4)] = (ny-2)*nx;
	}
	else
	{
		neighboursTable[idxConv(6,index,1)] = L;	// !!!
		neighboursTable[idxConv(6,index,2)] = L;		// !!!
		neighboursTable[idxConv(6,index,3)] = L;	// !!!
		neighboursTable[idxConv(6,index,4)] = L;		// !!!
	}
	neighboursTable[idxConv(6,index,5)] = index - nx;
}

void SysConf::SetInitialMoessner()
{

	Jindex.resize(3,1);
	complex<double> dummyCmplx = 1i*2.*M_PI/3.;
	Jindex[1] = exp(dummyCmplx);
	dummyCmplx = 1i*4.*M_PI/3.;
	Jindex[2] = exp(dummyCmplx);

	// TODO This version is kind of a cheat ... it creates the limit case for V/t -> - infty
	int dummyNbF = 0;
	int kkk;
	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int jjj = 0; jjj < ny/3; ++jjj)
		{
			for(int iii = 0; iii < nx/3; ++iii)
			{
				kkk = (3*jjj*nx + 3*iii)*N;
				spinConf[kkk + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + nnn;
// 				flippableSpinPositions[kkk + nnn] = dummyNbF;
// 				++dummyNbF;

				spinConf[kkk + N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + N + nnn;
// 				flippableSpinPositions[kkk + N + nnn] = dummyNbF;
// 				++dummyNbF;

				spinConf[kkk + 2*N  + nnn] = -1;

				kkk = ((3*jjj+1)*nx + 3*iii)*N;
				spinConf[kkk + nnn] = -1;

				spinConf[kkk + N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + N + nnn;
// 				flippableSpinPositions[kkk + N + nnn] = dummyNbF;
// 				++dummyNbF;

				spinConf[kkk + 2*N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + 2*N + nnn;
// 				flippableSpinPositions[kkk + 2*N + nnn] = dummyNbF;
// 				++dummyNbF;

				kkk = ((3*jjj+2)*nx + 3*iii)*N;
				spinConf[kkk + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + nnn;
// 				flippableSpinPositions[kkk + nnn] = dummyNbF;
// 				++dummyNbF;

				spinConf[kkk + N + nnn] = -1;

				spinConf[kkk + 2*N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + 2*N + nnn;
// 				flippableSpinPositions[kkk + 2*N + nnn] = dummyNbF;
// 				++dummyNbF;
			}
		}
		if(BorderType==0||BorderType==2)
		{
			kkk = L;
			spinConf[kkk + nnn] = 0;
		}
	}

	bool IsFlippable = false;
	int flippableCandidate = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		flippableCandidate = idxConv(N,0,iii);
		IsFlippable = TestFlippable(flippableCandidate,0);
		if(IsFlippable)
		{
			flippableSpinLists[dummyNbF] = flippableCandidate*N;
			flippableSpinPositions[flippableCandidate*N] = dummyNbF;
			++dummyNbF;
		}
	}

	// > Finally, copy informations over all the layers
	int pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[pos] = spin;
			flippableSpinPositions[spin] = pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;

// 	TotalOldNbF = dummyNbF;
// 	TotalNewNbF = dummyNbF;
}

void SysConf::SetNeighboursManual()
{
	ifstream inputFile(ConfigFile.c_str());
	int dummyPos;
	int kkk = 0;
	for(int iii = 0; iii < ny; ++iii)
	{
		for(int jjj = 0; jjj < nx; ++jjj)
		{
			for(int lll = 0; lll < NbOfNeights; ++lll)
			{
				inputFile >> dummyPos;
				neighboursTable[idxConv(6,kkk,lll)] = dummyPos;
			}
			Jump(inputFile,1);

			++kkk;
		}
	}
}

void SysConf::SetInitialManual(string& configName)
{
	SetConfigFileManual(configName);

	int dummyNbF = 0;
	int dumbSpin = 0;
	int dummyPos = 0;
	int kkk = 0;
	bool IsFlippable = false;
	int flippableCandidate = 0;

	ifstream inputFile(ConfigFile.c_str());

//	for(int iii = 0; iii < ny; ++iii)
//	{
//		for(int jjj = 0; jjj < nx; ++jjj)
//		{
//			for(int lll = 0; lll < NbOfNeights; ++lll)
//			{
//				inputFile >> dummyPos;
//				neighboursTable[idxConv(6,kkk,lll)] = dummyPos;
//			}
//			inputFile >> dumbSpin;
//			spinConf[kkk*N] = dumbSpin;
//
//			++kkk;
//		}
//		Jump(inputFile,1);
//	}

	for(int iii = 0; iii < ny; ++iii)
	{
		for(int jjj = 0; jjj < nx; ++jjj)
		{
			for(int lll = 0; lll < NbOfNeights; ++lll)
			{
				inputFile >> dummyPos;
				neighboursTable[idxConv(6,kkk,lll)] = dummyPos;
			}
			inputFile >> dumbSpin;
			spinConf[kkk*N] = dumbSpin;

			++kkk;
		}
//		Jump(inputFile,1);
	}


	// Extra spins
	spinConf[nx*ny*N] = 0;
	spinConf[(nx*ny+1)*N] = -1;
	spinConf[(nx*ny+2)*N] = 1;

	inputFile.close();

#if BORDER_TYPE == 4
// ANTI-PERIODIC BORDERS - SET WEIGHTS
	ReadAntiWeights();
#endif
// 	for(int iii = 0; iii < L; ++iii)
// 	{
// 		for(int lll = 0; lll < NbOfNeights; ++lll)
// 		{
// 			cout << neighboursTable[idxConv(6,iii,lll)] << " ";
// 		}
// 		cout << "| " << (int)spinConf[iii*N] << endl;
// 	}

	for(int iii = 0; iii < L; ++iii)
	{
		flippableCandidate = idxConv(N,0,iii);
		IsFlippable = TestFlippable(flippableCandidate,0);
		if(IsFlippable)
		{
			flippableSpinLists[dummyNbF] = flippableCandidate*N;
			flippableSpinPositions[flippableCandidate*N] = dummyNbF;
			++dummyNbF;
		}
	}

	// > Finally, copy informations over all the layers
	int pos = 0;
	int spin = 0;
	int index = 0;
	for(int iii = 0; iii < LTot; ++iii)
	{
		index = N*iii;
		spin = spinConf[index];
		for(int jjj = 1; jjj < N; ++jjj)
		{
			spinConf[index+jjj] = spin;
		}
	}

	for(int iii = 1; iii < N; ++iii)
	{
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[pos] = spin;
			flippableSpinPositions[spin] = pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;
}

void SysConf::ReadAntiWeights()
{
	weightTable.resize(nx*ny*NbOfNeights,1);
	string WeightFile = ConfigFile + "_weights";
	ifstream inputFile(WeightFile.c_str());

	int dummyWeight = 0;

	for(int kkk = 0; kkk < L; ++ kkk)
	{
		for(int lll = 0; lll < NbOfNeights; ++lll)
		{
			inputFile >> dummyWeight;
			weightTable[idxConv(6,kkk,lll)] = dummyWeight;
		}
	}
	inputFile.close();
}

void SysConf::SetAntiWeights()
{
	weightTable.resize(nx*ny*NbOfNeights,1);
	int dummyPos = 0;
	int idx = 0;
	for(int iii = 0; iii < nx/2; ++iii)
	{
		idx = 2*iii;
		dummyPos = neighboursTable[idxConv(NbOfNeights,idx,4)];
		weightTable[idxConv(NbOfNeights,idx,4)] = -1;
		weightTable[idxConv(NbOfNeights,dummyPos,1)] = -1;

		dummyPos = neighboursTable[idxConv(NbOfNeights,idx,5)];
		weightTable[idxConv(NbOfNeights,idx,5)] = -1;
		weightTable[idxConv(NbOfNeights,dummyPos,2)] = -1;

		dummyPos = neighboursTable[idxConv(NbOfNeights,idx,0)];
		weightTable[idxConv(NbOfNeights,idx,0)] = -1;
		weightTable[idxConv(NbOfNeights,dummyPos,3)] = -1;

		idx = 2*iii + 1;
		dummyPos = neighboursTable[idxConv(NbOfNeights,idx,5)];
		weightTable[idxConv(NbOfNeights,idx,5)] = -1;
		weightTable[idxConv(NbOfNeights,dummyPos,2)] = -1;
	}

	for(int iii = 0; iii < nx + 1; ++ iii)
	{
		cout << iii << "\t: ";
		for(int jjj = 0; jjj < 6; ++jjj)
		{
			cout << weightTable[idxConv(NbOfNeights,iii,jjj)] << " ";
		}
		cout << endl;
	}
	cout << "..." << endl;
	for(int iii = nx*(ny-1) - 2; iii < nx*ny; ++ iii)
	{
		cout << iii << "\t: ";
		for(int jjj = 0; jjj < 6; ++jjj)
		{
			cout << weightTable[idxConv(NbOfNeights,iii,jjj)];
		}
		cout << endl;
	}
}

// void SysConf::SetInitialManual(const input_params& input, string& ConfigFile)
// {
// 	int dummyNbF = 0;
// 	int dumbSpin = 0;
// 	int kkk = 0;
// 	bool IsFlippable = false;
// 	int flippableCandidate = 0;
//
// 	ifstream inputFile(ConfigFile.c_str());
//
// 	for(int iii = 0; iii < ny; ++iii)
// 	{
// 		for(int jjj = 0; jjj < nx; ++jjj)
// 		{
//
// 			inputFile >> dumbSpin;
// 			spinConf[kkk*N] = dumbSpin;
// 			++kkk;
// 		}
// 		Jump(inputFile,1);
// 	}
//
// 	inputFile.close();
//
// 	for(int iii = 0; iii < L; ++iii)
// 	{
// 		flippableCandidate = idxConv(N,0,iii);
// 		IsFlippable = TestFlippable(flippableCandidate,0);
// 		if(IsFlippable)
// 		{
// 			flippableSpinLists[dummyNbF] = flippableCandidate*N;
// 			flippableSpinPositions[flippableCandidate*N] = dummyNbF;
// 			++dummyNbF;
// 		}
// 	}
//
// 	// > Finally, copy informations over all the layers
// 	int pos = 0;
// 	int spin = 0;
// 	for(int iii = 1; iii < N; ++iii)
// 	{
// 		for(int jjj = 0; jjj < dummyNbF; ++jjj)
// 		{
// 			pos = iii*dummyNbF + jjj;
// 			spin = flippableSpinLists[jjj] + iii;
// 			flippableSpinLists[pos] = spin;
// 			flippableSpinPositions[spin] = pos;
// 		}
// 	}
//
// 	TotalOldNbF = N*dummyNbF;
// 	TotalNewNbF = N*dummyNbF;
// }

// *******************************************************************
// Initial state / SimType == "Partition"
// *******************************************************************
void SysConf::SetHexOffset()
{
	offset.resize(6);
	offset[0].x = -1; offset[0].y = -1;
	offset[1].x = 0; offset[1].y = -1;
	offset[2].x = 1; offset[2].y = 0;
	offset[3].x = 1; offset[3].y = 1;
	offset[4].x = 0; offset[4].y = 1;
	offset[5].x = -1; offset[5].y = 0;
};

void SysConf::SetNeighboursPart()
{
	SetHexOffset();
	// ---> Set up the neighbour relations of the spins.
	//	   First, we have to set up a translation table between the (x,y) coordinates and
	//	the list positions. This table will be implemented as a binary search tree ('map')

	// > Parameters
	int positionTable = 0;
	coord position = {0,0};
	rowStart.resize(p + nx + 1,-1);
	rowEnd.resize(p + nx + 1,-1);
	rowSize.resize(p + nx + 1,-1);

	// > Translation table
	map<coord,int> translate;

	// > Set limits
	for(int iii = 0; iii < p; ++iii)
	{
		rowStart[iii] = 0;
	}
	for(int iii = p; iii < p + nx + 1; ++iii)
	{
		rowStart[iii] = iii - p;
	}

	for(int iii = 0; iii < nx + 1; ++iii)
	{
		rowEnd[iii] = ny + iii;
	}
	for(int iii = nx + 1; iii < p + nx + 1; ++iii)
	{
		rowEnd[iii] = nx + ny;
	}

	for(int iii = 0; iii < p + nx + 1; ++iii)
	{
		rowSize[iii] = rowEnd[iii] - rowStart[iii];
	}
	// > Set up flippable spins
	for(int iii = 1; iii < p + nx; ++iii)
	{
		position.x = iii;
		for(int jjj = rowStart[iii] + 1; jjj < rowEnd[iii]; ++jjj)
		{
			position.y = jjj;
			translate.insert(pair<coord,int>(position,positionTable));
			++positionTable;
		}
	}

	// > Set up borders
	// First line
	position.x = 0;
	for(int jjj = rowStart[0]; jjj < rowEnd[0]; ++jjj)
	{
		position.y = jjj;
		translate.insert(pair<coord,int>(position,positionTable));
		++positionTable;
	}

	// Right border
	for(int iii = 0; iii < p + nx; ++iii)
	{
		position.x = iii;
		position.y = rowEnd[iii];
		translate.insert(pair<coord,int>(position,positionTable));
		++positionTable;
	}

	// Last line
	position.x = p + nx;
	for(int jjj = rowEnd[p + nx]; jjj > rowStart[p + nx]; --jjj)
	{
		position.y = jjj;
		translate.insert(pair<coord,int>(position,positionTable));
		++positionTable;
	}

	// Left border
	for(int iii = p + nx; iii > 0; --iii)
	{
		position.x = iii;
		position.y = rowStart[iii];
		translate.insert(pair<coord,int>(position,positionTable));
		++positionTable;
	}

	// > Now set neighbouring relations
	position.x = 0;
	position.y = 0;

	int index  = 0;
	for(int iii = 1; iii < p + nx; ++iii)
	{
		position.x = iii;
		for(int jjj = rowStart[iii] + 1; jjj < rowEnd[iii]; ++jjj)
		{
			position.y = jjj;
			index = translate[position];
			for(int kkk = 0; kkk < NbOfNeights; ++kkk)
			{
				neighboursTable[idxConv(NbOfNeights,index,kkk)] = translate[position + offset[kkk]];
			}
		}
	}

// 		cout << "ny == 1!" << endl;
	// > Finally, set the upper outer borders (needed to calculate the equivalent partitions)
	position.x = 0;
	for(int jjj = 1; jjj < ny; ++jjj)
	{
		position.y = jjj;
		index = translate[position];
		neighboursTable[idxConv(NbOfNeights,index,2)] = translate[position + offset[2]];
	}

	for(int iii = 0; iii < nx; ++iii)
	{
		position.x = iii;
		position.y = rowEnd[iii];
		index = translate[position];
		neighboursTable[idxConv(NbOfNeights,index,2)] = translate[position + offset[2]];
	}

// 	for(int iii = 0; iii < LTot; ++iii)
// 	{
// 		cout << iii << " : ";
// 		for(int jjj = 0; jjj < 6; ++jjj)
// 		{
// 			cout << neighboursTable[idxConv(6,iii,jjj)] << " ";
// 		}
// 		cout << endl;
// 	}
}
/* void SysConf::SetInitialPartStar3(input_params input)
{
	int dummyNbF = 0;
	int firstSpin = 0;
	int spinPos = firstSpin;

	for(int iii = 0; iii < nx/3; ++iii)
	{
		First column
		for(int jjj = 0; jjj < p/3; ++jjj)
		{
			First column
			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = -1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
		}

		firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,3)];
		spinPos = firstSpin;

		Second column
		for(int jjj = 0; jjj < p/3; ++jjj)
		{
			First column
			spinConf[idxConv(N,spinPos,0)] = -1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
		}

		firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,3)];
		spinPos = firstSpin;

		Third column
		for(int jjj = 0; jjj < p/3; ++jjj)
		{
			First column
			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = -1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			spinConf[idxConv(N,spinPos,0)] = 1;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
		}

		firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,3)];
		spinPos = firstSpin;
	}

	bool IsFlippable = false;
	int flippableCandidate = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		flippableCandidate = idxConv(N,0,iii);
		IsFlippable = TestFlippable(flippableCandidate,0);
		if(IsFlippable)
		{
			flippableSpinLists[dummyNbF] = flippableCandidate*N;
			flippableSpinPositions[flippableCandidate*N] = dummyNbF;
			++dummyNbF;
		}
	}

	> Finally, copy informations over all the layers
	int pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[pos] = spin;
			flippableSpinPositions[spin] = pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;
} */

void SysConf::SetInitialPart()
{
	// Set staggered order parameter
	random::bernoulli_distribution<> bernoulli(0.5);

	int dummyNbF = 0;

	int nbOfDoublesUp = 0;
	int nbOfDoublesDown = 0;

	int lastFlipped = 0;
	// ---> Set borders
	for(int iii = L; iii < LTot; ++iii)
	{
		spinConf[idxConv(N,iii,0)] = pow(-1,iii);
	}

	// --> Now set the other spins
	//	> Each triangle of the lattice, formed by spins 'spin0', 'spin1' and 'spin2',
	//	  must have a total magnetization equal to +-1/2. This way, we'll have only one
	//	  frustration per triangle and we'll obey the dimer relations in the hexagonal
	//	  lattice. We have to take this constraint into account when we choose a random
	//	  configuration.
	//	> Since each 'spin0' is inside 6 triangles, we have to do this test 6 times in the
	//	  worst case

	idxList	flippableCandidate;
	flippableCandidate.reserve(L);

	bool ThereIsAChoice = true;
	int	spin1 = 0, spin2 = 0;
	int 	spin0 = 0;

	while(spin0 < L)
	{
		ThereIsAChoice = true;
		nbOfDoublesUp = 0;
		nbOfDoublesDown = 0;

		// Test for the first triangle
		spin1 = neighboursTable[idxConv(NbOfNeights,spin0,0)];
		spin2 = neighboursTable[idxConv(NbOfNeights,spin0,NbOfNeights - 1)];
		if((spinConf[idxConv(N,spin1,0)]==spinConf[idxConv(N,spin2,0)])&&spinConf[idxConv(N,spin2,0)]!=0)
		{
			if(spinConf[idxConv(N,spin1,0)]==1)
			{
				++nbOfDoublesUp;
			}
			else
			{
				++nbOfDoublesDown;
			}

			// Then we're not free to choose the spin
			spinConf[idxConv(N,spin0,0)] = -spinConf[idxConv(N,spin1,0)];
			ThereIsAChoice = false;
		}

		for(int jjj = 1; jjj < NbOfNeights; ++jjj)
		{
			spin1 = neighboursTable[idxConv(NbOfNeights,spin0,jjj)];
			spin2 = neighboursTable[idxConv(NbOfNeights,spin0,jjj-1)];
			if((spinConf[idxConv(N,spin1,0)]==spinConf[idxConv(N,spin2,0)])&&spinConf[idxConv(N,spin2,0)]!=0)
			{
				if(spinConf[idxConv(N,spin1,0)]==1)
				{
					++nbOfDoublesUp;
				}
				else
				{
					++nbOfDoublesDown;
				}

				// Then we're not free to choose the spin
				spinConf[idxConv(N,spin0,0)] = -spinConf[idxConv(N,spin1,0)];
				ThereIsAChoice = false;
			}
		}

		if(ThereIsAChoice)
		{
			// Then we can choose an spin
			// ---> Also, this spin might be flippable, so append it to the
			//	   candidates list
			spinConf[idxConv(N,spin0,0)] = -1 + 2*bernoulli(m_rng);
			flippableCandidate.push_back(spin0);
			lastFlipped = spin0;
		}
		if(nbOfDoublesDown>0&&nbOfDoublesUp>0)
		{
			for(int iii = lastFlipped + 1; iii < spin0 + 1; ++iii)
			{
				spinConf[idxConv(N,iii,0)] = 0;
			}
			spinConf[idxConv(N,lastFlipped,0)] = -spinConf[idxConv(N,lastFlipped,0)];
			spin0 = lastFlipped + 1;
		}
		else
		{
			++spin0;
		}
	}

	// Find which spins are flippable
	if(BorderType==0||BorderType==2)
	{
		for(int iii = L; iii < LTot; ++iii)
		{
			spinConf[idxConv(N,iii,0)] = 0;
		}
	}

	bool IsFlippable = false;

	for(uint iii = 0; iii < flippableCandidate.size(); ++iii)
	{
		IsFlippable = TestFlippable(flippableCandidate[iii],0);
		if(IsFlippable)
		{
			flippableSpinLists[dummyNbF] = flippableCandidate[iii]*N;
			flippableSpinPositions[flippableCandidate[iii]*N] = dummyNbF;
			++dummyNbF;
		}
	}

	// > Finally, copy informations over all the layers
	int pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}

		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[pos] = spin;
			flippableSpinPositions[spin] = pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;
}

void SysConf::SetInitialPartEmpty()
{
	// ---> Set borders
	for(int iii = L; iii < LTot; ++iii)
	{
		spinConf[idxConv(N,iii,0)] = pow(-1,iii);
	}

    int lineSize = ny;
    int spinPos = 0;
    int spinNeigh = 0;

    // First part : Number of spins on each line raises
    for(int iii = 0; iii < p - 1; ++iii)
    {
        // Alternating spins
        for(int jjj = 0; jjj < ny;++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }

        // Constant line
        for(int jjj = ny; jjj < lineSize; ++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }
        ++lineSize;
    }

    // Second part : Number of spins on each line stays constant
    for(int iii = p - 1; iii < nx; ++iii)
    {
        // Alternating spins
        for(int jjj = 0; jjj < ny;++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }

        // Constant line
        for(int jjj = ny; jjj < lineSize; ++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }
    }

    // Third part : Number of spins on each line diminshes
    for(int iii = nx; iii < nx + p - 1; ++iii)
    {
        // Alternating spins
        --lineSize;
        for(int jjj = 0; jjj < ny;++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }

        // Constant line
        for(int jjj = ny; jjj < lineSize; ++jjj)
        {
            spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
            spinConf[idxConv(N,spinPos,0)] = spinConf[idxConv(N,spinNeigh,0)];
            ++spinPos;
        }
    }

	if(BorderType==0||BorderType==2)
	{
		for(int iii = L; iii < LTot; ++iii)
		{
			spinConf[idxConv(N,iii,0)] = 0;
		}
	}

	// ---> Test the flippable spins
	//      Actually, there's only 1 flippable spin, located at nx steps from "L", in the
	//    direction 3

	int dummyFlippable = 0;
	for(int iii = 0; iii < nx - 1; ++iii)
	{
	    dummyFlippable = neighboursTable[idxConv(NbOfNeights,dummyFlippable,3)];
	}

	flippableSpinLists[0] = dummyFlippable*N;
    flippableSpinPositions[dummyFlippable*N] = 0;

	// > Finally, copy informations over all the layers
	int spin = 0;

	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}

        spin = flippableSpinLists[0] + iii;
        flippableSpinLists[iii] = spin;
        flippableSpinPositions[spin] = iii;
	}

	TotalOldNbF = N;
	TotalNewNbF = N;
}
/* void SysConf::SetInitialPartOldAntiferro(input_params input)
{
	int dummyNbF = 0;

	CoreStart.resize(nx+p-1,-1);
	CoreEnd.resize(nx+p-1,-1);

	int CoreStartCount = 0;

	// ---> Set borders
	for(int iii = L; iii < LTot; ++iii)
	{
		spinConf[idxConv(N,iii,0)] = pow(-1,iii);
	}

	if(ny>1)
	{
		int sectorSize = nx/2*(nx/2+1)/2;
		int spinPos = 0;
		int spinNeigh = 0;

		vector<int> sectors(6*sectorSize,-1);

		// ---> First sector
		int firstPos = ny/2-1;
		int lineLength = 1;
		int counter = 0;
		for(int iii = 0; iii < ny/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,0,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
			}

			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
			++CoreStartCount;

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,1)];
		}

		// ---> Second sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < p/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,1,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,3)];
			}

			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
			++CoreStartCount;
			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			++CoreStartCount;


			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,2)];
		}

		--CoreStartCount;

		// ---> Third sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < nx/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,2,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
			}

			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			++CoreStartCount;

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,3)];
		}

		// ---> Fourth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii <ny/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,3,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,4)];
		}

		// ---> Fifth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < p/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,4,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,3)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,3)];
			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,4)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,5)];
		}

		++CoreStartCount;

		// ---> Sixth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < nx/2;++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				sectors[idxConv(sectorSize,5,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,4)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,0)];
		}

		int defaultSpin = (int)spinConf[idxConv(N,sectors[0],0)];
		int spinCompare = 0;

		// ---> Set other spins
		for(int pos = 0; pos < L; ++pos)
		{
			if(spinConf[idxConv(N,pos,0)]==0)
			{
				// ---> Then we have an empty cell
				spinConf[idxConv(N,pos,0)] = defaultSpin;
				for(int iii = 0; iii  < NbOfNeights; ++iii)
				{
					spinCompare = neighboursTable[idxConv(NbOfNeights,pos,iii)];

					// ---> Search for neighbour spins that are different
					if((int)spinConf[idxConv(N,spinCompare,0)]==-defaultSpin)
					{
						spinConf[idxConv(N,pos,0)] = defaultSpin;
						break;
					}
	//
	// 				// ---> If none where found, we have to flip this spin
					spinConf[idxConv(N,pos,0)] = -defaultSpin;
				}
			}
		}
	}
	else
	{
		// There's a easier way ...
		int posAbove = 0;
		int pos = 0;
		int baseHeight = 0;

		for(int iii = 0; iii < nx; ++iii)
		{
			baseHeight = max(p - iii,0);
			posAbove = L + iii + 1;

			// Spins above the horizontal dimer
			for(int jjj = 0; jjj < p - baseHeight; ++jjj)
			{
				pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
				spinConf[idxConv(N,pos,0)] = -spinConf[idxConv(N,posAbove,0)];
				posAbove = pos;
			}

			// Spin at the horizontal dimer
			pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
			spinConf[idxConv(N,pos,0)] = spinConf[idxConv(N,posAbove,0)];
			posAbove = pos;

			// Spins below the horizontal dimer
			for(int jjj = p - baseHeight + 1; jjj < p; ++jjj)
			{
				pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
				spinConf[idxConv(N,pos,0)] = -spinConf[idxConv(N,posAbove,0)];
				posAbove = pos;
			}
		}
	}

	if(BorderType==0||BorderType==2)
	{
		for(int iii = L; iii < LTot; ++iii)
		{
			spinConf[idxConv(N,iii,0)] = 0;
		}
	}

	// ---> Test the flippable spins
	for(int iii = 0; iii < L; ++iii)
	{
		if(TestFlippable(iii,0))
		{
			flippableSpinLists[dummyNbF] = iii*N;
			flippableSpinPositions[iii*N] = dummyNbF;
			++dummyNbF;
		}
	}
	// > Finally, copy informations over all the layers
	int Pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			Pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[Pos] = spin;
			flippableSpinPositions[spin] = Pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;

	// Count the size of the core
	CoreSize = 0;
	for(int iii = 0; iii < nx + p - 1; ++iii)
	{
		CoreSize += CoreEnd[iii] - CoreStart[iii] + 1;
	}
} */

void SysConf::SetInitialPartAntiFerro()
{
	int dummyNbF = 0;

	CoreStart.resize(nx+p-1,-1);
	CoreEnd.resize(nx+p-1,-1);

	int CoreStartCount = 0;

	// ---> Set borders
	for(int iii = L; iii < LTot; ++iii)
	{
		spinConf[idxConv(N,iii,0)] = pow(-1,iii);
	}

	if(ny>1)
	{
		// Each sector corresponds to a staggered corner
		vector<int> sectorSides(6,-1);
// 		vector<int> sectorSize(6,-1);
// 		int totalSectors = 0;

		// Sector sides
		sectorSides[5] = ny/2;
		sectorSides[0] = ny - sectorSides[5];
		sectorSides[1] = p  - sectorSides[0];
		sectorSides[2] = nx - sectorSides[1];
		sectorSides[3] = ny - sectorSides[2];
		sectorSides[4] = p -  sectorSides[3];

		// Sector SIZES
// 		for (int iii = 0; iii < 6; ++iii)
// 		{
// 			sectorSize[iii] = (sectorSides[iii]*(sectorSides[iii]+1))/2;
// 			totalSectors += sectorSize[iii];
// 		}

		// Initialization
		int spinPos = 0;
		int spinNeigh = 0;
		int firstPos = sectorSides[0]-1;
		int lineLength = 1;
		int counter = 0;

		// ---> First sector
		for(int iii = 0; iii < sectorSides[0];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
			}

			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
			++CoreStartCount;

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,1)];
		}

		// ---> Second sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < sectorSides[1];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,3)];
			}
			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
			++CoreStartCount;
			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			++CoreStartCount;


			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,2)];
		}

		--CoreStartCount;
		// ---> Third sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < sectorSides[2];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
// 				sectors[idxConv(sectorSize,2,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
			}

			CoreStart[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			++CoreStartCount;

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,3)];
		}

		// ---> Fourth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii <sectorSides[3];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
// 				sectors[idxConv(sectorSize,3,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,2)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,5)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,2)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,4)];
		}

		// ---> Fifth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < sectorSides[4];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
// 				sectors[idxConv(sectorSize,4,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,3)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,0)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,3)];

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,4)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,5)];
		}

		++CoreStartCount;

		// ---> Sixth sector
		lineLength = 1;
		counter = 0;
		firstPos = spinPos;
		for(int iii = 0; iii < sectorSides[5];++iii)
		{
			spinPos = firstPos;
			for(int jjj = 0; jjj < lineLength; ++jjj)
			{
				// Fill a line
// 				sectors[idxConv(sectorSize,5,counter)] = spinPos;
				spinNeigh = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
				spinConf[idxConv(N,spinPos,0)] = - spinConf[idxConv(N,spinNeigh,0)];

				++counter;
				spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,1)];
			}

			--CoreStartCount;
			CoreEnd[CoreStartCount] = neighboursTable[idxConv(NbOfNeights,spinPos,4)];

			// Go to next line
			++lineLength;
			firstPos = neighboursTable[idxConv(NbOfNeights,firstPos,0)];
		}

		int spinCompare = 0;
		int defaultSpin = (int)spinConf[idxConv(N,ny/2-1,0)];

		// ---> Set other spins
		for(int pos = 0; pos < L; ++pos)
		{
			if(spinConf[idxConv(N,pos,0)]==0)
			{
				// ---> Then we have an empty cell
				spinConf[idxConv(N,pos,0)] = defaultSpin;
				for(int iii = 0; iii  < NbOfNeights; ++iii)
				{
					spinCompare = neighboursTable[idxConv(NbOfNeights,pos,iii)];

					// ---> Search for neighbour spins that are different
					if((int)spinConf[idxConv(N,spinCompare,0)]==-defaultSpin)
					{
						spinConf[idxConv(N,pos,0)] = defaultSpin;
						break;
					}
	//
	// 				// ---> If none where found, we have to flip this spin
					spinConf[idxConv(N,pos,0)] = -defaultSpin;
				}
			}
		}
	}
	else
	{
		// There's a easier way ...
		int posAbove = 0;
		int pos = 0;
		int baseHeight = 0;

		for(int iii = 0; iii < nx; ++iii)
		{
			baseHeight = max(p - iii,0);
			posAbove = L + iii + 1;

			// Spins above the horizontal dimer
			for(int jjj = 0; jjj < p - baseHeight; ++jjj)
			{
				pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
				spinConf[idxConv(N,pos,0)] = -spinConf[idxConv(N,posAbove,0)];
				posAbove = pos;
			}

			// Spin at the horizontal dimer
			pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
			spinConf[idxConv(N,pos,0)] = spinConf[idxConv(N,posAbove,0)];
			posAbove = pos;

			// Spins below the horizontal dimer
			for(int jjj = p - baseHeight + 1; jjj < p; ++jjj)
			{
				pos = neighboursTable[idxConv(NbOfNeights,posAbove,2)];
				spinConf[idxConv(N,pos,0)] = -spinConf[idxConv(N,posAbove,0)];
				posAbove = pos;
			}
		}
	}

	if(BorderType==0||BorderType==2)
	{
		for(int iii = L; iii < LTot; ++iii)
		{
			spinConf[idxConv(N,iii,0)] = 0;
		}
	}

	// ---> Test the flippable spins
	for(int iii = 0; iii < L; ++iii)
	{
		if(TestFlippable(iii,0))
		{
			flippableSpinLists[dummyNbF] = iii*N;
			flippableSpinPositions[iii*N] = dummyNbF;
			++dummyNbF;
		}
	}
	// > Finally, copy informations over all the layers
	int Pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			Pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[Pos] = spin;
			flippableSpinPositions[spin] = Pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;

	// Count the size of the core
	CoreSize = 0;
	for(int iii = 0; iii < nx + p - 1; ++iii)
	{
		CoreSize += CoreEnd[iii] - CoreStart[iii] + 1;
	}
}

void SysConf::SetInitialPartStar3()
{
	int dummyNbF = 0;
	int firstSpin = 0;
	int spinPos = firstSpin;

	int idx = 0;

	for(int iii = 1; iii < p; ++iii)
	{
		idx = iii - 1;
		for(int jjj = 0; jjj <rowSize[iii] - 1; ++jjj)
		{
			if(idx%3==0)
				spinConf[idxConv(N,spinPos,0)] = -1;
			else
				spinConf[idxConv(N,spinPos,0)] = 1;
			++idx;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
		}
		firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,2)];
		spinPos = firstSpin;
	}

	idx = p - 1;
	for(int jjj = 0; jjj <rowSize[p] - 1; ++jjj)
	{
		if(idx%3==0)
			spinConf[idxConv(N,spinPos,0)] = -1;
		else
			spinConf[idxConv(N,spinPos,0)] = 1;
		++idx;
		spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
	}

	firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,3)];
	spinPos = firstSpin;

	idx = p - 1;

	for(int iii = 1; iii < nx; ++iii)
	{
		idx = p - 1 + 2*iii;
		for(int jjj = 0; jjj <rowSize[iii+p] - 1; ++jjj)
		{
			if(idx%3==0)
				spinConf[idxConv(N,spinPos,0)] = -1;
			else
				spinConf[idxConv(N,spinPos,0)] = 1;
			++idx;
			spinPos = neighboursTable[idxConv(NbOfNeights,spinPos,4)];
		}
		firstSpin = neighboursTable[idxConv(NbOfNeights,firstSpin,3)];
		spinPos = firstSpin;
	}

	bool IsFlippable = false;
	int flippableCandidate = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		flippableCandidate = idxConv(N,0,iii);
		IsFlippable = TestFlippable(flippableCandidate,0);
		if(IsFlippable)
		{
			flippableSpinLists[dummyNbF] = flippableCandidate*N;
			flippableSpinPositions[flippableCandidate*N] = dummyNbF;
			++dummyNbF;
		}
	}

// 	> Finally, copy informations over all the layers
	int pos = 0;
	int spin = 0;
	for(int iii = 1; iii < N; ++iii)
	{
		for(int kkk = 0; kkk < LTot; ++kkk)
		{
			spinConf[idxConv(N,kkk,iii)] = spinConf[idxConv(N,kkk,0)];
		}
		for(int jjj = 0; jjj < dummyNbF; ++jjj)
		{
			pos = iii*dummyNbF + jjj;
			spin = flippableSpinLists[jjj] + iii;
			flippableSpinLists[pos] = spin;
			flippableSpinPositions[spin] = pos;
		}
	}

	TotalOldNbF = N*dummyNbF;
	TotalNewNbF = N*dummyNbF;
}

// *******************************************************************
// RNG setup
// *******************************************************************
void SysConf::SetRNG(double Seed)
{
	m_rng.seed(Seed);
}

void SysConf::ExportRNG(const char outfile[])
{
	ofstream out(outfile,ios::trunc);
	out << m_rng;
	out.close();
}

void SysConf::ImportRNG(const char infile[])
{
	ifstream inp(infile);
	inp >> m_rng;
	inp.close();
}

void SysConf::PrintRNG_DEBUG(int nnn)
{
	random::uniform_real_distribution<> distUnif(0,1);
	for(int iii = 0; iii < nnn; ++iii)
		cout << distUnif(m_rng) << endl;
}

// *******************************************************************
// MC methods - Common
// *******************************************************************

// bool SysConf::TestFlippable(int spin, int layer)
// {
// 	int test = 0;
// 	int index = 0;
// 	for(int iii = 0; iii < NbOfNeights;++iii)
// 	{
// 		index = idxConv(NbOfNeights,spin,iii);
// 		test += spinConf[idxConv(N,neighboursTable[index],layer)];
// 	}
//
// 	return (test==0);
// };

// bool SysConf::TestFlippable(int spin, int layer)
// {
// 	int test = 0;
// 	int index = 0;
//
// 	bool output = true;
// 	int numberOfZeros = 0;
//
// 	if(BorderType!=0)
// 	{
// 		for(int iii = 0; iii < NbOfNeights;++iii)
// 		{
// 			index = idxConv(NbOfNeights,spin,iii);
// 			test += spinConf[idxConv(N,neighboursTable[index],layer)];
// 		}
//
// 		output = (test==0);
// 	}
// 	else
// 	{
// 		for(int iii = 0; iii < NbOfNeights;++iii)
// 		{
// 			index = idxConv(NbOfNeights,spin,iii);
// 			neighsSpin[iii] = spinConf[idxConv(N,neighboursTable[index],layer)];
//
// 			test += neighsSpin[iii];
//
// 			if(neighsSpin[iii]==0)
// 			{
// 				++numberOfZeros;
// 			}
// 		}
//
// 		// Have to do a detailed test ...
// 		for(int iii = 1; iii < NbOfNeights; ++iii)
// 		{
// 			if(neighsSpin[iii]*neighsSpin[iii-1]==1)
// 			{
// 				output = false;
// 				break;
// 			}
// 		}
// 		if(neighsSpin[0]*neighsSpin[NbOfNeights-1]==1)
// 		{
// 			output = false;
// 		}
// 	}
//
//
// 	return output;
// };

bool SysConf::TestFlippable(int spin, int layer)
{
	int test = 0;
	int index = 0;

	bool output = true;

#if BORDER_TYPE == 4
// ANTI-PERIODIC BORDERS
	for(int iii = 0; iii < NbOfNeights;++iii)
	{
		index = idxConv(NbOfNeights,spin,iii);
		test += weightTable[index]*spinConf[idxConv(N,neighboursTable[index],layer)];
	}

	output = (test==0);
#else
// NORMAL BORDERS
	int numberOfZeros = 0;

	switch (BorderType)
	{
		case 0:	// Free flips!
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index = idxConv(NbOfNeights,spin,iii);
					neighsSpin[iii] = spinConf[idxConv(N,neighboursTable[index],layer)];
				}

				// Have to do a detailed test ...
				for(int iii = 1; iii < NbOfNeights; ++iii)
				{
					if(neighsSpin[iii]*neighsSpin[iii-1]==1)
					{
						output = false;
						break;
					}
				}
				if(neighsSpin[0]*neighsSpin[NbOfNeights-1]==1)
				{
					output = false;
				}
				break;
			}
		case 2:	// Semi-free flips
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index = idxConv(NbOfNeights,spin,iii);
					neighsSpin[iii] = spinConf[idxConv(N,neighboursTable[index],layer)];

					if(neighsSpin[iii]==0)
					{
						++numberOfZeros;
					}
				}

				if(numberOfZeros%2==0)
				{
					// Have to do a detailed test ...
					for(int iii = 1; iii < NbOfNeights; ++iii)
					{
						if(neighsSpin[iii]*neighsSpin[iii-1]==1)
						{
							output = false;
							break;
						}
					}
					if(neighsSpin[0]*neighsSpin[NbOfNeights-1]==1)
					{
						output = false;
					}
				}
				else
				{
					output = false;
				}
				break;
			}
		default:
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index = idxConv(NbOfNeights,spin,iii);
					test += spinConf[idxConv(N,neighboursTable[index],layer)];
				}

				output = (test==0);

				break;
			}
	}
#endif
	return output;
};

void SysConf::ClusterUpdate()
{
	// > Parameters
	int nbToExchange = 0;
	int posToChange = 0;
	int spinToMove = 0;
	int LastPosition = TotalOldNbF - 1;

	int nnn = (clusterStart) % N;

	for(int iii = 0; iii < clusterSize; ++iii)
	{
		++nnn;
		if(nnn==N)
		{
			nnn = 0;
		}

		nbToExchange = min(nbToRemove[nnn],nbToAdd[nnn]);
		posToChange = 0;
		spinToMove = 0;


		// > Spin conf update
		spinConf[idxConv(N,kkkCluster,nnn)] = -spinConf[idxConv(N,kkkCluster,nnn)];

		// > Spin tables update
		// > We have 3 situations regarding 'flippableSpins'
		//   1) nbToAdd == nbToRemove : no holes, no additional terms
		//   2) nbToAdd < nbToRemove  : we'll have some holes
		//   3) nbToAdd > nbToRemove  : no holes, but we'll have to add some terms
		//
		// > We can exchange the spins over the first 'Min(nbToRemove,nbToAdd)' spins
		//   to add/remove. Then :
		//   *  If we have the case c1), there's nothing more to do
		//   *  Case 2), we have to compact the lists
		//   *  Case 3), we have to add elements

		// > Exchange positions
		// * If we're only adding or removing spins, we won't enter this loop
		//   because nbToExchange == 0

		for(int iii = 0; iii < nbToExchange; ++iii)
		{
			// > Get position to exchange
			posToChange = flippableSpinPositions[toRemove[nnn][iii]];

			// > Exchange the spins
			flippableSpinLists[posToChange] = toAdd[nnn][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toRemove[nnn][iii]] = -1;
			flippableSpinPositions[toAdd[nnn][iii]] = posToChange;
		}

		// > Add elements
		// * If we're on the case 2), we won't enter this loop
		//   because nbToExchange == nbToAdd
		for(int iii = nbToExchange; iii < nbToAdd[nnn]; ++iii)
		{
			// > Update the LastPosition
			++LastPosition;

			// > Add the spin
			flippableSpinLists[LastPosition] = toAdd[nnn][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toAdd[nnn][iii]] = LastPosition;
		}

		// > Remove elements and fill holes
		// * If we're on the case 3), we won't enter this loop
		//   because nbToExchange == nbToRemove
		for(int iii = nbToExchange; iii < nbToRemove[nnn]; ++iii)
		{
			// > Get position to remove
			posToChange = flippableSpinPositions[toRemove[nnn][iii]];

			// > Remove spin
			flippableSpinLists[posToChange] = -1;
			flippableSpinPositions[toRemove[nnn][iii]] = -1;

			if(posToChange < LastPosition)
			{
				// Then we created a hole -> we have to update!
				spinToMove = flippableSpinLists[LastPosition];

				// Move the element
				flippableSpinLists[posToChange] = spinToMove;
				flippableSpinLists[LastPosition] = -1;

				// Update position
				flippableSpinPositions[spinToMove] = posToChange;
			}

			// > Update the LastPosition
			--LastPosition;
		}
	}

	TotalOldNbF = TotalNewNbF;
}

bool SysConf::GetRandomAccept(double inAccept)
{
	random::bernoulli_distribution<> bernoulliMC(inAccept);
	return bernoulliMC(m_rng);
}

int SysConf::GetAbsIndex(int InTotalOldNbF)
{
	random::uniform_int_distribution<> dist(0, InTotalOldNbF - 1);
	return dist(m_rng);
}

void SysConf::TestUpdate(int layer)
{

	int ChosenNeigh = 0;
	int spin = 0;
	bool IsFlippable = false;

	nbToRemove[layer] = 0;
	nbToAdd[layer] = 0;

	spinConf[idxConv(N,kkkCluster,layer)] = -spinConf[idxConv(N,kkkCluster,layer)];
	for(int iii = 0; iii < NbOfNeights; ++iii)
	{

		ChosenNeigh = neighboursTable[idxConv(NbOfNeights,kkkCluster,iii)];

		if(ChosenNeigh<L)
		{
			spin = idxConv(N,ChosenNeigh,layer);
			if(flippableSpinPositions[spin]!=-1)
			{
				// Then the spin is flippable in the old configuration, but won't be
				//	anymore if the new conf is accepted
				// >>> Mark for removal
				toRemove[layer][nbToRemove[layer]] = spin;
				++nbToRemove[layer];
			}
			else
			{
				// Then the spin can't be flipped in the old configuration, but MIGHT
				//	be in the new one
				// >>> Have to test
				IsFlippable = TestFlippable(ChosenNeigh,layer);

				if(IsFlippable)
				{
					toAdd[layer][nbToAdd[layer]] = spin;
					++nbToAdd[layer];
				}
			}
			if(nbToAdd[layer]>3||nbToRemove[layer]>3)
			{
				cout << "Red alert!" << endl;
			}
		}
	}
	spinConf[idxConv(N,kkkCluster,layer)] = -spinConf[idxConv(N,kkkCluster,layer)];
}

void SysConf::RaiseMeasures(int UpdateInterval)
{
	measurementsDone += UpdateInterval;
}


// *******************************************************************
// MC methods - Diagonal = N3
// *******************************************************************
void SysConf::MC_Prepare(MCParameters& inParams)
{
	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	vector<int> dummy(NbOfNeights/2,-1);

	toRemove.resize(N,dummy);
	toAdd.resize(N,dummy);

	IsAccepted = true;
	accept = 0;
	totalAccept = 0;

	Kz = inParams.Kz;
	Kt = inParams.Kt;
	q = inParams.q;
	Vt = inParams.Vt;

	baseProb = 1 - q;

	CalculateE();

	clusterStart = 0;
	clusterEnd = 0;
	clusterSize = 1;

	alpha0 = 0;
	alpha3 = 1;

	shift0 = 0;
	shift3 = 0;

	expZ.resize(NbOfNeights+1,0);
	cout << " --- " << endl;
	for(int iii = 0; iii < NbOfNeights + 1; ++iii)
	{
		expZ[iii]= exp(-Kz*Vt*(iii - NbOfNeights/2));
		cout << expZ[iii] << endl;
	}
	cout << " --- " << endl;
	neighsSpin.resize(NbOfNeights,0);

	minExp = 0;
	maxExp = 0;
}

void SysConf::CalculateE()
{
	energy = Kz * TotalOldNbF * Vt;

	for(int nnn = 0; nnn < N-1; ++nnn)
		for(int iii = 0; iii < L; ++iii)
		{
			energy -= Kt*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn)+1];
		}

	for(int iii = 0; iii < L; ++iii)
	{
		energy -= Kt*spinConf[idxConv(N,iii,0) + N-1]*spinConf[idxConv(N,iii,0)];
	}
}

void SysConf::ClusterIter()
{
	BuildCluster();
	GetAcceptanceCluster();

	IsAccepted = GetRandomAccept(accept);
	if(IsAccepted)
	{
		++totalAccept;
		ClusterUpdate();
		energy += DE;
	}
}

void SysConf::BuildCluster()
{
	// >>>> Get a random spin
	absIndex 	= GetAbsIndex(TotalOldNbF/*,m_rng*/);
	spinToFlip 	= flippableSpinLists[absIndex];
	nnnCluster	= spinToFlip % N;
	kkkCluster 	= spinToFlip / N;
	firstSpin 	= spinConf[spinToFlip];

	// Initialize parameters
	TotalNewNbF 		= TotalOldNbF;

// 	bool 	CanBuildCluster	= true;
// 	bool	IsSpinAdded 	= true;
//
// 	double	addProb		= 0;
//
	clusterStart 		= nnnCluster;
	clusterEnd 		= nnnCluster;
	clusterSize		= 1;

	UpperBorderIsF		= 1;
	LowerBorderIsF		= 1;

	// > Calculate the changes created by choosing this first spin
	TestUpdate(nnnCluster);

	TotalNewNbF += nbToAdd[nnnCluster]-nbToRemove[nnnCluster];

	// > Create cluster
	// First upwards
	GrowCluster(clusterEnd,UpperBorderIsF,1,N,0);

	// Then downwards
	GrowCluster(clusterStart,LowerBorderIsF,-1,-1,N-1);
}

void SysConf::GrowCluster(int& workPoint, int& borderStatus, int step, int border, int borderRedirect)
{
	bool CanBuildCluster 	= true;
	bool	IsSpinAdded 	= true;

	borderStatus		= 1;
	double	addProb		= 0;

	while(CanBuildCluster)
	{
		workPoint += step;

		if(workPoint==border)
		{
			workPoint = borderRedirect;
		}

		if(clusterSize==N)
		{
			// Then the cluster has ended -> Finish
			CanBuildCluster=false;
		}
		else
		{
			if(flippableSpinPositions[idxConv(N,kkkCluster,workPoint)]==-1)
			{	// Then we can't flip the spin
				// Also, the cluster has ended -> Finish
				borderStatus = 0;
				CanBuildCluster = false;
			}
			else
			{
				// Set up would-be new configuration
				TestUpdate(workPoint);

				if(spinConf[idxConv(N,kkkCluster,workPoint)]!=firstSpin)
				{
					// Then we've hit a different spin
					// Also, the cluster has ended -> Finish
					CanBuildCluster = false;
				}
				else
				{
					// Add spin to cluster with probability addProb
					addProb = baseProb*min(1,expZ[nbToAdd[workPoint]-nbToRemove[workPoint]+3]);
					IsSpinAdded = GetRandomAccept(addProb);
					if(IsSpinAdded)
					{
						// Update spin configuration
						++clusterSize;
						TotalNewNbF += nbToAdd[workPoint]-nbToRemove[workPoint];
					}
					else
					{
						// We didn't want this new spin

						// Also, the cluster has ended -> Finish
						CanBuildCluster = false;
					}
				}
			}
		}
	}
}

void SysConf::GetAcceptanceCluster()
{
	double upperSpin = (int)spinConf[idxConv(N,kkkCluster,clusterEnd)];
	double lowerSpin = (int)spinConf[idxConv(N,kkkCluster,clusterStart)];

	DE = Kz*Vt*(TotalNewNbF - TotalOldNbF) +
				2*Kt*firstSpin*(upperSpin + lowerSpin);

	// Extra weight given by the cluster update
	double weight = expZ[nbToAdd[nnnCluster] - nbToRemove[nnnCluster]+3]*exp(-2*Kt*firstSpin*(upperSpin + lowerSpin));
								// Spin flip - OK
	weight = weight*pow(1-(1-q)*min(1,expZ[nbToAdd[clusterEnd] - nbToRemove[clusterEnd] + 3]),-UpperBorderIsF*firstSpin*upperSpin);
								// Cluster boreder m
	weight = weight*pow(1-(1-q)*min(1,expZ[nbToAdd[clusterStart] - nbToRemove[clusterStart] + 3]),-LowerBorderIsF*firstSpin*lowerSpin);
								// Cluster boreder m'

	// Calculate the acceptance
	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);
}

// *******************************************************************
// MC methods - Diagonal = N0
// *******************************************************************

void SysConf::MC_Prepare_Zero(MCParameters& inParams)
{
	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	vector<int> dummy(NbOfNeights/2,-1);
	toRemove.resize(N,dummy);
	toAdd.resize(N,dummy);

	DeltaN0.resize(N,0);

	IsAccepted = true;
	accept = 0;
	totalAccept = 0;

	Kz = inParams.Kz;

	Kt = inParams.Kt;
	q = inParams.q;
	Vt = inParams.Vt;

	baseProb = 1 - q;

	CalculateE_Zero();

	clusterStart = 0;
	clusterEnd = 0;
	clusterSize = 1;

	alpha0 = 1;
	alpha3 = 0;

	shift0 = 0;
	shift3 = 0;

	expZ.resize(NbOfNeights+1,0);
	for(int iii = 0; iii < NbOfNeights + 1; ++iii)
	{
		expZ[iii]= exp(-Kz*Vt*(iii - NbOfNeights/2));
	}

	neighsSpin.resize(NbOfNeights,0);
}

void SysConf::CalculateE_Zero()
{
	int N0 = GetN0();
	energy = Kz * N0 * Vt;

	for(int nnn = 0; nnn < N-1; ++nnn)
		for(int iii = 0; iii < L; ++iii)
		{
			energy -= Kt*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn)+1];
		}

	for(int iii = 0; iii < L; ++iii)
	{
		energy -= Kt*spinConf[idxConv(N,iii,0) + N-1]*spinConf[idxConv(N,iii,0)];
	}
}

void SysConf::ClusterIter_Zero()
{
	BuildCluster_Zero();
	GetAcceptanceCluster_Zero();

	IsAccepted = GetRandomAccept(accept);
	if(IsAccepted)
	{
		++totalAccept;
		ClusterUpdate();
		energy += DE;
	}
}

void SysConf::BuildCluster_Zero()
{
	// >>>> Get a random spin
	absIndex 		= GetAbsIndex(TotalOldNbF/*,m_rng*/);
	spinToFlip 	= flippableSpinLists[absIndex];
	nnnCluster	= spinToFlip % N;
	kkkCluster 	= spinToFlip / N;
	firstSpin 	= spinConf[spinToFlip];

	// Initialize parameters
	TotalNewNbF 		= TotalOldNbF;

// 	bool 	CanBuildCluster	= true;
// 	bool	IsSpinAdded 	= true;
//
// 	double	addProb		= 0;
//
	clusterStart 		= nnnCluster;
	clusterEnd 		= nnnCluster;
	clusterSize		= 1;

	UpperBorderIsF		= 1;
	LowerBorderIsF		= 1;

	// > Calculate the changes created by choosing this first spin
	TestUpdate(nnnCluster);
	Update_DeltaN0(nnnCluster);

	Total_DeltaN0 = DeltaN0[nnnCluster];
	TotalNewNbF += nbToAdd[nnnCluster]-nbToRemove[nnnCluster];

	// > Create cluster

	// Then downwards
	GrowCluster_Zero(clusterStart,LowerBorderIsF,-1,-1,N-1);

	// First upwards
	GrowCluster_Zero(clusterEnd,UpperBorderIsF,1,N,0);
}

void SysConf::GrowCluster_Zero(int& workPoint, int& borderStatus, int step, int border, int borderRedirect)
{
	bool CanBuildCluster 	= true;
	bool	IsSpinAdded 	= true;

	borderStatus		= 1;
	double	addProb		= 0;

	while(CanBuildCluster)
	{
		workPoint += step;

		if(workPoint==border)
		{
			workPoint = borderRedirect;
		}

		if(clusterSize==N)
		{
			// Then the cluster has ended -> Finish
			CanBuildCluster=false;
		}
		else
		{
			if(flippableSpinPositions[idxConv(N,kkkCluster,workPoint)]==-1)
			{	// Then we can't flip the spin
				// Also, the cluster has ended -> Finish
				borderStatus = 0;
				CanBuildCluster = false;
			}
			else
			{
				// Set up would-be new configuration
				TestUpdate(workPoint);
				Update_DeltaN0(workPoint);

				if(spinConf[idxConv(N,kkkCluster,workPoint)]!=firstSpin)
				{
					// Then we've hit a different spin
					// Also, the cluster has ended -> Finish
					CanBuildCluster = false;
				}
				else
				{
					// Add spin to cluster with probability addProb
					addProb = baseProb*min(1,exp(-Kz*Vt*DeltaN0[workPoint]));
					IsSpinAdded = GetRandomAccept(addProb);
					if(IsSpinAdded)
					{
						// Update spin configuration
						++clusterSize;
						TotalNewNbF += nbToAdd[workPoint]-nbToRemove[workPoint];
						Total_DeltaN0 += DeltaN0[workPoint];
					}
					else
					{
						// We didn't want this new spin

						// Also, the cluster has ended -> Finish
						CanBuildCluster = false;
					}
				}
			}
		}
	}
}

void SysConf::GetAcceptanceCluster_Zero()
{
	double upperSpin = (int)spinConf[idxConv(N,kkkCluster,clusterEnd)];
	double lowerSpin = (int)spinConf[idxConv(N,kkkCluster,clusterStart)];

	DE = Kz*Vt*Total_DeltaN0 +
				2*Kt*firstSpin*(upperSpin + lowerSpin);

	// Extra weight given by the cluster update
	double weight = exp(-Kz*Vt*DeltaN0[nnnCluster]-2*Kt*firstSpin*(upperSpin + lowerSpin));
								// Spin flip - OK
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*DeltaN0[clusterEnd])),-UpperBorderIsF*firstSpin*upperSpin);
								// Cluster boreder m
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*DeltaN0[clusterStart])),-LowerBorderIsF*firstSpin*lowerSpin);
								// Cluster boreder m'

	// Calculate the acceptance
	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);
}

void SysConf::Update_DeltaN0(int layer)
{
	int ChosenNeigh = 0;
	int localMag = 0;
	int index = 0;

	DeltaN0[layer] = 0;

	for(int iii = 0; iii < NbOfNeights; ++iii)
	{
		localMag = 0;
		index = idxConv(NbOfNeights,kkkCluster,iii);
		ChosenNeigh = neighboursTable[index];

		if(ChosenNeigh<L)
		{
			localMag = abs(GetLocalField(ChosenNeigh,layer));

			if(localMag==6)
			{
				// Then the spin has zero frustrated relations in the old configuration,
				//   but that won't be the case anymore if the new conf is accepted
				// >>> Remove 1 from \DeltaN0
				--DeltaN0[layer];
			}
// 			else if(localMag==4)
// 			{
// 				// Then the spin has 1 frustrated relation in the old configuration,
// 				//   but MIGHT lose it in the new one
// 				// >>> Have to test
//
// 				spinConf[idxConv(N,kkkCluster,layer)] = -spinConf[idxConv(N,kkkCluster,layer)];
// 				dummyMag = abs(GetLocalField(ChosenNeigh,layer));
// 				if(dummyMag==6)
// 				{
// 					++DeltaN0[layer];
// 				}
// 				spinConf[idxConv(N,kkkCluster,layer)] = spinConf[idxConv(N,kkkCluster,layer)];
// 			}
			else if(localMag==4)
			{
				// Then the spin has 1 frustrated relation in the old configuration,
				//   but MIGHT lose it in the new one
				// >>> Have to test
			#if BORDER_TYPE == 4
			// ANTI-PERIODIC BORDERS
				if(spinConf[idxConv(N,kkkCluster,layer)]==weightTable[index]*spinConf[idxConv(N,ChosenNeigh,layer)])
				{
					++DeltaN0[layer];
				}
			#else
				if(spinConf[idxConv(N,kkkCluster,layer)]==spinConf[idxConv(N,ChosenNeigh,layer)])
				{
					++DeltaN0[layer];
				}
			#endif
			}
		}
	}

}

// *******************************************************************
// MC methods - Diagonal = Mixed
// *******************************************************************

void SysConf::MC_Prepare_Mixed(MCParameters& inParams)
{
	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	vector<int> dummy(NbOfNeights/2,-1);
	toRemove.resize(N,dummy);
	toAdd.resize(N,dummy);

	if(SimType.compare("Part")==0)
	{
		equivPart.resize(nx*ny*N,0);
		horizontalDimerPos.resize(ny*N,0);
		if(ny==1)
		{
			spinChain.resize((nx+p)*N,1);
		}
	}

	DeltaN0.resize(N,0);

	IsAccepted = true;
	accept = 0;
	totalAccept = 0;

	Kz = inParams.Kz;

	Kt = inParams.Kt;
	q = inParams.q;
	Vt = inParams.Vt;

	double dummyAlpha = sqrt(inParams.alpha0*inParams.alpha0 + inParams.alpha3*inParams.alpha3);

	alpha0 = inParams.alpha0/dummyAlpha;
	alpha3 = inParams.alpha3/dummyAlpha;

	shift0 = inParams.shift0;
	shift3 = inParams.shift3;

	baseProb = 1 - q;

	CalculateE_Mixed();

	clusterStart = 0;
	clusterEnd = 0;
	clusterSize = 1;

	expZ.resize(NbOfNeights+1,0);
	for(int iii = 0; iii < NbOfNeights + 1; ++iii)
	{
		expZ[iii]= exp(-Kz*(alpha3*Vt + shift3)*(iii - NbOfNeights/2));
	}

	neighsSpin.resize(NbOfNeights,0);
	DimerDensityProfile.resize(N*L,0);
	N3_Map.resize(N*L,0);

	Jindex.resize(3,1);
	complex<double> dummyCmplx = 1i*2.*M_PI/3.;
	Jindex[1] = exp(dummyCmplx);
	dummyCmplx = 1i*4.*M_PI/3.;
	Jindex[2] = exp(dummyCmplx);

	dummySubs.resize(3);
}

void SysConf::CalculateE_Mixed()
{
	int N0 = GetN0();
	int N3 = TotalOldNbF;

	energy = Kz * ((Vt  * alpha0 + shift0) * N0 + (Vt  * alpha3 + shift3) * N3);

	for(int nnn = 0; nnn < N-1; ++nnn)
	{
		for(int iii = 0; iii < L; ++iii)
		{
			energy -= Kt*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn)+1];
		}
	}
	for(int iii = 0; iii < L; ++iii)
	{
		energy -= Kt*spinConf[idxConv(N,iii,0) + N-1]*spinConf[idxConv(N,iii,0)];
	}
}

void SysConf::ClusterIter_Mixed()
{
	BuildCluster_Mixed();
	GetAcceptanceCluster_Mixed();

	IsAccepted = GetRandomAccept(accept);
	if(IsAccepted)
	{
		++totalAccept;
		ClusterUpdate();
		energy += DE;
	}
}

void SysConf::BuildCluster_Mixed()
{
	// >>>> Get a random spin
	absIndex 	= GetAbsIndex(TotalOldNbF/*,m_rng*/);
	spinToFlip 	= flippableSpinLists[absIndex];
	nnnCluster	= spinToFlip % N;
	kkkCluster 	= spinToFlip / N;
	firstSpin 	= spinConf[spinToFlip];

	// Initialize parameters
	TotalNewNbF 		= TotalOldNbF;

// 	bool 	CanBuildCluster	= true;
// 	bool	IsSpinAdded 	= true;
//
// 	double	addProb		= 0;
//
	clusterStart 		= nnnCluster;
	clusterEnd 		= nnnCluster;
	clusterSize		= 1;

	UpperBorderIsF		= 1;
	LowerBorderIsF		= 1;

	// > Calculate the changes created by choosing this first spin
	TestUpdate(nnnCluster);

	if(abs(alpha0)>10E-8||abs(shift0)>10E-8)
	{
		Update_DeltaN0(nnnCluster);
		Total_DeltaN0 = DeltaN0[nnnCluster];
	}

	TotalNewNbF += nbToAdd[nnnCluster]-nbToRemove[nnnCluster];

	// > Create cluster
	// First upwards
	GrowCluster_Mixed(clusterEnd,UpperBorderIsF,1,N,0);

	// Then downwards
	GrowCluster_Mixed(clusterStart,LowerBorderIsF,-1,-1,N-1);
}

void SysConf::GrowCluster_Mixed(int& workPoint, int& borderStatus, int step, int border, int borderRedirect)
{
	bool CanBuildCluster 	= true;
	bool	IsSpinAdded 	= true;

	borderStatus		= 1;
	double	addProb		= 0;

	while(CanBuildCluster)
	{
		workPoint += step;

		if(workPoint==border)
		{
			workPoint = borderRedirect;
		}

		if(clusterSize==N)
		{
			// Then the cluster has ended -> Finish
			CanBuildCluster=false;
		}
		else
		{
			if(flippableSpinPositions[idxConv(N,kkkCluster,workPoint)]==-1)
			{	// Then we can't flip the spin
				// Also, the cluster has ended -> Finish
				borderStatus = 0;
				CanBuildCluster = false;
			}
			else
			{
				// Set up would-be new configuration
				TestUpdate(workPoint);

				if(abs(alpha0)>10E-8||abs(shift0)>10E-8)
				{
					Update_DeltaN0(workPoint);
				}

				if(spinConf[idxConv(N,kkkCluster,workPoint)]!=firstSpin)
				{
					// Then we've hit a different spin
					// Also, the cluster has ended -> Finish
					CanBuildCluster = false;
				}
				else
				{
					// Add spin to cluster with probability addProb
					addProb = baseProb*min(1,exp(-Kz*(Vt*alpha0+shift0)*DeltaN0[workPoint])*
											expZ[nbToAdd[workPoint]-nbToRemove[workPoint]+3]);
					IsSpinAdded = GetRandomAccept(addProb);
					if(IsSpinAdded)
					{
						// Update spin configuration
						++clusterSize;
						TotalNewNbF += nbToAdd[workPoint]-nbToRemove[workPoint];
						if(abs(alpha0)>10E-8||abs(shift0)>10E-8)
						{
							Total_DeltaN0 += DeltaN0[workPoint];
						}
					}
					else
					{
						// We didn't want this new spin

						// Also, the cluster has ended -> Finish
						CanBuildCluster = false;
					}
				}
			}
		}
	}
}

void SysConf::GetAcceptanceCluster_Mixed()
{
	double upperSpin = (int)spinConf[idxConv(N,kkkCluster,clusterEnd)];
	double lowerSpin = (int)spinConf[idxConv(N,kkkCluster,clusterStart)];


	// Original codes !
	// N0 :

	//	DE = Kz*Vt*Total_DeltaN0 +
	//				2*Kt*firstSpin*(upperSpin + lowerSpin);
	//
	//	// Extra weight given by the cluster update
	//	double weight = exp(-Kz*Vt*DeltaN0[nnnCluster]-2*Kt*firstSpin*(upperSpin + lowerSpin));
	//								// Spin flip - OK
	//	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*DeltaN0[clusterEnd])),-UpperBorderIsF*firstSpin*upperSpin);
	//								// Cluster boreder m
	//	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*Vt*DeltaN0[clusterStart])),-LowerBorderIsF*firstSpin*lowerSpin);
	//								// Cluster boreder m'
	//
	//	// Calculate the acceptance
	//	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);

	// N3 :

	//	DE = Kz*Vt*(TotalNewNbF - TotalOldNbF) +
	//				2*Kt*firstSpin*(upperSpin + lowerSpin);
	//
	//	// Extra weight given by the cluster update
	//	double weight = expZ[nbToAdd[nnnCluster] - nbToRemove[nnnCluster]+3]*exp(-2*Kt*firstSpin*(upperSpin + lowerSpin));
	//								// Spin flip - OK
	//	weight = weight*pow(1-(1-q)*min(1,expZ[nbToAdd[clusterEnd] - nbToRemove[clusterEnd] + 3]),-UpperBorderIsF*firstSpin*upperSpin);
	//								// Cluster boreder m
	//	weight = weight*pow(1-(1-q)*min(1,expZ[nbToAdd[clusterStart] - nbToRemove[clusterStart] + 3]),-LowerBorderIsF*firstSpin*lowerSpin);
	//								// Cluster boreder m'
	//
	//	// Calculate the acceptance
	//	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);

	// Extra weight given by the cluster update
	DE = Kz*((Vt*alpha0 + shift0)*Total_DeltaN0 +
				(Vt*alpha3 + shift3)*(TotalNewNbF - TotalOldNbF)) +
				2*Kt*firstSpin*(upperSpin + lowerSpin);

	double weight = expZ[nbToAdd[nnnCluster]-nbToRemove[nnnCluster]+3]*
						exp(-Kz*(Vt*alpha0+shift0)*DeltaN0[nnnCluster]-
						2*Kt*firstSpin*(upperSpin + lowerSpin));
								// Spin flip - OK
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*(Vt*alpha0+shift0)*DeltaN0[clusterEnd])*
									expZ[nbToAdd[clusterEnd]-nbToRemove[clusterEnd]+3]),
									-UpperBorderIsF*firstSpin*upperSpin);
								// Cluster boreder m
	weight = weight*pow(1-(1-q)*min(1,exp(-Kz*(Vt*alpha0+shift0)*DeltaN0[clusterStart])*
									expZ[nbToAdd[clusterStart]-nbToRemove[clusterStart]+3]),
									-LowerBorderIsF*firstSpin*lowerSpin);
								// Cluster boreder m'

	// Calculate the acceptance
	accept = min(1,weight*(double)TotalOldNbF/TotalNewNbF);
}

// *******************************************************************
// Order Parameters / getters
// *******************************************************************
double SysConf::GetEnergy()
{
	return energy;
}

double SysConf::GetClusterSize()
{
	return clusterSize;
}

//double SysConf::GetMagnetization()
//{
//	double Mag = 0;
//
//	for(int iii = 0; iii < N; ++iii)
//	{
//		for(int jjj = 0; jjj < L; ++jjj)
//		{
//			Mag += (double)spinConf[idxConv(N,jjj,iii)];
//		}
//	}
//	return pow(Mag/(L*N),2);
//}

double SysConf::GetMagnetization()
{
	double Mag = 0;
	double outMag = 0;

	for(int iii = 0; iii < N; ++iii)
	{
		Mag = 0;
		for(int jjj = 0; jjj < L; ++jjj)
		{
			Mag += (double)spinConf[idxConv(N,jjj,iii)];
		}
		outMag += pow(Mag/L,2)/N;
	}
	return sqrt(outMag);
}


double SysConf::GetCoreMagnetization()
{
	double Mag = 0;

	for(int iii = 0; iii < N; ++iii)
	{
		for(int jjj = 0; jjj < nx + p - 1; ++jjj)
		{
			for(int kkk = CoreStart[jjj]; kkk <= CoreEnd[jjj]; ++kkk)
			{
				Mag += (double)spinConf[idxConv(N,kkk,iii)];
			}
		}
	}
	return pow(Mag/(CoreSize*N),2);
}

ulong SysConf::GetTotalAccept()
{
	return totalAccept;
}

void SysConf::GetNf(vector<double> &Ncount)
{
	int test_Ncount = 0;

	int spinPos = 0;
	int spinValue = 0;

	int neighPos = 0;
	int neighValue = 0;

	int index = 0;
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int lll = 0; lll < L; ++lll)
		{
			// Count the number of dimers
			test_Ncount	= 0;
			spinPos 	= idxConv(N,lll,nnn);
			spinValue 	= (int)spinConf[spinPos];

			if(spinValue!=0)
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index 		= idxConv(NbOfNeights,lll,iii);
					neighPos	= idxConv(N,neighboursTable[index],nnn);

				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					neighValue	= weightTable[index]*(int)spinConf[neighPos];
				#else
					neighValue	= (int)spinConf[neighPos];
				#endif

					if(neighValue==spinValue)
					{
						++test_Ncount;
					}
				}

				Ncount[test_Ncount] 	+= 1./N;
			}
		}
	}

//	for(int lll = 0; lll < L; ++lll)
//	{
//		for(int nnn = 0; nnn < N; ++nnn)
//		{
//			test = 0;
//			for(int iii = 0; iii < NbOfNeights;++iii)
//			{
//				index = idxConv(NbOfNeights,lll,iii);
//				test += (int)spinConf[idxConv(N,neighboursTable[index],nnn)];
//			}
//
//			switch(test)
//			{
//				case 0:
//					Ncount[3] += 1;
//					break;
//				case 2:
//					Ncount[2]  += 1;
//					break;
//				case -2:
//					Ncount[2]  += 1;
//					break;
//				case 4:
//					Ncount[1]  += 1;
//					break;
//				case -4:
//					Ncount[1]  += 1;
//					break;
//				case 6:
//					Ncount[0]  += 1;
//					break;
//				case -6:
//					Ncount[0]  += 1;
//					break;
//			}
//		}
//	}
//
//	for(int iii = 0; iii < 4; ++iii)
//		Ncount[iii] = Ncount[iii]/N;
}

int SysConf::GetN0()
{
	bool test = true;
	int output = 0;

	for(int lll = 0; lll < L; ++lll)
	{
		for(int nnn = 0; nnn < N; ++nnn)
		{
			test = HasZeroDimers(lll,nnn);
			if(test)
				++output;
		}
	}

	return output;
}

int SysConf::HasZeroDimers(int spin, int layer)
{
	int test = 0;
	bool output = false;
	int index = 0;
	int numberOfZeros = 0;

	for(int iii = 0; iii < NbOfNeights;++iii)
	{
		index = idxConv(NbOfNeights,spin,iii);
	#if BORDER_TYPE == 4
	// ANTI-PERIODIC BORDERS
		test += weightTable[index]*spinConf[idxConv(N,neighboursTable[index],layer)];
	#else
		test += spinConf[idxConv(N,neighboursTable[index],layer)];
	#endif

		if(spinConf[idxConv(N,neighboursTable[index],layer)]==0)
		{
			++numberOfZeros;
		}
	}

	if(abs(test)== 6 - numberOfZeros)
	{
		output = true;
	}

	return output;
};

int SysConf::GetLocalField(int lll, int nnn)
{
	int output = 0;
	int index = 0;

	for(int iii = 0; iii < NbOfNeights;++iii)
	{
		index = idxConv(NbOfNeights,lll,iii);
	#if BORDER_TYPE == 4
	// ANTI-PERIODIC BORDERS
		output += weightTable[index]*(int)spinConf[idxConv(N,neighboursTable[index],nnn)];
	#else
		output += (int)spinConf[idxConv(N,neighboursTable[index],nnn)];
	#endif
	}

	return output;
}

void SysConf::BuildEquivalentSpinChain()
{
	int pos = 0;
	int index = 0;
//	spinChain.resize((nx+p)*N,1);	// spinChain[(nx+p)*iii + jjj] = stack iii, spin jjj
//	equivPart.resize(nx*N,0);	// equivPart[(nx)*iii + jjj] = stack iii, spin jjj

	for(int nnn = 0; nnn < N; ++nnn)
	{
		// Set the first nx spins as positive
		for(int jjj = 0; jjj < nx + p; ++jjj)
		{
			spinChain[idxConv(nx+p,nnn,jjj)] = -1;
		}

		// Recover equivalent partition
		for(int jjj = 0; jjj < nx;++jjj)
		{
			equivPart[idxConv(nx,nnn,jjj)] = p;
			pos = L + jjj + 1;

			for(int kkk = 0; kkk < p; ++kkk)
			{

				// The neighbour below SHOULD be neighboursTable[L][2] (see hex offset)
				if(spinConf[idxConv(N,pos,nnn)]==spinConf[idxConv(N,neighboursTable[idxConv(6,pos,2)],nnn)])
				{
					// Then we got the horizontal dimer / frustrated relations
					break;
				}
				else
				{
//					--equivPart[idxConv(nx,nnn,jjj)];
					--equivPart[idxConv(nx,jjj,nnn)];
					pos = neighboursTable[idxConv(6,pos,2)];
				}
			}
		}

		// Set up spin chain
		for(int jjj = 0; jjj < nx; ++jjj)
		{
//			index = nx - 1 + equivPart[idxConv(nx,nnn,jjj)] - jjj;
			index = nx - 1 + equivPart[idxConv(N,jjj,nnn)] - jjj;
			spinChain[idxConv(nx+p,nnn,index)] = 1;
		}
	}
};

double SysConf::CalculateStagMag()
{
	double stagMag = 0;

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int jjj = 0; jjj < nx + p; ++jjj)
		{
			stagMag += pow(-1,jjj)*spinChain[idxConv(nx+p,nnn,jjj)];
		}
	}

	return abs(stagMag)/((nx+p)*N);
};

double SysConf::CalculateChainEnergy()
{
	double chainEnergy = 0;

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int jjj = 0; jjj < nx + p - 1; ++jjj)
		{
			chainEnergy += spinChain[idxConv(nx+p,nnn,jjj)]*spinChain[idxConv(nx+p,nnn,jjj+1)];
		}
	}

	return abs(chainEnergy)/(N);
};

void SysConf::CalculateChainCorrelations(vector<double>& CorrSz, vector<double>& CorrStag)
{
	for(int jjj = 0; jjj < nx + p; ++jjj)
	{
		CorrSz[jjj] = 0;
		CorrStag[jjj] = 0;
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int jjj = 0;jjj < nx + p; ++jjj)
		{
			CorrSz[jjj] 	+= spinChain[idxConv(nx+p,nnn,0)]*spinChain[idxConv(nx+p,nnn,jjj)];
			CorrStag[jjj] 	+= pow(-1,jjj)*spinChain[idxConv(nx+p,nnn,0)]*spinChain[idxConv(nx+p,nnn,jjj)];
		}
	}

	for(int jjj = 0; jjj < nx + p; ++jjj)
	{
		CorrSz[jjj] = CorrSz[jjj]/N;
	}
};

void SysConf::CalculateMeanPartition()
{
	// Each column of the spins hexagon is associated to one or more partition parts, depending on its position. "nbOfParts" keeps track of this number.
	// !!! WARNING : The code suppose ny (leq) nx !!!
	int nbOfParts =  1;
	int nxInitCol =  0;
	int nyInitRow =  0;

	// First step : nbOfParts grows
	for(int lll = L + 1; lll < L + ny; ++lll)
	{
		nxInitCol = 0;
		nyInitRow = (ny - 1) - (lll - L - 1);

		GetHeights(nbOfParts,lll);
		SetParts(nbOfParts,nxInitCol,nyInitRow);

		++nbOfParts;
	}

	// Second step : nbOfParts stay constant
	for(int lll = L + ny; lll < L + nx + 1; ++lll)
	{
		nxInitCol = (lll - L - ny);
		nyInitRow = 0;

		GetHeights(nbOfParts,lll);
		SetParts(nbOfParts,nxInitCol,nyInitRow);
	}

	// Third step : nbOfParts decreases
	for(int lll = L + nx + 1; lll < L + nx + ny; ++lll)
	{
		--nbOfParts;

		nxInitCol = (lll - L - ny);
		nyInitRow = 0;

		GetHeights(nbOfParts,lll);
		SetParts(nbOfParts,nxInitCol,nyInitRow);
	}
};

void SysConf::GetHeights(int nbOfParts, int lll)
{
	int spinIdx = -1;
	int spinPos = -1;
	int nextIdx = -1;
	int nextSpinPos = -1;

	int distance = 0;
	int dimersFound = 0;
	for(int nnn = 0; nnn < N; ++nnn)
	{
		spinIdx 	= lll;
		nextIdx		= neighboursTable[idxConv(NbOfNeights,lll,2)];

		spinPos     = idxConv(N,spinIdx,nnn);
		nextSpinPos = idxConv(N,nextIdx,nnn);

		dimersFound = 0;
		distance    = 0;
		while(dimersFound < nbOfParts)
		{
			if(spinConf[spinPos]==spinConf[nextSpinPos])
			{
				// We've got a dimer!
				horizontalDimerPos[idxConv(N,dimersFound,nnn)] = p - distance + dimersFound;
				++dimersFound;
			}

			// Anyway, go to the next spin!
			++distance;
			spinIdx = nextIdx;
			nextIdx = neighboursTable[idxConv(NbOfNeights,spinIdx,2)];

			if(nextIdx==-1)
			{
				// !!! It should never get here! !!!
				break;
			}

			spinPos     = idxConv(N,spinIdx,nnn);
			nextSpinPos = idxConv(N,nextIdx,nnn);
		}
	}
};

void SysConf::SetParts(int nbOfParts, int nxInitCol, int nyInitRow)
{
	int nxCol = nxInitCol;
	int nyRow = nyInitRow;

	int partIdx   = -1;
	int heightIdx = -1;

	for(int kkk = 0; kkk < nbOfParts; ++kkk)
	{
		for(int nnn = 0; nnn < N; ++nnn)
		{
			heightIdx = idxConv(N,kkk,nnn);
			partIdx   = idxConv(N,idxConv(nx,nyRow,nxCol),nnn);
			equivPart[partIdx] = horizontalDimerPos[heightIdx];
		}
		++nxCol;
		++nyRow;
	}
};

void SysConf::GetMeanPart(vector<double> &output)
{
    for(int iii = 0; iii < ny; ++iii)
    {
		for(int jjj = 0; jjj < nx; ++jjj)
		{
			output[idxConv(nx,iii,jjj)] = 0;
		}
    }

    for(int iii = 0; iii < ny; ++iii)
    {
		for(int jjj = 0; jjj < nx; ++jjj)
		{
			for(int nnn = 0; nnn < N; ++nnn)
			{
				output[idxConv(nx,iii,jjj)] += (double)equivPart[idxConv(N,idxConv(nx,iii,jjj),nnn)]/N;
			}
		}
    }
};

void SysConf::SetSublatticeManual()
{
	int idx = 0;

	bool countA = true;
	bool countB = true;
	bool countC = true;

	SublatticeSize.resize(3,0);

	SublatticePositions.resize(L,0);

	SublatticeA.resize(0);
	SublatticeA.reserve(L/3);

	SublatticeB.resize(0);
	SublatticeB.reserve(L/3);

	SublatticeC.resize(0);
	SublatticeC.reserve(L/3);

	for(int lll = 0; lll < L/3; ++lll)
	{
		countA = true;
		countB = true;
		countC = true;

		idx = 3*lll;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighboursTable[idxConv(6,idx,kkk)]>L-1)
			{
				countA = false;
				break;
			}
		}

		if(countA)
		{
			SublatticeA.push_back(idx);
			++SublatticeSize[0];
			SublatticePositions[idx] = 0;
		}

		idx = 3*lll+1;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighboursTable[idxConv(6,idx,kkk)]>L-1)
			{
				countB = false;
				break;
			}
		}

		if(countB)
		{
			SublatticeB.push_back(idx);
			++SublatticeSize[1];
			SublatticePositions[idx] = 1;
		}

		idx = 3*lll+2;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighboursTable[idxConv(6,idx,kkk)]>L-1)
			{
				countC = false;
				break;
			}
		}

		if(countC)
		{
			SublatticeC.push_back(idx);
			++SublatticeSize[2];
			SublatticePositions[idx] = 2;
		}
	}
}

void SysConf::SetSublatticeMoessner()
{
	int idx = 0;

	SublatticeSize.resize(3,L/3);

	SublatticePositions.resize(L,0);

	SublatticeA.resize(0);
	SublatticeA.reserve(L/3);

	SublatticeB.resize(0);
	SublatticeB.reserve(L/3);

	SublatticeC.resize(0);
	SublatticeC.reserve(L/3);

	for(int iii = 0; iii < ny/3; ++iii)
	{
		for(int jjj = 0; jjj < nx/3; ++jjj)
		{
			// First line
			// A
			idx = idxConv(nx,3*iii,3*jjj);
			SublatticeA.push_back(idx);
			SublatticePositions[idx] = 0;

			// B
			idx = idxConv(nx,3*iii,3*jjj+1);
			SublatticeB.push_back(idx);
			SublatticePositions[idx] = 1;

			// C
			idx = idxConv(nx,3*iii,3*jjj+2);
			SublatticeC.push_back(idx);
			SublatticePositions[idx] = 2;

			// Second line
			// C
			idx = idxConv(nx,3*iii+1,3*jjj);
			SublatticeC.push_back(idx);
			SublatticePositions[idx] = 2;

			// A
			idx = idxConv(nx,3*iii+1,3*jjj+1);
			SublatticeA.push_back(idx);
			SublatticePositions[idx] = 0;

			// B
			idx = idxConv(nx,3*iii+1,3*jjj+2);
			SublatticeB.push_back(idx);
			SublatticePositions[idx] = 1;

			// Third line
			// B
			idx = idxConv(nx,3*iii+2,3*jjj);
			SublatticeB.push_back(idx);
			SublatticePositions[idx] = 1;

			// C
			idx = idxConv(nx,3*iii+2,3*jjj+1);
			SublatticeC.push_back(idx);
			SublatticePositions[idx] = 2;

			// A
			idx = idxConv(nx,3*iii+2,3*jjj+2);
			SublatticeA.push_back(idx);
			SublatticePositions[idx] = 0;
		}
	}
}


void SysConf::GetSublatticePeriodic(vector<double>& output, vector<double> & NSiteCount)
{
	int idx = -1;
	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 0;
	}

	for(int iii = 0; iii < ny/3; ++iii)
	{
		for(int jjj = 0; jjj < nx/3; ++jjj)
		{
			// First line
			// A
			idx = idxConv(nx,3*iii,3*jjj);
			output[0] += NSiteCount[idx];

			// B
			idx = idxConv(nx,3*iii,3*jjj+1);
			output[1] += NSiteCount[idx];

			// C
			idx = idxConv(nx,3*iii,3*jjj+2);
			output[2] += NSiteCount[idx];

			// Second line
			// C
			idx = idxConv(nx,3*iii+1,3*jjj);
			output[2] += NSiteCount[idx];

			// A
			idx = idxConv(nx,3*iii+1,3*jjj+1);
			output[0] += NSiteCount[idx];

			// B
			idx = idxConv(nx,3*iii+1,3*jjj+2);
			output[1] += NSiteCount[idx];

			// Third line
			// B
			idx = idxConv(nx,3*iii+2,3*jjj);
			output[1] += NSiteCount[idx];

			// C
			idx = idxConv(nx,3*iii+2,3*jjj+1);
			output[2] += NSiteCount[idx];

			// A
			idx = idxConv(nx,3*iii+2,3*jjj+2);
			output[0] += NSiteCount[idx];
		}
	}

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 3.*output[iii]/L;
	}

// 	sort(output.begin(), output.end());
}

void SysConf::GetSublatticeOpen(vector<double>& output, vector<double> & NSiteCount)
{

	int idx;

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 0;
	}

	// Sublattice A
	for(int iii = 0; iii < SublatticeSize[0]; ++iii)
	{
		idx = SublatticeA[iii];
		output[0] += NSiteCount[idx];
	}

	// Sublattice B
	for(int iii = 0; iii < SublatticeSize[1]; ++iii)
	{
		idx = SublatticeB[iii];
		output[1] += NSiteCount[idx];
	}

	// Sublattice C
	for(int iii = 0; iii < SublatticeSize[2]; ++iii)
	{
		idx = SublatticeC[iii];
		output[2] += NSiteCount[idx];
	}

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = output[iii]/SublatticeSize[iii];
	}

	double SubAux = -1;
	vector<double>::iterator maxSub = max_element(output.begin(),output.end());
	int distMax = distance(output.begin(),maxSub);
	vector<double>::iterator minSub = min_element(output.begin(),output.end());
	int distMin = distance(output.begin(),minSub);

	if(*maxSub - 2 > 2 - *minSub)
	{
		// Then we're in the plaquette domain
		if(distMax == 1)   // ( C A B ) -> ( A B C )
		{
			SubAux = output[0];
			output[0] = output[1];	// ( A A B ) | C
			output[1] = output[2];	// ( A B B ) | C
			output[2] = SubAux;				// ( A B C ) | C
		}

		if(distMax == 2)   // ( B C A ) -> ( A B C )
		{
			SubAux = output[2];				// ( B C A ) | A
			output[2] = output[1];  // ( B C C ) | A
			output[1] = output[0];  // ( B B C ) | A
			output[0] = SubAux;				// ( A B C ) | A
		}
	}
	else
	{
		// Then we're in the star domain
		if(distMin == 0)   // ( C A B ) -> ( A B C )
		{
			SubAux = output[0];
			output[0] = output[1];	// ( A A B ) | C
			output[1] = output[2];	// ( A B B ) | C
			output[2] = SubAux;				// ( A B C ) | C
		}

		if(distMin == 1)   // ( B C A ) -> ( A B C )
		{
			SubAux = output[2];				// ( B C A ) | A
			output[2] = output[1];  // ( B C C ) | A
			output[1] = output[0];  // ( B B C ) | A
			output[0] = SubAux;				// ( A B C ) | A
		}
	}
}

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount, vector<double> & Ndimer, vector<double> & MeanSublattice, vector<double>& DimerIndex,complex<double>& complexPhase, double& realPhase)
{
	// ---> Dimensions
	//      Ncount		: 4
	//	NSiteCount	: L
	//	MeanDimer	: L * NbOfNeights

	int test_Ncount = 0;
	int N3Count = 0;

	int index = 0;

	int spinPos = 0;
	int neighPos = 0;

	int spinValue = 0;
	int neighValue = 0;

	// Reset vectors
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		NSiteCount[iii] = 0;
		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
		{
			Ndimer[idxConv(NbOfNeights,iii,jjj)] = 0;
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int lll = 0; lll < L; ++lll)
		{
			// Count the number of dimers
			test_Ncount	= 0;
			spinPos 	= idxConv(N,lll,nnn);
			spinValue 	= (int)spinConf[spinPos];

			if(spinValue!=0)
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index 		= idxConv(NbOfNeights,lll,iii);
					neighPos	= idxConv(N,neighboursTable[index],nnn);

				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					neighValue	= weightTable[index]*(int)spinConf[neighPos];
				#else
					neighValue	= (int)spinConf[neighPos];
				#endif

					if(neighValue==spinValue)
					{
						++test_Ncount;
						Ndimer[index] += 1./N;
					}
				}

				DimerDensityProfile[spinPos] = test_Ncount;
				Ncount[test_Ncount] 	+= 1./N;
				NSiteCount[lll]		+= (double)test_Ncount/N;

				if(test_Ncount==3)
				{
					N3_Map[N3Count] = spinPos;
					++N3Count;
				}
			}

		}
	}



//	if(SimType.compare("Moessner")==0)
//	{
//		if(BorderType==1)
//		{
//			GetSublatticePeriodic(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//		else
//		{
//			GetSublatticeOpen(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//	}
//
//	if(SimType.compare("Manual")==0)
//	{
//		GetSublatticeOpen(MeanSublattice,NSiteCount);
//		std::sort(MeanSublattice.begin(),MeanSublattice.end());
//	}

//	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
//	{
//		GetSublatticeOpen(MeanSublattice,NSiteCount);
//		std::sort(MeanSublattice.begin(),MeanSublattice.end());
//	}

	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
	{
		GetSublatticeOpen(MeanSublattice,NSiteCount);
	}

//	if(SimType.compare("Manual")==0)
//	{
//		GetSublatticePeriodic(MeanSublattice,NSiteCount);
//		std::sort(MeanSublattice.begin(),MeanSublattice.end());
//	}

	// Calculate the plaquette and columnar indexes
//	PlaquetteIndex = 0;
//	Star3Index = 0;

	int Sub0Value = 0;
	int Sub1Value = 0;
	int spinNeigh1 = 0;
	int spinNeigh0 = 0;
	int index0 = 0;
	int index1 = 0;


	complexPhase = 0.;
	int chosenSublattice = -1;

	for(uint iii = 0; iii < DimerIndex.size(); ++iii)
		DimerIndex[iii] = 0;

	int lll = 0;
	int nnn = 0;
	for(int iii = 0; iii < N3Count; ++iii)
	{
		spinPos = N3_Map[iii];
		lll = spinPos/N;
		nnn = spinPos%N;

		// Complex phase
		chosenSublattice = SublatticePositions[lll];
		complexPhase += Jindex[chosenSublattice];


		// Indexes
		for(int jjj = 0; jjj < NbOfNeights/2; ++jjj)
		{
			index0 		= idxConv(NbOfNeights,lll,2*jjj);
			index1 		= idxConv(NbOfNeights,lll,(2*jjj+1));

			spinNeigh0  = idxConv(N,neighboursTable[index0],nnn);
			spinNeigh1  = idxConv(N,neighboursTable[index1],nnn);

			Sub0Value   = DimerDensityProfile[spinNeigh0];
			Sub1Value   = DimerDensityProfile[spinNeigh1];

			// DimerIndex[0] = 3,2,0
			// DimerIndex[1] = 3,1,1
			// DimerIndex[2] = 3,3,0
			// DimerIndex[3] = 3,2,1
			// DimerIndex[4] = 3,3,1

			DimerIndex[0]  += (Sub0Value==2)*(Sub1Value==0);
			DimerIndex[1]  += (Sub0Value==0)*(Sub1Value==2);

			DimerIndex[2]  += (Sub0Value==1)*(Sub1Value==1);

			DimerIndex[3]  += (Sub0Value==3)*(Sub1Value==0);
			DimerIndex[4]  += (Sub0Value==0)*(Sub1Value==3);

			DimerIndex[5]  += (Sub0Value==2)*(Sub1Value==1);
			DimerIndex[6]  += (Sub0Value==1)*(Sub1Value==2);

			DimerIndex[7]  += (Sub0Value==3)*(Sub1Value==1);
			DimerIndex[8]  += (Sub0Value==1)*(Sub1Value==3);
		}
	}

	complexPhase = complexPhase/((double)N*L);
	realPhase = cos(3*arg(complexPhase));

	for(uint iii = 0; iii < DimerIndex.size(); ++iii)
		DimerIndex[iii] = 2*DimerIndex[iii]/(N3Count*NbOfNeights);
//
//	PlaquetteIndex = PlaquetteIndex/N3Count;
//	Star3Index = Star3Index/N3Count;
}

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount, vector<double> & Ndimer, vector<double> & MeanSublattice, vector<double>& DimerIndex)
{
	// ---> Dimensions
	//      Ncount		: 4
	//	NSiteCount	: L
	//	MeanDimer	: L * NbOfNeights

	int test_Ncount = 0;
	int N3Count = 0;

	int index = 0;

	int spinPos = 0;
	int neighPos = 0;

	int spinValue = 0;
	int neighValue = 0;

	// Reset vectors
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		NSiteCount[iii] = 0;
		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
		{
			Ndimer[idxConv(NbOfNeights,iii,jjj)] = 0;
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int lll = 0; lll < L; ++lll)
		{
			// Count the number of dimers
			test_Ncount	= 0;
			spinPos 	= idxConv(N,lll,nnn);
			spinValue 	= (int)spinConf[spinPos];

			if(spinValue!=0)
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index 		= idxConv(NbOfNeights,lll,iii);
					neighPos	= idxConv(N,neighboursTable[index],nnn);

				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					neighValue	= weightTable[index]*(int)spinConf[neighPos];
				#else
					neighValue	= (int)spinConf[neighPos];
				#endif

					if(neighValue==spinValue)
					{
						++test_Ncount;
						Ndimer[index] += 1./N;
					}
				}

				DimerDensityProfile[spinPos] = test_Ncount;
				Ncount[test_Ncount] 	+= 1./N;
				NSiteCount[lll]		+= (double)test_Ncount/N;

				if(test_Ncount==3)
				{
					N3_Map[N3Count] = spinPos;
					++N3Count;
				}
			}

		}
	}

//	if(SimType.compare("Moessner")==0)
//	{
//		if(BorderType==1)
//		{
//			GetSublatticePeriodic(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//		else
//		{
//			GetSublatticeOpen(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//	}
//
//	if(SimType.compare("Manual")==0)
//	{
//		GetSublatticeOpen(MeanSublattice,NSiteCount);
//		std::sort(MeanSublattice.begin(),MeanSublattice.end());
//	}

	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
	{
		GetSublatticeOpen(MeanSublattice,NSiteCount);
	}

	// Calculate the plaquette and columnar indexes
//	PlaquetteIndex = 0;
//	Star3Index = 0;

	int Sub0Value = 0;
	int Sub1Value = 0;
	int spinNeigh1 = 0;
	int spinNeigh0 = 0;
	int index0 = 0;
	int index1 = 0;

	for(uint iii = 0; iii < DimerIndex.size(); ++iii)
		DimerIndex[iii] = 0;

	int lll = 0;
	int nnn = 0;
	for(int iii = 0; iii < N3Count; ++iii)
	{
		spinPos = N3_Map[iii];
		lll = spinPos/N;
		nnn = spinPos%N;

		for(int jjj = 0; jjj < NbOfNeights/2; ++jjj)
		{
			index0 		= idxConv(NbOfNeights,lll,2*jjj);
			index1 		= idxConv(NbOfNeights,lll,(2*jjj+1));

			spinNeigh0  = idxConv(N,neighboursTable[index0],nnn);
			spinNeigh1  = idxConv(N,neighboursTable[index1],nnn);

			Sub0Value   = DimerDensityProfile[spinNeigh0];
			Sub1Value   = DimerDensityProfile[spinNeigh1];

			// DimerIndex[0] = 3,2,0
			// DimerIndex[1] = 3,1,1
			// DimerIndex[2] = 3,3,0
			// DimerIndex[3] = 3,2,1
			// DimerIndex[4] = 3,3,1

			DimerIndex[0]  += (Sub0Value==2)*(Sub1Value==0);
			DimerIndex[1]  += (Sub0Value==0)*(Sub1Value==2);

			DimerIndex[2]  += (Sub0Value==1)*(Sub1Value==1);

			DimerIndex[3]  += (Sub0Value==3)*(Sub1Value==0);
			DimerIndex[4]  += (Sub0Value==0)*(Sub1Value==3);

			DimerIndex[5]  += (Sub0Value==2)*(Sub1Value==1);
			DimerIndex[6]  += (Sub0Value==1)*(Sub1Value==2);

			DimerIndex[7]  += (Sub0Value==3)*(Sub1Value==1);
			DimerIndex[8]  += (Sub0Value==1)*(Sub1Value==3);
		}
	}

	for(uint iii = 0; iii < DimerIndex.size(); ++iii)
		DimerIndex[iii] = 2*DimerIndex[iii]/(N3Count*NbOfNeights);
//
//	PlaquetteIndex = PlaquetteIndex/N3Count;
//	Star3Index = Star3Index/N3Count;
}

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount, vector<double> & Ndimer, vector<double> & MeanSublattice)
{
	// ---> Dimensions
	//      Ncount		: 4
	//	NSiteCount	: L
	//	MeanDimer	: L * NbOfNeights

	int test_Ncount = 0;

	int index = 0;

	int spinPos = 0;
	int neighPos = 0;

	int spinValue = 0;
	int neighValue = 0;

	// Reset vectors
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		NSiteCount[iii] = 0;
		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
		{
			Ndimer[idxConv(NbOfNeights,iii,jjj)] = 0;
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int lll = 0; lll < L; ++lll)
		{
			// Count the number of dimers
			test_Ncount	= 0;
			spinPos 	= idxConv(N,lll,nnn);
			spinValue 	= (int)spinConf[spinPos];

			if(spinValue!=0)
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index 		= idxConv(NbOfNeights,lll,iii);
					neighPos	= idxConv(N,neighboursTable[index],nnn);

				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					neighValue	= weightTable[index]*(int)spinConf[neighPos];
				#else
					neighValue	= (int)spinConf[neighPos];
				#endif

					if(neighValue==spinValue)
					{
						++test_Ncount;
						Ndimer[index] += 1./N;
					}
				}

				DimerDensityProfile[spinPos] = test_Ncount;
				Ncount[test_Ncount] 	+= 1./N;
				NSiteCount[lll]		+= (double)test_Ncount/N;
			}

		}
	}

//	if(SimType.compare("Moessner")==0)
//	{
//		if(BorderType==1)
//		{
//			GetSublatticePeriodic(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//		else
//		{
//			GetSublatticeOpen(MeanSublattice,NSiteCount);
//			std::sort(MeanSublattice.begin(),MeanSublattice.end());
//		}
//	}
//
//	if(SimType.compare("Manual")==0)
//	{
//		GetSublatticeOpen(MeanSublattice,NSiteCount);
//		std::sort(MeanSublattice.begin(),MeanSublattice.end());
//	}

	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
	{
		GetSublatticeOpen(MeanSublattice,NSiteCount);
	}
}

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount, vector<double> & Ndimer)
{
	// ---> Dimensions
	//      Ncount		: 4
	//	NSiteCount	: L
	//	MeanDimer	: L * NbOfNeights

	int test_Ncount = 0;

	int index = 0;

	int spinPos = 0;
	int neighPos = 0;

	int spinValue = 0;
	int neighValue = 0;

	// Reset vectors
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		NSiteCount[iii] = 0;
		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
		{
			Ndimer[idxConv(NbOfNeights,iii,jjj)] = 0;
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		for(int lll = 0; lll < L; ++lll)
		{
			// Count the number of dimers
			test_Ncount	= 0;
			spinPos 	= idxConv(N,lll,nnn);
			spinValue 	= (int)spinConf[spinPos];

			if(spinValue!=0)
			{
				for(int iii = 0; iii < NbOfNeights;++iii)
				{
					index 		= idxConv(NbOfNeights,lll,iii);
					neighPos	= idxConv(N,neighboursTable[index],nnn);

				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					neighValue	= weightTable[index]*(int)spinConf[neighPos];
				#else
					neighValue	= (int)spinConf[neighPos];
				#endif

					if(neighValue==spinValue)
					{
						++test_Ncount;
						Ndimer[index] += 1./N;
					}
				}

				Ncount[test_Ncount] 	+= 1./N;
				NSiteCount[lll]		+= (double)test_Ncount/N;
			}

		}
	}
}

void SysConf::GetN3Correlation(double& corr, vector<double>& Star3Network)
{
	int 		spinPos 	= 0;
	int 		neighPos 	= 0;

	int lll = 0;
	int nnn = 0;

	int index = 0;

	corr = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		Star3Network[iii] = 0;
	}

	for(int fff = 0; fff < TotalOldNbF; ++fff)
	{
		spinPos = flippableSpinLists[fff];
		nnn = spinPos % N;
		lll = spinPos / N;

		for(int iii = 0; iii < NbOfNeights; ++iii)
		{
			index = idxConv(NbOfNeights,lll,iii);
			neighPos = idxConv(N,neighboursTable[index],nnn);

			if(flippableSpinPositions[neighPos]!=-1)
			{
				corr += 1./N;
				Star3Network[lll] += 1./N;

				break;
			}
		}
	}
}


void	SysConf::GetN3N3SpacialCorrelation(vector<double>& N3SpatialCorrelation, vector<double>& N3LocalDensities)
{
	int 	spinPos 		= 0;
	int 	distance	 	= 0;
	int		hashPosition	= 0;

	int lll = 0;
	int jjj = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		N3LocalDensities[iii] = 0;
	}

	for(int iii = 0; iii < L*nbOfDists; ++iii)
	{
		N3SpatialCorrelation[iii] = 0;
	}

	// For all flippable spins
	for(int fff = 0; fff < TotalOldNbF; ++fff)
	{
		spinPos = flippableSpinLists[fff];
		lll = spinPos / N;

		N3LocalDensities[lll] += 1./N;
		for(int ggg = 0; ggg < TotalOldNbF; ++ggg)
		{
			jjj = flippableSpinLists[ggg] / N;

			distance = SqrDistancesTable[idxConv(L,lll,jjj)];
			hashPosition = SqrDistancesPos[distance];
			N3SpatialCorrelation[idxConv(nbOfDists,lll,hashPosition)] += 1./N;
		}
	}
}

//
//void	SysConf::GetSzSzCorrelation(vector<double>& corr)
//{
//	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] = 0;
//	}
//
//	double firstSpin = 0;
//	for(int mmm = 0; mmm < N; ++mmm)
//	{
//		for(int iii = 0; iii < L; ++iii)
//		{
//			firstSpin = spinConf[idxConv(N,iii,mmm)];
//			for(int nnn = 0; nnn < mmm; ++nnn)
//			{
//				corr[N-mmm+nnn] += firstSpin*spinConf[idxConv(N,iii,nnn)];
//			}
//
//			for(int nnn = mmm; nnn < N; ++nnn)
//			{
//				corr[nnn-mmm] 	+= firstSpin*spinConf[idxConv(N,iii,nnn)];
//			}
//		}
//	}
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] 		= corr[nnn]/(L*N);
//	}
//}
//
//void	SysConf::GetDimerDimerCorrelation(vector<double>& corr)
//{
//	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
//	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] = 0;
//	}
//
//	double 	firstSpin = 0;
//	double 	dimer = 0;
//	double 	firstDimer = 0;
//	int 	posNeigh = 0;
//	int 	index = 0;
//
//	for(int mmm = 0; mmm < N; ++mmm)
//	{
//		for(int iii = 0; iii < L; ++iii)
//		{
//			firstSpin = spinConf[idxConv(N,iii,mmm)];
//
//			for(int jjj = 0; jjj < NbOfNeights; ++jjj)
//			{
//				index = idxConv(NbOfNeights,iii,jjj);
//				posNeigh  = neighboursTable[index];
//
//
//			#if BORDER_TYPE == 4
//			// ANTI-PERIODIC BORDERS
//				firstDimer = (firstSpin*weightTable[index]*spinConf[idxConv(N,posNeigh,mmm)] + 1)/2.;
//			#else
//				firstDimer = (firstSpin*spinConf[idxConv(N,posNeigh,mmm)] + 1)/2.;
//			#endif
//
//
//				for(int nnn = 0; nnn < mmm; ++nnn)
//				{
//				#if BORDER_TYPE == 4
//				// ANTI-PERIODIC BORDERS
//					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#else
//					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#endif
//					corr[N-mmm+nnn] += firstDimer*dimer;
//				}
//
//				for(int nnn = mmm; nnn < N; ++nnn)
//				{
//				#if BORDER_TYPE == 4
//				// ANTI-PERIODIC BORDERS
//					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#else
//					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#endif
//					corr[nnn-mmm] 	+= firstDimer*dimer;
//				}
//			}
//		}
//	}
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] 		= corr[nnn]/(NbOfNeights*L*N);
//	}
//}

void	SysConf::GetSzSzCorrelation(vector<double>& corr)
{
	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] = 0;
	}

	double firstSpin = 0;

	for(int iii = 0; iii < L; ++iii)
	{
		firstSpin = spinConf[idxConv(N,iii,N/2)];
		for(int nnn = 0; nnn < N/2; ++nnn)
		{
			corr[N/2+nnn] += firstSpin*spinConf[idxConv(N,iii,nnn)];
		}

		for(int nnn = N/2; nnn < N; ++nnn)
		{
			corr[nnn-N/2] 	+= firstSpin*spinConf[idxConv(N,iii,nnn)];
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] 		= corr[nnn]/L;
	}
}

void	SysConf::GetDimerDimerCorrelation(vector<double>& corr)
{
	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] = 0;
	}

	double 	firstSpin = 0;
	double 	dimer = 0;
	double 	firstDimer = 0;
	int 	posNeigh = 0;
	int 	index = 0;

		for(int iii = 0; iii < L; ++iii)
		{
			firstSpin = spinConf[idxConv(N,iii,N/2)];

			for(int jjj = 0; jjj < NbOfNeights; ++jjj)
			{
				index = idxConv(NbOfNeights,iii,jjj);
				posNeigh  = neighboursTable[index];


			#if BORDER_TYPE == 4
			// ANTI-PERIODIC BORDERS
				firstDimer = (firstSpin*weightTable[index]*spinConf[idxConv(N,posNeigh,N/2)] + 1)/2.;
			#else
				firstDimer = (firstSpin*spinConf[idxConv(N,posNeigh,N/2)] + 1)/2.;
			#endif


				for(int nnn = 0; nnn < N/2; ++nnn)
				{
				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#else
					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#endif
					corr[N/2+nnn] += firstDimer*dimer;
				}

				for(int nnn = N/2; nnn < N; ++nnn)
				{
				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#else
					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#endif
					corr[nnn-N/2] 	+= firstDimer*dimer;
				}
			}
		}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] 		= corr[nnn]/(NbOfNeights*L);
	}
}

double  SysConf::GetQEnergy_N3()
{
	double energyOut = GetKinEnergy() + GetPotEnergy_N3();

	return energyOut;
};

double  SysConf::GetQEnergy_N0()
{
	double energyOut = GetKinEnergy() + GetPotEnergy_N0();

	return energyOut;
};

double  SysConf::GetQEnergy_Mixed()
{
	double energyOut = GetKinEnergy() + GetPotEnergy_Mixed();

	return energyOut;
};

double  SysConf::GetQNewEnergy_N3()
{
        double energyOut = GetNewKinEnergy() + GetPotEnergy_N3();

        return energyOut;
};
 
double  SysConf::GetQNewEnergy_N0()
{
        double energyOut = GetNewKinEnergy() + GetPotEnergy_N0();

        return energyOut;
};

double  SysConf::GetQNewEnergy_Mixed()
{
        double energyOut = GetNewKinEnergy() + GetPotEnergy_Mixed(); 
	return energyOut;
};

void  SysConf::GetAllQEnergy_Mixed(vector<double> & dummyVec)
{
	dummyVec[0] = GetKinEnergy();
	dummyVec[1] = GetPotEnergy_Mixed();
};

void  SysConf::GetAllNewQEnergy_Mixed(vector<double> & dummyVec)
{
	dummyVec[0] = GetNewKinEnergy();
	dummyVec[1] = GetPotEnergy_Mixed();
};

double  SysConf::GetPotEnergy_N3()
{
	double energyOut = ((alpha3*Vt + shift3)*TotalOldNbF)/N;

	return energyOut;
};

double  SysConf::GetPotEnergy_N0()
{
	double energyOut = ((alpha0*Vt + shift0)*GetN0())/N;

	return energyOut;
};

double  SysConf::GetPotEnergy_Mixed()
{
	double energyOut = GetPotEnergy_N0() + GetPotEnergy_N3();

	return energyOut;
};

double  SysConf::GetKinEnergy()
{
	double energyOut = 0;
	double cte = 1./(N*sinh(2.*Kz));  	// !!! Kz = Delta Beta

	if(N>1)
	{
        for(int nnn = 0; nnn < N - 1; ++nnn)
        {
            for(int iii = 0; iii < L; ++iii)
            {
               	energyOut += cte*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn+1)];
            }
        }

		for(int iii = 0; iii < L; ++iii)
		{
			energyOut += cte*spinConf[idxConv(N,iii,N-1)]*spinConf[idxConv(N,iii,0)];
		}
	}
	else
	{
		for(int iii = 0; iii < L; ++iii)
		{
			energyOut += cte*spinConf[iii]*spinConf[iii];
		}
	}

//	double cte = 1./(sinh(2.*Kz));  	// !!! Kz = Delta Beta
//	for(int iii = 0; iii < L; ++iii)
//	{
//		energyOut += cte*spinConf[idxConv(N,iii,0)]*spinConf[idxConv(N,iii,1)];
//	}

	energyOut -= L/tanh(2.*Kz);
	return energyOut;
};

double  SysConf::GetNewKinEnergy()
{
	double energyOut = 0;
	double cte = 1./(N*Kz);        // !!! Kz = Delta Beta

	if(N>1)
	{
		for(int nnn = 0; nnn < N - 1; ++nnn)
		{
			for(int iii = 0; iii < L; ++iii)
			{
				if(spinConf[idxConv(N,iii,nnn)]==-spinConf[idxConv(N,iii,nnn+1)])
				{
					energyOut -= cte;
				}
			}
		}

		for(int iii = 0; iii < L; ++iii)
		{
			if(spinConf[idxConv(N,iii,N-1)]==-spinConf[idxConv(N,iii,0)])
			{
				energyOut -= cte;
			}
		}
	}
	else
	{
		energyOut = 0;
	}

	return energyOut;
};

void  SysConf::GetKinEnergyDensity(vector<double> & KinOut)
{
	double cte = 1./(N*sinh(2.*Kz));  	// !!! Kz = Delta Beta

	if(N>1)
	{
		for(int iii = 0; iii < L; ++iii)
		{
			KinOut[iii] = 0;
			for(int nnn = 0; nnn < N - 1; ++nnn)
			{
				KinOut[iii] += cte*spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,iii,nnn+1)];
			}
		}

		for(int iii = 0; iii < L; ++iii)
		{
			KinOut[iii] += cte*spinConf[idxConv(N,iii,N-1)]*spinConf[idxConv(N,iii,0)];
		}
	}
	else
	{
		for(int iii = 0; iii < L; ++iii)
		{
			KinOut[iii] += cte*spinConf[iii]*spinConf[iii];
		}
	}

//	double cte = 1./(sinh(2.*Kz));  	// !!! Kz = Delta Beta
//	for(int iii = 0; iii < L; ++iii)
//	{
//		energyOut += cte*spinConf[idxConv(N,iii,0)]*spinConf[idxConv(N,iii,1)];
//	}

	for(int iii = 0; iii < L; ++iii)
	{
		KinOut[iii] -= 1./tanh(2.*Kz);
	}
};

// *******************************************************************
// Conf randomizers
// *******************************************************************
void SysConf::RandomSpinFlips(int flips)
{
	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	toRemove.resize(N,vector<int>(NbOfNeights/2,-1));
	toAdd.resize(N,vector<int>(NbOfNeights/2,-1));

	int nbToExchange;
	int posToChange;
	int spinToMove;
	int LastPosition = TotalOldNbF - 1;

	for(int iii = 0; iii < flips; ++iii)
	{
		absIndex 		= GetAbsIndex(TotalOldNbF/*,m_rng*/);
		spinToFlip 	= flippableSpinLists[absIndex];
		nnnCluster	= spinToFlip % N;
		kkkCluster	= spinToFlip / N;

		TestUpdate(nnnCluster);
		TotalNewNbF += nbToAdd[nnnCluster]-nbToRemove[nnnCluster];

		nbToExchange = min(nbToRemove[nnnCluster],nbToAdd[nnnCluster]);
		posToChange = 0;
		spinToMove = 0;

		// > Spin conf update
		spinConf[idxConv(N,kkkCluster,nnnCluster)] = -spinConf[idxConv(N,kkkCluster,nnnCluster)];

		// > Spin tables update
		// > We have 3 situations regarding 'flippableSpins'
		//   1) nbToAdd == nbToRemove : no holes, no additional terms
		//   2) nbToAdd < nbToRemove  : we'll have some holes
		//   3) nbToAdd > nbToRemove  : no holes, but we'll have to add some terms
		//
		// > We can exchange the spins over the first 'Min(nbToRemove,nbToAdd)' spins
		//   to add/remove. Then :
		//   *  If we have the case c1), there's nothing more to do
		//   *  Case 2), we have to compact the lists
		//   *  Case 3), we have to add elements

		// > Exchange positions
		// * If we're only adding or removing spins, we won't enter this loop
		//   because nbToExchange == 0

		for(int iii = 0; iii < nbToExchange; ++iii)
		{
			// > Get position to exchange
			posToChange = flippableSpinPositions[toRemove[nnnCluster][iii]];

			// > Exchange the spins
			flippableSpinLists[posToChange] = toAdd[nnnCluster][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toRemove[nnnCluster][iii]] = -1;
			flippableSpinPositions[toAdd[nnnCluster][iii]] = posToChange;
		}

		// > Add elements
		// * If we're on the case 2), we won't enter this loop
		//   because nbToExchange == nbToAdd
		for(int iii = nbToExchange; iii < nbToAdd[nnnCluster]; ++iii)
		{
			// > Update the LastPosition
			++LastPosition;

			// > Add the spin
			flippableSpinLists[LastPosition] = toAdd[nnnCluster][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toAdd[nnnCluster][iii]] = LastPosition;
		}

		// > Remove elements and fill holes
		// * If we're on the case 3), we won't enter this loop
		//   because nbToExchange == nbToRemove
		for(int iii = nbToExchange; iii < nbToRemove[nnnCluster]; ++iii)
		{
			// > Get position to remove
			posToChange = flippableSpinPositions[toRemove[nnnCluster][iii]];

			// > Remove spin
			flippableSpinLists[posToChange] = -1;
			flippableSpinPositions[toRemove[nnnCluster][iii]] = -1;

			if(posToChange < LastPosition)
			{
				// Then we created a hole -> we have to update!
				spinToMove = flippableSpinLists[LastPosition];

				// Move the element
				flippableSpinLists[posToChange] = spinToMove;
				flippableSpinLists[LastPosition] = -1;

				// Update position
				flippableSpinPositions[spinToMove] = posToChange;
			}

			// > Update the LastPosition
			--LastPosition;
		}

		TotalOldNbF = TotalNewNbF;
	}
}

void SysConf::SetRandomColFlips()
{
	nbToRemove.resize(N,0);
	nbToAdd.resize(N,0);

	toRemove.resize(N,vector<int>(NbOfNeights/2,-1));
	toAdd.resize(N,vector<int>(NbOfNeights/2,-1));

	dummyToAdd.resize(NbOfNeights/2,-1);
	dummyToRemove.resize(NbOfNeights/2,-1);
}
void SysConf::RandomColFlip()
{
	int nbToExchange;
	int posToChange;
	int spinToMove;
	int LastPosition = TotalOldNbF - 1;

	int layerNeigh = -1;

	absIndex 		= GetAbsIndex(TotalOldNbF/*,m_rng*/);
	spinToFlip 	= flippableSpinLists[absIndex];
	nnnCluster	= spinToFlip % N;
	kkkCluster	= spinToFlip / N; // Important !!!

	TestUpdate(nnnCluster);

	for(int jjj = 0; jjj < nbToRemove[nnnCluster]; ++jjj)
	{
		dummyToRemove[jjj] = toRemove[nnnCluster][jjj]/N;
	}
	for(int jjj = 0; jjj < nbToAdd[nnnCluster]; ++jjj)
	{
		dummyToAdd[jjj] = toAdd[nnnCluster][jjj]/N;
	}

	for(int nnn = 0; nnn < nnnCluster; ++nnn)
	{
		nbToRemove[nnn]=nbToRemove[nnnCluster];
		nbToAdd[nnn]=nbToAdd[nnnCluster];

		for(int jjj = 0; jjj < nbToRemove[nnn]; ++jjj)
		{
			layerNeigh = idxConv(N,dummyToRemove[jjj],nnn);
			toRemove[nnn][jjj] = layerNeigh;
		}

		for(int jjj = 0; jjj < nbToAdd[nnn]; ++jjj)
		{
			layerNeigh = idxConv(N,dummyToAdd[jjj],nnn);
			toAdd[nnn][jjj] = layerNeigh;
		}
	}

	for(int nnn = nnnCluster+1; nnn < N; ++nnn)
	{
		nbToRemove[nnn]=nbToRemove[nnnCluster];
		nbToAdd[nnn]=nbToAdd[nnnCluster];

		for(int jjj = 0; jjj < nbToRemove[nnn]; ++jjj)
		{
			layerNeigh = idxConv(N,dummyToRemove[jjj],nnn);
			toRemove[nnn][jjj] = layerNeigh;
		}

		for(int jjj = 0; jjj < nbToAdd[nnn]; ++jjj)
		{
			layerNeigh = idxConv(N,dummyToAdd[jjj],nnn);
			toAdd[nnn][jjj] = layerNeigh;
		}
	}

	TotalNewNbF += N*(nbToAdd[nnnCluster]-nbToRemove[nnnCluster]);

	nbToExchange = min(nbToRemove[nnnCluster],nbToAdd[nnnCluster]);
	posToChange = 0;
	spinToMove = 0;

	// > Spin conf update
	for(int nnn = 0; nnn < N; ++nnn)
	{
		spinConf[idxConv(N,kkkCluster,nnn)] = -spinConf[idxConv(N,kkkCluster,nnn)];

		// > Spin tables update
		// > We have 3 situations regarding 'flippableSpins'
		//   1) nbToAdd == nbToRemove : no holes, no additional terms
		//   2) nbToAdd < nbToRemove  : we'll have some holes
		//   3) nbToAdd > nbToRemove  : no holes, but we'll have to add some terms
		//
		// > We can exchange the spins over the first 'Min(nbToRemove,nbToAdd)' spins
		//   to add/remove. Then :
		//   *  If we have the case c1), there's nothing more to do
		//   *  Case 2), we have to compact the lists
		//   *  Case 3), we have to add elements

		// > Exchange positions
		// * If we're only adding or removing spins, we won't enter this loop
		//   because nbToExchange == 0


		for(int iii = 0; iii < nbToExchange; ++iii)
		{
			// > Get position to exchange
			posToChange = flippableSpinPositions[toRemove[nnn][iii]];

			// > Exchange the spins
			flippableSpinLists[posToChange] = toAdd[nnn][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toRemove[nnn][iii]] = -1;
			flippableSpinPositions[toAdd[nnn][iii]] = posToChange;
		}

		// > Add elements
		// * If we're on the case 2), we won't enter this loop
		//   because nbToExchange == nbToAdd
		for(int iii = nbToExchange; iii < nbToAdd[nnn]; ++iii)
		{
			// > Update the LastPosition
			++LastPosition;

			// > Add the spin
			flippableSpinLists[LastPosition] = toAdd[nnn][iii];

			// > Update 'flippableSpinPositions'
			flippableSpinPositions[toAdd[nnn][iii]] = LastPosition;
		}

		// > Remove elements and fill holes
		// * If we're on the case 3), we won't enter this loop
		//   because nbToExchange == nbToRemove
		for(int iii = nbToExchange; iii < nbToRemove[nnn]; ++iii)
		{
			// > Get position to remove
			posToChange = flippableSpinPositions[toRemove[nnn][iii]];

			// > Remove spin
			flippableSpinLists[posToChange] = -1;
			flippableSpinPositions[toRemove[nnn][iii]] = -1;

			if(posToChange < LastPosition)
			{
				// Then we created a hole -> we have to update!
				spinToMove = flippableSpinLists[LastPosition];

				// Move the element
				flippableSpinLists[posToChange] = spinToMove;
				flippableSpinLists[LastPosition] = -1;

				// Update position
				flippableSpinPositions[spinToMove] = posToChange;
			}

			// > Update the LastPosition
			--LastPosition;
		}
	}
	TotalOldNbF = TotalNewNbF;
}

// *******************************************************************
// Output/debug
// *******************************************************************
void SysConf::PrintConfigStatus()
{
	string DummyInitCond = "Unrecognized";
	string DummyBorderCond = "Unrecognized";

	if(BorderType==0||BorderType==2)
	{
		DummyBorderCond = "Free";
	}
	else if(BorderType == 1)
	{
		if(SimType.compare("Part")==0)
		{
			DummyBorderCond = "Closed";
		}
		else if(SimType.compare("Moessner")==0)
		{
			DummyBorderCond = "Periodic";
		}
	}

	if(initCondType == 0)
	{
		DummyInitCond = "Arctic";
	}
	else if(initCondType == 1)
	{
		DummyInitCond = "Random";
	}
	else if(initCondType == 2)
	{
		DummyInitCond = "Star 3";
	}

	cout << " Geometry " << endl;
	cout << "     Type                 : " << SimType << endl;
	cout << "     Boundary             : " << DummyBorderCond << endl;
	cout << "     Initial cond         : " << DummyInitCond << endl;
	cout << "     Sizes (nx X ny X p)  : " << nx << " " << ny << " " << p << endl;
	cout << "     Bulk sites           : " << L << endl;
	cout << "     Total sites          : " << LTot << endl << endl;

	cout << " Monte Carlo parameters" << endl;
	cout << "     Stacks               : " << N << endl;
	cout << "     Kz                   : " << Kz << endl;
	cout << "     Kt                   : " << Kt << endl;
	cout << "     Vt                   : " << Vt << endl << endl;

	cout << " Order parameters" << endl;
	CalculateE();
	cout << "     Energy               : " << energy << endl;
	double PrintedMag = GetMagnetization();
	cout << "     Magnetization        : " << PrintedMag << endl;
	vector<double> PrintedNf(4,0);
	GetNf(PrintedNf);
	cout << "     N_f (f = 0, 1, 2, 3) : " << PrintedNf[0] << " " << PrintedNf[1] << " " << PrintedNf[2] << " " << PrintedNf[3] << endl << endl;


//
// 	int		N;
// 	int 		NDimers;
// 	int		NbOfNeights;
//
// 	int 		TotalOldNbF;		// Number of flippable sites of the actual configuration
// 	int 		TotalNewNbF;		// Number of flippable sites of the proposed configuration
//
// 	// ---> MC parameters
// 	double		Kz;
// 	double		Kt;
// 	double		q;
// 	double		Vt;
//
// 	// ---> Simulation type
// 	string		SimType;
// 	int 		initCondType;
};

void SysConf::PrintDebug(ostream& out)
{
	out << " ### DEBUG Print SysConf ###" << endl;
	out << " Type                   = " << SimType << endl;
	out << " L                      = " << L << endl;
	out << " LTot                   = " << LTot << endl;
	out << " LBorder                = " << LBorder << endl;
	out << " NDimers                = " << NDimers << endl;
	out << " NbOfNeights            = " << NbOfNeights << endl;
	out << " N                      = " << N << endl << endl;

	out << " TotalOldNbF            = " << TotalOldNbF << endl;
	out << " TotalNewNbF            = " << TotalNewNbF << endl;

	out << " neighboursTable        = " << endl;

	for(int iii = 0; iii < LTot; ++iii)
	{
		for(int jjj = 0; jjj < NbOfNeights; ++jjj)
		{
			out << " " << neighboursTable[idxConv(NbOfNeights,iii,jjj)];
		}
		out << endl;
	}

	out << " spinConf               = " << endl;
	for(int jjj = 0; jjj < N; ++jjj)
	{
		for(int iii = 0; iii < L; ++iii)
		{
			out << " " << (int)spinConf[idxConv(N,iii,jjj)];
		}
		out << endl;
	}
	out << endl;

	out << " flippableSpinLists     = " << endl;
	for(uint iii = 0; iii < flippableSpinLists.size(); ++iii)
	{
		out << " " << flippableSpinLists[iii];
	}
	out << endl;

	out << " flippableSpinPositions = " << endl;
	for(int iii = 0; iii < LTot; ++iii)
	{
		for(int jjj = 0; jjj < N; ++jjj)
		{
			out << " " << flippableSpinPositions[idxConv(N,iii,jjj)];
		}
		out << endl;
	}
	out << endl;
};

void SysConf::SetPositions()
{
	m_x.resize(L,0);
	m_y.resize(L,0);

	if(SimType.compare("Part")==0)
	{
		int counter = 0;
		vector<double> init_x(p + nx - 1,-1);
		vector<double> init_y(p + nx - 1,-1);

		vector<int> lineStart(p + nx - 1,-1);
		vector<int> lineSize(p + nx - 1,-1);

		// > Set limits
		for(int iii = 0; iii < p; ++iii)
		{
			lineStart[iii] = 0;
			init_x[iii] = 0;
			init_y[iii] = (p - 1 - iii);
		}
		for(int iii = p; iii < p + nx -1; ++iii)
		{
			lineStart[iii] = (iii - p + 1);
			init_x[iii] = (iii - p + 1)*sqrt(3)/2;
			init_y[iii] = -(iii - p + 1)*0.5;
		}

		for(int iii = 0; iii < nx; ++iii)
		{
			lineSize[iii] = (ny + iii) - lineStart[iii];
		}
		for(int iii = nx; iii < p + nx - 1; ++iii)
		{
			lineSize[iii] = (nx - 1) + ny - lineStart[iii];
		}

		for(int iii = 0; iii < p + nx - 1; ++iii)
		{
			for(int jjj = 0; jjj < lineSize[iii]; ++jjj)
			{
				m_x[counter] = init_x[iii] + jjj*sqrt(3)/2;
				m_y[counter] = init_y[iii] + 0.5*jjj;
				++counter;
			}
		}
//		for(int iii = 0; iii < counter; ++iii)
//		{
//			cout << m_x[iii] << " " << m_y[iii] << endl;
//		}
	}
	else if(SimType.compare("Moessner")==0)
	{
		int counter = 0;
		double init_x = 0;
		double init_y = 0;

		for(int jjj = 0; jjj < ny; ++jjj)
		{
			init_x = jjj*0.5;
			init_y = -jjj*sqrt(3)/2;

			for(int iii = 0; iii < nx; ++iii)
			{
				m_x[counter] = init_x + iii;
				m_y[counter] = init_y;
				++counter;
			}
		}
	}
}

int SysConf::MinimalSqrDistance(int pointA, int pointB)
{
	double h_shift_x = nx;
	double h_shift_y = 0;

	double v_shift_x = -ny/2.;
	double v_shift_y = sqrt(3)*ny/2.;

	double dx = m_x[pointA]-m_x[pointB];
	double dy = m_y[pointA]-m_y[pointB];

	double minimalDist = pow(dx,2)+pow(dy,2);

	int wh = 0;
	int wv = 0;

	if(BorderType==1)
	{
		// Periodic - must check borders
		double tempDist = -1;

		// (-,+)
		wh = -1;
		wv = 1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (-,0)
		wh = -1;
		wv = 0;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (-,-)
		wh = -1;
		wv = -1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (0,+)
		wh = 0;
		wv = 1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (0,-)
		wh = 0;
		wv = -1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (+,+)
		wh = 1;
		wv = 1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (+,0)
		wh = 1;
		wv = 0;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}

		// (+,-)
		wh = 1;
		wv = -1;
		tempDist = pow(dx + wh*h_shift_x + wv*v_shift_x,2)+pow(dy + wh*h_shift_y + wv*v_shift_y,2);
		if(tempDist < minimalDist)
		{
			minimalDist = tempDist;
		}
	}

	return round(minimalDist);
}

void SysConf::SetDistances()
{
	int dummy = 0;

	SqrDistancesTable.resize(L*L,0);

	set<int> temporarySqrDistancesKey;

	// Build distance table
	for(int iii = 0;iii < L; ++iii)
	{
		for(int jjj = 0; jjj <= iii; ++jjj)
		{
			dummy = MinimalSqrDistance(iii,jjj);
			SqrDistancesTable[idxConv(L,iii,jjj)] = dummy;
			SqrDistancesTable[idxConv(L,jjj,iii)] = dummy;
			temporarySqrDistancesKey.insert(dummy);
		}
	}

	// Set maximum distance
	if(BorderType==1)
	{
		maximumSqrDist = nx*nx/3;
	}
	else
	{
		maximumSqrDist = nx*nx+ny*ny+nx*ny;
	}

	// Convert the temporary key on a permanent one
	SqrDistancesKey.resize(temporarySqrDistancesKey.size());
	SqrDistancesPos.resize(maximumSqrDist+1,-1);

	nbOfDists = SqrDistancesKey.size();
	copy(temporarySqrDistancesKey.begin(),temporarySqrDistancesKey.end(),SqrDistancesKey.begin());

	for(int iii = 0; iii < nbOfDists; ++iii)
	{
		SqrDistancesPos[SqrDistancesKey[iii]] = iii;
	}
}
