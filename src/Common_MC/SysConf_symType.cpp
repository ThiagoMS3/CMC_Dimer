#include "SysConf.h"

/*		Instructions of the SysConf class dependent on simulation type
 *
 * 		-- Periodic (Moessner)
 * 		-- Manual
 * 		-- Partition
 *
 * 		TODO : properly break it in parts
 *
 */

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
