#include "MCIsing.h"

lagged_fibonacci607 rng(42);
// lagged_fibonacci607 rng(time(0));

// --- Build spin configuration ---
// --------------------------------
// void rngLoop()
// {
// 	int count = 0;
// 	for(int iii = 0; iii < 500; ++iii)
// 	{
// 		count += GetRandomSpin(9)==6;
// 	}
// 	cout << (double)count/500 << endl;
// };

// void SetNeighbours(int L,					// Lattice dimensions
// 			    int LTot,
// 			    int NbOfNeights,			//
// 			    int nx,				//
// 			    int ny,				//
// 			    int p,				//
// 			    vector<coord> &offset,
// 			    idxTable &neighboursTable		// Output
// 			   )
// {
// 	// ---> Set up the neighbour relations of the spins.
// 	//	   First, we have to set up a translation table between the (x,y) coordinates and
// 	//	the list positions. This table will be implemented as a binary search tree ('map')
// 	
// 	// > Parameters
// 	int positionTable = 0;
// 	coord position = {0,0};
// 	vector<int> rowStart(p + nx + 1,-1);
// 	vector<int> rowEnd(p + nx + 1,-1);
// 
// 	// > Translation table
// 	map<coord,int> translate;
// 
// 	// > Set limits
// 	for(int iii = 0; iii < p; ++iii)
// 	{
// 		rowStart[iii] = 0;
// 	}
// 	for(int iii = p; iii < p + nx + 1; ++iii)
// 	{
// 		rowStart[iii] = iii - p;
// 	}
// 	
// 	for(int iii = 0; iii < nx + 1; ++iii)
// 	{
// 		rowEnd[iii] = ny + iii;
// 	}
// 	for(int iii = nx + 1; iii < p + nx + 1; ++iii)
// 	{
// 		rowEnd[iii] = nx + ny;
// 	}
// 	
// 	// > Set up flippable spins
// 	for(int iii = 1; iii < p + nx; ++iii)
// 	{
// 		position.x = iii;
// 		for(int jjj = rowStart[iii] + 1; jjj < rowEnd[iii]; ++jjj)
// 		{
// 			position.y = jjj;
// 			translate.insert(pair<coord,int>(position,positionTable));
// 			++positionTable;
// 		}
// 	}
// 	
// 	// > Set up borders
// 	// First line
// 	position.x = 0;
// 	for(int jjj = rowStart[0]; jjj < rowEnd[0]; ++jjj)
// 	{
// 		position.y = jjj;
// 		translate.insert(pair<coord,int>(position,positionTable));
// 		++positionTable;
// 	}
// 	
// 	// Right border
// 	for(int iii = 0; iii < p + nx; ++iii)
// 	{
// 		position.x = iii;
// 		position.y = rowEnd[iii];
// 		translate.insert(pair<coord,int>(position,positionTable));
// 		++positionTable;
// 	}
// 	
// 	// Last line
// 	position.x = p + nx;
// 	for(int jjj = rowEnd[p + nx]; jjj > rowStart[p + nx]; --jjj)
// 	{
// 		position.y = jjj;
// 		translate.insert(pair<coord,int>(position,positionTable));
// 		++positionTable;
// 	}
// 	
// 	// Left border
// 	for(int iii = p + nx; iii > 0; --iii)
// 	{
// 		position.x = iii;
// 		position.y = rowStart[iii];
// 		translate.insert(pair<coord,int>(position,positionTable));
// 		++positionTable;
// 	}
// 
// 	// > Now set neighbouring relations
// 	position.x = 0;
// 	position.y = 0;
// 	
// 	int index  = 0;
// 	for(int iii = 1; iii < p + nx; ++iii)
// 	{
// 		position.x = iii;
// 		for(int jjj = rowStart[iii] + 1; jjj < rowEnd[iii]; ++jjj)
// 		{
// 			position.y = jjj;
// 			index = translate[position];
// 			for(int kkk = 0; kkk < NbOfNeights; ++kkk)
// 			{
// 				neighboursTable[index][kkk] = translate[position + offset[kkk]];
// 			}
// 		}
// 	}
// 	
// 	if(ny==1)
// 	{
// 		for(int iii = 0; iii < nx; ++iii)
// 		{
// 			position.x = iii;
// 			position.y = rowEnd[iii];
// 			index = translate[position];
// 			neighboursTable[index][2] = translate[position + offset[2]];
// 		}
// 	}
// };

// void SetNeighboursMoessner(int nx,				// lines
// 			    int ny,				// cols
// 			    vector<int> &neighboursTable		// Output
// 			   )
// {
// 	// ---> Set up the neighbour relations of the spins.
// 	//	   For the periodic boundary conditions, the setup of the spins is a lot easier ...
// 	
// 	// Parameters
// 	int index = 0;
// 	
// 	// Set bulk
// 	for(int iii = 1; iii < nx-1; ++iii)
// 	{
// 		for(int jjj = 1; jjj < ny-1; ++jjj)
// 		{
// 			index = nx*jjj + iii;
// 			neighboursTable[idxConv(6,index,0)] = index - 1;
// 			neighboursTable[idxConv(6,index,1)] = index + nx - 1;
// 			neighboursTable[idxConv(6,index,2)] = index + nx;
// 			neighboursTable[idxConv(6,index,3)] = index + 1;
// 			neighboursTable[idxConv(6,index,4)] = index - nx + 1;
// 			neighboursTable[idxConv(6,index,5)] = index - nx;
// 		}
// 	}
// 	
// 	// Set left border --> change 4 & 5
// 	for(int iii = 1; iii < nx-1; ++iii)
// 	{
// 		index = iii;
// 		neighboursTable[idxConv(6,index,0)] = index - 1;
// 		neighboursTable[idxConv(6,index,1)] = index + nx - 1;
// 		neighboursTable[idxConv(6,index,2)] = index + nx;
// 		neighboursTable[idxConv(6,index,3)] = index + 1;
// 		neighboursTable[idxConv(6,index,4)] = index + (ny-1)*nx + 1;	// !!!
// 		neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;		// !!!
// 	}
// 	
// 	// Set right border --> change 1 & 2
// 	for(int iii = 1; iii < nx-1; ++iii)
// 	{
// 		index = ((ny-1)*nx + iii);
// 		neighboursTable[idxConv(6,index,0)] = index - 1;
// 		neighboursTable[idxConv(6,index,1)] = iii -1;			// !!!
// 		neighboursTable[idxConv(6,index,2)] = iii;			// !!!
// 		neighboursTable[idxConv(6,index,3)] = index + 1;
// 		neighboursTable[idxConv(6,index,4)] = index - nx + 1;
// 		neighboursTable[idxConv(6,index,5)] = index - nx;
// 	}
// 	
// 	// Set upper border --> change 0 & 1
// 	for(int jjj = 1; jjj < ny-1; ++jjj)
// 	{
// 		index = jjj*nx;
// 		neighboursTable[idxConv(6,index,0)] = index + nx - 1;		// !!!
// 		neighboursTable[idxConv(6,index,1)] = index + 2*nx - 1;		// !!!
// 		neighboursTable[idxConv(6,index,2)] = index + nx;
// 		neighboursTable[idxConv(6,index,3)] = index + 1;
// 		neighboursTable[idxConv(6,index,4)] = index - nx + 1;
// 		neighboursTable[idxConv(6,index,5)] = index - nx;
// 	}
// 
// 	// Set lower border --> change 3 & 4
// 	for(int jjj = 1; jjj < ny-1; ++jjj)
// 	{
// 		index = ((jjj + 1)*nx - 1);
// 		neighboursTable[idxConv(6,index,0)] = index - 1;
// 		neighboursTable[idxConv(6,index,1)] = index + nx - 1;
// 		neighboursTable[idxConv(6,index,2)] = index + nx;
// 		neighboursTable[idxConv(6,index,3)] = index - nx + 1;		// !!!
// 		neighboursTable[idxConv(6,index,4)] = index - 2*nx + 1;		// !!!
// 		neighboursTable[idxConv(6,index,5)] = index - nx;
// 	}
// 	
// 	// > Finally, set corners
// 	// Upper left --> change 0, 1, 4 & 5 
// 	index = 0;
// 	neighboursTable[idxConv(6,index,0)] = index + nx - 1;			// !!!
// 	neighboursTable[idxConv(6,index,1)] = index + 2*nx - 1;			// !!!
// 	neighboursTable[idxConv(6,index,2)] = index + nx;
// 	neighboursTable[idxConv(6,index,3)] = index + 1;
// 	neighboursTable[idxConv(6,index,4)] = index + (ny-1)*nx + 1;		// !!!
// 	neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;			// !!!
// 	
// 	// Lower left --> change 3, 4 & 5 
// 	index = (nx - 1);
// 	neighboursTable[idxConv(6,index,0)] = index - 1;
// 	neighboursTable[idxConv(6,index,1)] = index + nx - 1;
// 	neighboursTable[idxConv(6,index,2)] = index + nx;
// 	neighboursTable[idxConv(6,index,3)] = index - nx + 1;			// !!!
// 	neighboursTable[idxConv(6,index,4)] = index + (ny-2)*nx + 1;		// !!! !!!
// 	neighboursTable[idxConv(6,index,5)] = index + (ny-1)*nx;			// !!!
// 	
// 	// Upper right --> change 0, 1 & 2
// 	index = (ny-1)*nx;;
// 	neighboursTable[idxConv(6,index,0)] = index + nx - 1;			// !!!
// 	neighboursTable[idxConv(6,index,1)] = index - (ny-2)*nx - 1;		// !!! !!!
// 	neighboursTable[idxConv(6,index,2)] = index - (ny-1)*nx;			// !!!
// 	neighboursTable[idxConv(6,index,3)] = index + 1;
// 	neighboursTable[idxConv(6,index,4)] = index - nx + 1;
// 	neighboursTable[idxConv(6,index,5)] = index - nx;
// 	
// 	// Lower right --> change 1, 2, 3 & 4
// 	index = (nx*ny-1);
// 	neighboursTable[idxConv(6,index,0)] = index - 1;
// 	neighboursTable[idxConv(6,index,1)] = nx - 2;
// 	neighboursTable[idxConv(6,index,2)] = nx - 1;
// 	neighboursTable[idxConv(6,index,3)] = (ny-1)*nx;
// 	neighboursTable[idxConv(6,index,4)] = (ny-2)*nx;
// 	neighboursTable[idxConv(6,index,5)] = index - nx;
// };

// void SetInitialSpinConf(int N,
// 				    int L,
// 				    int LTot,
// 				    int NbOfNeights,
// 				    idxTable &neighboursTable,
// 				    
// 				    spinStack &spinConf,
// 				    idxList &flippableSpinLists,
// 				    idxList &flippableSpinPositions,
// 				    int &TotalOldNbF
// 				   )
// {
// 	// TODO Change these to global !!!
// 	bernoulli_distribution<> bernoulli(0.5);
// 	
// 	int dummyNbF = 0;
// 	
// 	int nbOfDoublesUp = 0;
// 	int nbOfDoublesDown = 0;
// 	
// 	int lastFlipped = 0;
// 	// ---> Set borders
// 	for(int iii = L; iii < LTot; ++iii)
// 	{
// 		spinConf[iii][0] = pow(-1,iii);
// 	}
// 	
// 	// --> Now set the other spins
// 	//	> Each triangle of the lattice, formed by spins 'spin0', 'spin1' and 'spin2',
// 	//	  must have a total magnetization equal to +-1/2. This way, we'll have only one
// 	//	  frustration per triangle and we'll obey the dimer relations in the hexagonal
// 	//	  lattice. We have to take this constraint into account when we choose a random
// 	//	  configuration.
// 	//	> Since each 'spin0' is inside 6 triangles, we have to do this test 6 times in the
// 	//	  worst case
// 	
// 	idxList	flippableCandidate;
// 	flippableCandidate.reserve(L);
// 	
// 	bool ThereIsAChoice = true;
// 	int	spin1 = 0, spin2 = 0;
// 	int 	spin0 = 0;
// 	
// 	while(spin0 < L)
// 	{
// 		ThereIsAChoice = true;
// 		nbOfDoublesUp = 0;
// 		nbOfDoublesDown = 0;
// 		
// 		// Test for the first triangle
// 		spin1 = neighboursTable[spin0][0];
// 		spin2 = neighboursTable[spin0][NbOfNeights - 1];
// 		if((spinConf[spin1][0]==spinConf[spin2][0])&&spinConf[spin2][0]!=0)
// 		{
// 			if(spinConf[spin1][0]==1)
// 			{
// 				++nbOfDoublesUp;
// 			}
// 			else
// 			{
// 				++nbOfDoublesDown;
// 			}
// 			
// 			// Then we're not free to choose the spin
// 			spinConf[spin0][0] = -spinConf[spin1][0];
// 			ThereIsAChoice = false;
// 		}
// 
// 		for(int jjj = 1; jjj < NbOfNeights; ++jjj)
// 		{
// 			spin1 = neighboursTable[spin0][jjj];
// 			spin2 = neighboursTable[spin0][jjj-1];
// 			if((spinConf[spin1][0]==spinConf[spin2][0])&&spinConf[spin2][0]!=0)
// 			{
// 				if(spinConf[spin1][0]==1)
// 				{
// 					++nbOfDoublesUp;
// 				}
// 				else
// 				{
// 					++nbOfDoublesDown;
// 				}
// 				
// 				// Then we're not free to choose the spin
// 				spinConf[spin0][0] = -spinConf[spin1][0];
// 				ThereIsAChoice = false;
// 			}
// 		}
// 
// 		if(ThereIsAChoice)
// 		{
// 			// Then we can choose an spin
// 			// ---> Also, this spin might be flippable, so append it to the
// 			//	   candidates list
// 			spinConf[spin0][0] = -1 + 2*bernoulli(rng);
// 			flippableCandidate.push_back(spin0);
// 			lastFlipped = spin0;
// 		}
// 		if(nbOfDoublesDown>0&&nbOfDoublesUp>0)
// 		{
// 			for(int iii = lastFlipped + 1; iii < spin0 + 1; ++iii)
// 			{
// 				spinConf[iii][0] = 0;
// 			}
// 			spinConf[lastFlipped][0] = -spinConf[lastFlipped][0];
// 			spin0 = lastFlipped + 1;
// 		}
// 		else
// 		{
// 			++spin0;
// 		}
// 	}
// 	
// 	// Find which spins are flippable
// 	bool IsFlippable = false;
// 	
// 	for(uint iii = 0; iii < flippableCandidate.size(); ++iii)
// 	{
// 		IsFlippable = TestFlippable(spinConf,0,neighboursTable[flippableCandidate[iii]],NbOfNeights);
// 		if(IsFlippable)
// 		{
// 			flippableSpinLists[dummyNbF] = flippableCandidate[iii];
// 			flippableSpinPositions[flippableCandidate[iii]] = dummyNbF;
// 			++dummyNbF;
// 		}
// 	}
// 	
// 	// > Finally, copy informations over all the layers
// 	int pos = 0;
// 	int offset = 0;
// 	int spin = 0;
// 	for(int iii = 1; iii < N; ++iii)
// 	{
// 		for(int kkk = 0; kkk < LTot; ++kkk)
// 		{
// 			spinConf[kkk][iii] = spinConf[kkk][0];
// 		}
// 		
// 		offset = iii*L;
// 		for(int jjj = 0; jjj < dummyNbF; ++jjj)
// 		{
// 			pos = iii*dummyNbF + jjj;
// 			spin = offset + flippableSpinLists[jjj];
// 			flippableSpinLists[pos] = spin;
// 			flippableSpinPositions[spin] = pos;
// 		}
// 	}
// 	
// 	TotalOldNbF = N*dummyNbF;
// };

// void SetInitialSpinConfMoessner(int N,
// 				    int L,
// 				int nx,
// 				int ny,
// 				    int NbOfNeights,
// 				    vector<int> const&neighboursTable,
// 				    
// 				    vector<char> &spinConf,
// 				    vector<int> &flippableSpinLists,
// 				    vector<int> &flippableSpinPositions,
// 				    int &TotalOldNbF
// 				   )
// {
// 	// TODO This version is kind of a cheat ... it creates the limit case for V/t -> - infty
// 	
// 	int dummyNbF = 0;
// 	int kkk;
// 	for(int nnn = 0; nnn < N; ++nnn)
// 	{
// 		for(int jjj = 0; jjj < ny/3; ++jjj)
// 		{
// 			for(int iii = 0; iii < nx/3; ++iii)
// 			{
// 				kkk = (3*jjj*nx + 3*iii)*N;
// 				spinConf[kkk + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + nnn;
// 				flippableSpinPositions[kkk + nnn] = dummyNbF;
// 				++dummyNbF;
// 				
// 				spinConf[kkk + N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + N + nnn;
// 				flippableSpinPositions[kkk + N + nnn] = dummyNbF;
// 				++dummyNbF;
// 				
// 				spinConf[kkk + 2*N  + nnn] = -1;
// 				
// 				kkk = ((3*jjj+1)*nx + 3*iii)*N;
// 				spinConf[kkk + nnn] = -1;
// 				
// 				spinConf[kkk + N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + N + nnn;
// 				flippableSpinPositions[kkk + N + nnn] = dummyNbF;
// 				++dummyNbF;
// 				
// 				spinConf[kkk + 2*N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + 2*N + nnn;
// 				flippableSpinPositions[kkk + 2*N + nnn] = dummyNbF;
// 				++dummyNbF;
// 				
// 				kkk = ((3*jjj+2)*nx + 3*iii)*N;
// 				spinConf[kkk + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + nnn;
// 				flippableSpinPositions[kkk + nnn] = dummyNbF;
// 				++dummyNbF;
// 				
// 				spinConf[kkk + N + nnn] = -1;
// 				
// 				spinConf[kkk + 2*N + nnn] = 1;
// 				flippableSpinLists[dummyNbF] = kkk + 2*N + nnn;
// 				flippableSpinPositions[kkk + 2*N + nnn] = dummyNbF;
// 				++dummyNbF;
// 			}
// 		}
// 	}
// 
// 	TotalOldNbF = dummyNbF;
// };

// --- Observables / energy ---
// ----------------------------
void CountNf(vector<char> 	const &spinConf,
	     int 	 	N,
	     int		L,
	     
	     vector<int> 	const &neighbours,
	     
	     int 	 	NbOfNeights,
	     vector<double> 	&Ncount
	    )
{
	int test = 0;
	int kkk;
	int index = 0;
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = 0;
	
	for(int lll = 0; lll < L; ++lll)
	{
		for(int nnn = 0; nnn < N; ++nnn)
		{
			test = 0;
			for(int iii = 0; iii < NbOfNeights;++iii)
			{
				index = idxConv(NbOfNeights,lll,iii);
				test += (int)spinConf[idxConv(N,neighbours[index],nnn)];
			}
			
			switch(test)
			{
				case 0: 
					Ncount[3] += 1;
					break;
				case 2:
					Ncount[2]  += 1;
					break;
				case -2:
					Ncount[2]  += 1;
					break;
				case 4:
					Ncount[1]  += 1;
					break;
				case -4:
					Ncount[1]  += 1;
					break;
				case 6:
					Ncount[0]  += 1;
					break;
				case -6:
					Ncount[0]  += 1;
					break;
			}
		}
	}
	
	for(int iii = 0; iii < 4; ++iii)
		Ncount[iii] = Ncount[iii]/N;
};

void CalculateStagMag(spinStack& spinConf,
		      idxTable& neighboursTable,
		      int LTot, int LBorder,int L, int N,
		      int nx, int ny, int p,
		      double& stagMag)
{
	// For the sake of simplicity, we suppose ny = 1;
	// Idea : we can get the equivalent XXZ spin chain by searching the horizontal frustrations of our CMC Ising model.
	//   

	vector<int> spinChain(nx+p,1);
	vector<char> equivPart(nx,0);
	
	int pos = 0;
	int index = 0;
	stagMag = 0;
	
	for(int nnn = 0; nnn < N; ++nnn)
	{
		// Set the first nx spins as positive
		for(int jjj = 0; jjj < nx + p; ++jjj)
		{
			spinChain[jjj] = -1;
		}
		
				// Recover equivalent partition
		for(int jjj = 0; jjj < nx;++jjj)
		{
			equivPart[jjj] = p;
			pos = L + ny + jjj;
			
			for(int kkk = 0; kkk < p; ++kkk)
			{
				// The neighbour below SHOULD be neighboursTable[L][2] (see hex offset)
				if(spinConf[pos][nnn]==spinConf[neighboursTable[pos][2]][nnn])
				{
					// Then we got the horizontal dimer / frustrated relations
					break;
				}
				else
				{
					--equivPart[jjj];
					pos = neighboursTable[pos][2];
				}
			}
		}
		
		// Set up spin chain
		for(int jjj = 0; jjj < nx; ++jjj)
		{
			index = nx - 1 + equivPart[jjj] - jjj;
			spinChain[index] = 1;
		}
		
		// Calculate the staggered magnetization
		for(int jjj = 0; jjj < nx + p; ++jjj)
		{
			stagMag += pow(-1,jjj)*spinChain[jjj];
		}
	}
	
	stagMag = stagMag/((nx+p)*N);
};

// --- Get random spin / build cluster ---
// ---------------------------------------
void SpinUpdate(spinStack& spinConf, int nnn, int spinToFlip, int TotalNewNbF, int& TotalOldNbF,
			 idxList& flippableSpinLists,idxList& flippableSpinPositions, idxList& toAdd, idxList& toRemove,
			 int nbToAdd, int nbToRemove)
{
	// > Parameters
	int nbToExchange = min(nbToRemove,nbToAdd);
	int posToChange = 0;
	int spinToMove = 0;
	int LastPosition = TotalOldNbF - 1;
	
	// > Spin conf update
	spinConf[spinToFlip][nnn] = -spinConf[spinToFlip][nnn];
	TotalOldNbF = TotalNewNbF;
	
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
		posToChange = flippableSpinPositions[toRemove[iii]];

		// > Exchange the spins
		flippableSpinLists[posToChange] = toAdd[iii];

		// > Update 'flippableSpinPositions'
		flippableSpinPositions[toRemove[iii]] = -1;
		flippableSpinPositions[toAdd[iii]] = posToChange;
	}

	// > Add elements
	// * If we're on the case 2), we won't enter this loop
	//   because nbToExchange == nbToAdd
	for(int iii = nbToExchange; iii < nbToAdd; ++iii)
	{
		// > Update the LastPosition
		++LastPosition;
		
		// > Add the spin
		flippableSpinLists[LastPosition] = toAdd[iii];

		// > Update 'flippableSpinPositions'
		flippableSpinPositions[toAdd[iii]] = LastPosition;
	}

	// > Remove elements and fill holes
	// * If we're on the case 3), we won't enter this loop
	//   because nbToExchange == nbToRemove
	for(int iii = nbToExchange; iii < nbToRemove; ++iii)
	{
		// > Get position to remove
		posToChange = flippableSpinPositions[toRemove[iii]];
		
		// > Remove spin
		flippableSpinLists[posToChange] = -1;
		flippableSpinPositions[toRemove[iii]] = -1;
		
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
};

void ClusterUpdate(vector<char>& spinConf,
		   int kkk,
		   int clusterSize,
		   int clusterStart,
		   int clusterEnd,
		   int N,

		   int TotalNewNbF,
		   int& TotalOldNbF,

		   idxList& flippableSpinLists,
		   idxList& flippableSpinPositions, 
		   
		   idxTable& toAdd, 
		   idxTable& toRemove,
		   idxList& nbToAdd,
		   idxList& nbToRemove
		  )
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
		spinConf[idxConv(N,kkk,nnn)] = -spinConf[idxConv(N,kkk,nnn)];		
		
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
};

void TestUpdate(vector<char>& spinConf, int nnn, int N, int kkk, int L, int NbOfNeights, vector<int> const& neighboursTable, idxList& flippableSpinPositions, vector<int>& toRemove, int& nbToRemove, vector<int>& toAdd, int& nbToAdd)
{
	int ChosenNeigh = 0;
	int spin = 0;
	bool IsFlippable = false;
	
	nbToRemove = 0;
	nbToAdd = 0;
	spinConf[idxConv(N,kkk,nnn)] = -spinConf[idxConv(N,kkk,nnn)];
	for(int iii = 0; iii < NbOfNeights; ++iii)
	{
		
		ChosenNeigh = neighboursTable[idxConv(NbOfNeights,kkk,iii)];
		
		if(ChosenNeigh<L)
		{
			spin = idxConv(N,ChosenNeigh,nnn);
			if(flippableSpinPositions[spin]!=-1)
			{
				// Then the spin is flippable in the old configuration, but won't be
				//	anymore if the new conf is accepted
				// >>> Mark for removal
				toRemove[nbToRemove] = spin;
				++nbToRemove;
			}
			else
			{
				// Then the spin can't be flipped in the old configuration, but MIGHT
				//	be in the new one
				// >>> Have to test
				IsFlippable = TestFlippable(spinConf,nnn,N,neighboursTable,ChosenNeigh,NbOfNeights);
				
				if(IsFlippable)
				{
					toAdd[nbToAdd] = spin;
					++nbToAdd;
				}
			}
		}
	}
	spinConf[idxConv(N,kkk,nnn)] = -spinConf[idxConv(N,kkk,nnn)];
};

void BuildCluster(// Chosen spin
		  int nnn,
		  int kkk,
		  char firstSpin,
		  
		  // Energy parameters
		  double Kz,
		  double Vt,
		  double baseProb,
		  
		  // Geometry
		  int L,
		  int N, 
		  int NbOfNeights,
		  int TotalOldNbF,
		  vector<int> const & neighboursTable,
		  
		  // Spin configuration
		  idxList& flippableSpinLists,
		  idxList& flippableSpinPositions,
		  vector<char>& spinConf,

		  // Updated parameters
		  int& TotalNewNbF,
		  idxTable& toAdd,
		  idxList& nbToAdd,
		  idxTable& toRemove,
		  idxList& nbToRemove,
		  
		  int& clusterStart,
		  int& clusterEnd, 
		  int& clusterSize,
		  int& UpperBorderIsF,
		  int& LowerBorderIsF
		 )
{
	// Initialize parameters
	TotalNewNbF 		= TotalOldNbF;

	bool 	CanBuildCluster	= true;
	bool	IsSpinAdded 	= true;
	
	double	addProb		= 0;
	
	clusterStart 		= nnn;
	clusterEnd 		= nnn;
	clusterSize		= 1;
	UpperBorderIsF		= 1;
	LowerBorderIsF		= 1;
	
	// > Calculate the changes created by choosing this first spin
	TestUpdate(spinConf,nnn,N,kkk,L,NbOfNeights,neighboursTable,flippableSpinPositions,
		   toRemove[nnn],nbToRemove[nnn],toAdd[nnn],nbToAdd[nnn]
		  );

	TotalNewNbF += nbToAdd[nnn]-nbToRemove[nnn];
	
	// > Create cluster
	// First upwards
	while(CanBuildCluster)
	{
		++clusterEnd;
		if(clusterEnd==N)
		{
			// Then we looped the boundary conditions
			clusterEnd = 0;
		}
		
		// Set up would-be new configuration
		TestUpdate(spinConf,clusterEnd,N,kkk,L,NbOfNeights,neighboursTable,flippableSpinPositions,toRemove[clusterEnd],nbToRemove[clusterEnd],
			toAdd[clusterEnd],nbToAdd[clusterEnd]);
		
		if(flippableSpinPositions[idxConv(N,kkk,clusterEnd)]==-1)
		{	// Then we can't flip the spin
			UpperBorderIsF = 0;
		}
		
		if(clusterSize==N||spinConf[idxConv(N,kkk,clusterEnd)]!=firstSpin||UpperBorderIsF==0)
		{
			// Then we looped the boundary conditions and met the start of the cluster, 
			//	OR we've hit a different spin
			//	OR the spin can't be flipped
			CanBuildCluster = false;
		}
		else
		{
			// Add spin to cluster with probability addProb
// 			addProb = min(1,baseProb*exp(-Kz*Vt*(nbToAdd[clusterEnd]-nbToRemove[clusterEnd])/2.));
			addProb = baseProb*min(1,exp(-Kz*Vt*(nbToAdd[clusterEnd]-nbToRemove[clusterEnd])));
			IsSpinAdded = GetRandomAccept(addProb);
			if(IsSpinAdded)
			{
				// Update spin configuration
				++clusterSize;
				TotalNewNbF += nbToAdd[clusterEnd]-nbToRemove[clusterEnd];
			}
			else
			{
				// We didn't want this new spin
				CanBuildCluster = false;
			}
		}
	}
	
	// Then downwards
	CanBuildCluster = true;
	while(CanBuildCluster)
	{
		--clusterStart;
		
		if(clusterStart==-1)
		{
			// Then we looped the boundary conditions
			clusterStart = N - 1;
		}

		// Set up would-be new configuration
		TestUpdate(spinConf,clusterStart,N,kkk,L,NbOfNeights,neighboursTable,flippableSpinPositions,toRemove[clusterStart],nbToRemove[clusterStart],
			toAdd[clusterStart],nbToAdd[clusterStart]);
		
		if(flippableSpinPositions[idxConv(N,kkk,clusterStart)]==-1)
		{	// Then we can't flip the spin
			LowerBorderIsF = 0;
		}
		
		if(clusterSize==N||spinConf[idxConv(N,kkk,clusterStart)]!=firstSpin||LowerBorderIsF==0)
		{
			// Then we looped the boundary conditions and met the end of the cluster, 
			//	OR we've hit an different spin
			//	OR the spin can't be flipped
			CanBuildCluster = false;
		}
		else
		{
			// Add spin to cluster with probability addProb
// 			addProb = min(1,baseProb*exp(-Kz*Vt*(nbToAdd[clusterStart]-nbToRemove[clusterStart])/2.));
			addProb = baseProb*min(1,exp(-Kz*Vt*(nbToAdd[clusterStart]-nbToRemove[clusterStart])));
			IsSpinAdded = GetRandomAccept(addProb);
			if(IsSpinAdded)
			{
				// Update spin configuration
				++clusterSize;
				TotalNewNbF += nbToAdd[clusterStart]-nbToRemove[clusterStart];
			}
			else
			{
				// We didn't want this new spin
				CanBuildCluster = false;
			}
		}
	}
};

// // --- Bunching / binning ---
// // --------------------------
// void BunchingErrors(vector<double>& data,
// 		    int order,
// 		    int NbMeasures,
// 		    vector<double>& errData
// 		   )
// {
// 	int NbOfPoints = NbMeasures;
// 	vector<double> dataIn(data);
// 	vector<double> dataOut;
// 	
// 	for(int iii = 0; iii < order;++iii)
// 	{
// 		NbOfPoints = NbOfPoints/2;
// 		dataOut.resize(NbOfPoints);
// 		
// 		Bunch(dataIn,dataOut,errData[iii],NbOfPoints);
// 		dataIn.resize(NbOfPoints);
// 		dataIn = dataOut;
// 	}
// };
// 
// void Bunch(vector<double>& dataIn,
// 	   vector<double>& dataOut,
// 	   double& errData,
// 	   int NbOfPoints
// 	  )
// {
// 	double Sum = 0;
// 	double SumSqr = 0;
// 	
// 	for(int iii = 0; iii < NbOfPoints; ++iii)
// 	{
// 		Sum += dataIn[2*iii]+dataIn[2*iii + 1];
// 		SumSqr += dataIn[2*iii]*dataIn[2*iii]+dataIn[2*iii + 1]*dataIn[2*iii + 1];
// 		dataOut[iii] = (dataIn[2*iii]+dataIn[2*iii + 1])/2;
// 	}
// 
// 	errData = sqrt(abs((SumSqr/(2.*NbOfPoints)-pow(Sum/(2.*NbOfPoints),2))/(2.*NbOfPoints)));
// };

// void DataBunch(double meas,vector<double> & dataBin, int & binSize, int minBin, int & binCounter, int & fillBin)
// {
// 	// Add measurement to one of the bins
// 	dataBin[binCounter] += meas/binSize;
// 	++fillBin;
// 	
// 	if(fillBin==binSize)
// 	{	// ----> Then we filled a bin and have to go to the next one
// 		++binCounter;
// 		fillBin = 0;
// 		
// 		if(binCounter==2*minBin)
// 		{	// ----> Then there's no 'next bin' -> we have to compress the data!
// 		
// 			for(int iii = 0; iii < minBin; ++iii)
// 			{
// 				// > Compress
// 				dataBin[iii] = (dataBin[2*iii] + dataBin[2*iii+1])/2;
// 			}
// 			
// 			// > Resize bins
// 			binSize = 2*binSize;
// 				
// 			// Reset bins
// 			binCounter = minBin;
// 			for(int iii = minBin; iii < 2*minBin; ++iii)
// 			{
// 				dataBin[iii] = 0;
// 			}
// 		}
// 	}
// };

// // --- Autocorrelations ---
// // ------------------------
// void CalculateAutocorrelation(vector<double>& autocorr, double& integAutocorr, vector<double>& data, int dataSize, int acSize)
// {
// 	double Qti = 0;
// 	double Qi = 0;
// 	double Qt = 0;
// 	
// 	int N;
// 	double var = 0;
// 	integAutocorr = 0.5;
// 	// VERIFIED
// 	// Autocorrelation for the point 0 = variance
// 	for(int kkk = 0; kkk < dataSize; ++kkk)
// 	{
// 		Qti += data[kkk]*data[kkk];
// 		Qi += data[kkk];
// 	}
// 	Qti = Qti/dataSize;
// 	Qi = Qi/dataSize;
// 	
// 	var = (Qti - Qi*Qi);
// 	autocorr[0] = 1;
// 	
// 	// Now, autocorrelation for the other points
// 	for(int offset = 1; offset < acSize; ++offset)
// 	{
// 		N = dataSize - offset;
// 		autocorr[offset] = 0;
// 		Qti = 0;
// 		Qi = 0;
// 		Qt = 0;
// 		
// 		for(int kkk = 0; kkk < N; ++kkk)
// 		{
// 			Qti += data[kkk]*data[kkk+offset];
// 			Qi += data[kkk];
// 			Qt += data[kkk+offset];
// 		}
// 		Qti = Qti/N;
// 		Qi = Qi/N;
// 		Qt = Qt/N;
// 		
// 		// Autocorrelation
// 		autocorr[offset] = (Qti - Qt*Qi)/var;
// 		
// 		// Integrated autocorrelation
// 		integAutocorr += autocorr[offset];
// 	}
// };

// --- Merge files ---
// -------------------
void MergeMCOutput(const char filenameBase[], string& commonFile, int nodes, int VtN, int NValue, int erase)
{
	string inputFilename;
	string outputFilename = commonFile + "Final_" + filenameBase;
	ifstream inputs;
	ofstream output;
	
	// Set up V/t
	vector<int> localVtN(nodes,VtN/nodes);
	vector<int> localIndexStart(nodes,0);
	
	double Vt = 0;
	double value = 0;
	
	int extraVtN = VtN%nodes;
	for(int iii = 0; iii < extraVtN; ++iii)
	{
		++localVtN[iii];
	}
	
	output.open(outputFilename.c_str(),ios::trunc);
	
	for(int rank = 0; rank < nodes; ++rank)
	{
		inputFilename = commonFile + "r" + ToString(rank)  + "_" + filenameBase;
		inputs.open(inputFilename.c_str());

		for(int iii = 0; iii <localVtN[rank];++iii)
		{
			inputs >> Vt;
			output << Vt;
			for(int jjj = 0; jjj < NValue; ++jjj)
			{
				inputs >> value;
				output << " " << value;
			}
			output << endl;
		}
		
		output.flush();
		inputs.close();
	}
	output.close();
	
	if(erase)
	{
		for(int rank = 0; rank < nodes; ++rank)
		{
			inputFilename = commonFile + "r" + ToString(rank)  + "_" + filenameBase;
			remove(inputFilename.c_str());
		}
	}
};

// void MergeIntegACOutput(char filenameBase[], string& commonFile, int nodes, double VtStart, int VtN, double dVt, int data, int erase)
// {
// 	string inputFilename;
// 	string outputFilename = commonFile + "Final_" + filenameBase;
// 	ifstream inputs;
// 	ofstream output;
// 	
// 	vector<vector<double> > dataBuffer(VtN,vector<double>(data,0));
// 	
// 	// Set up V/t
// 	vector<int> localVtN(nodes,VtN/nodes);
// 	vector<int> localIndexStart(nodes,0);
// 	
// 	double Vt;
// 	double buffer;
// 	
// 	int extraVtN = VtN%nodes;
// 	for(int iii = 0; iii < extraVtN; ++iii)
// 		++localVtN[iii];
// 	
// 	for(int rank = 0; rank < nodes; ++rank)
// 	{
// 		if(rank < extraVtN)
// 		{
// 			localIndexStart[rank] = localVtN[rank]*rank;
// 		}
// 		else
// 		{
// 			localIndexStart[rank] = extraVtN + localVtN[rank]*rank;
// 		}
// 	}
// 
// 	// ---> Read each file
// 	for(int rank = 0; rank < nodes; ++rank)
// 	{
// 		inputFilename = commonFile + "r" + ToString(rank) + "_" + filenameBase;
// 		inputs.open(inputFilename.c_str());
// 		for(int iii = 0; iii < localVtN[rank]; ++iii)
// 		{
// 			inputs >> buffer;
// 			for(int jjj = 0; jjj < data; ++jjj)
// 			{
// 				inputs >> dataBuffer[localIndexStart[rank]+iii][jjj];
// 			}
// 			Jump(inputs,1);
// 		}
// 		inputs.close();
// 	}
// 	
// 	// ---> Print it all in a single file
// 	output.open(outputFilename.c_str(),ios::trunc);
// 	for(int iii = 0; iii < VtN; ++iii)
// 	{
// 		Vt = VtStart + iii*dVt;
// 		output << Vt;
// 		for(int jjj = 0; jjj < data; ++jjj)
// 		{
// 			 output << " " << dataBuffer[iii][jjj];
// 		}
// 		output << endl;
// 	}
// 	
// 	output.close();
// 	
// 	if(erase)
// 	{
// 		for(int rank = 0; rank < nodes; ++rank)
// 		{
// 			inputFilename = commonFile + "r" + ToString(rank) + "_" + filenameBase;
// 			remove(inputFilename.c_str());
// 		}
// 	}
// };

void MergeACOutput(const char filenameBase[], string& commonFile, int nodes, int VtN, int AutocorrLength, int erase)
{
	string inputFilename;
	string outputFilename = commonFile + "Final_" + filenameBase;
	ifstream inputs;
	ofstream output;
	
	vector<vector<double> > dataBuffer(AutocorrLength,vector<double>(VtN,0));
	
	// Set up V/t
	vector<int> localVtN(nodes,VtN/nodes);
	vector<int> localIndexStart(nodes,0);

	int buffer;
	int index = 0;
	
	int extraVtN = VtN%nodes;
	for(int iii = 0; iii < extraVtN; ++iii)
		++localVtN[iii];
	
	for(int rank = 0; rank < nodes; ++rank)
	{
		if(rank < extraVtN)
		{
			localIndexStart[rank] = localVtN[rank]*rank;
		}
		else
		{
			localIndexStart[rank] = extraVtN + localVtN[rank]*rank;
		}
	}
	
	// ---> Read each file
	for(int rank = 0; rank < nodes; ++rank)
	{
		inputFilename = commonFile + "r" + ToString(rank) + "_" + filenameBase;
		inputs.open(inputFilename.c_str());
		for(int kkk = 0; kkk < localVtN[rank]; ++kkk)
		{
			for(int jjj = 0; jjj < AutocorrLength; ++jjj)
			{
				inputs >> dataBuffer[jjj][localIndexStart[rank] + kkk];
			}
			Jump(inputs,1);
		}
		inputs.close();
	}
	
	// ---> Print it all in a single file
	output.open(outputFilename.c_str(),ios::trunc);
	for(int iii = 0; iii < AutocorrLength; ++iii)
	{
		output << iii;
		for(int jjj =0; jjj < VtN; ++jjj)
		{
			 output << " " << dataBuffer[iii][jjj];
		}
		output << endl;
	}
	
	output.close();
	
	if(erase)
	{
		for(int rank = 0; rank < nodes; ++rank)
		{
			inputFilename = commonFile + "r" + ToString(rank) + "_" + filenameBase;
			remove(inputFilename.c_str());
		}
	}
};

void MergeBunchOutput(const char filenameBase[], string& commonFile, int nodes, double VtStart, int VtN, double dVt, int OrderBunch, int erase)
{
	string inputFilenameBase;
	string inputFilename;
	string outputFilename = commonFile + "Final_" + filenameBase + ".dat";
	ifstream inputs;
	ofstream output;
	
	vector<vector<double> > dataBuffer(OrderBunch,vector<double>(VtN,0));
	
	// Set up V/t
	vector<int> localVtN(nodes,VtN/nodes);
	vector<double> localVtStart(nodes,0);
	vector<int> localIndexStart(nodes,0);
	double Vt;
	vector<int> buffer(OrderBunch,0);
	int index = 0;
	
	int extraVtN = VtN%nodes;
	for(int iii = 0; iii < extraVtN; ++iii)
		++localVtN[iii];
	
	for(int rank = 0; rank < nodes; ++rank)
	{
		if(rank < extraVtN)
		{
			localVtStart[rank] = VtStart + localVtN[rank]*rank*dVt;
			localIndexStart[rank] = localVtN[rank]*rank;
		}
		else
		{
			localVtStart[rank] = VtStart + localVtN[rank]*rank*dVt + extraVtN*dVt;
			localIndexStart[rank] = extraVtN + localVtN[rank]*rank;
		}
	}
	
	// ---> Read each file
	for(int rank = 0; rank < nodes; ++rank)
	{
		inputFilenameBase = commonFile + "r" + ToString(rank) + "_" + filenameBase + "_";
		for(int iii = 0; iii < localVtN[rank]; ++iii)
		{
			Vt = localVtStart[rank] + iii*dVt;
			if(abs(Vt)>1E-4)
			{
				inputFilename = inputFilenameBase + ToString(Vt) + ".dat";
			}
			else
			{
				inputFilename = inputFilenameBase + "0.dat";
			}
			inputs.open(inputFilename.c_str());
			for(int jjj = 0; jjj < OrderBunch; ++jjj)
			{
				inputs >> buffer[jjj];
				inputs >> dataBuffer[jjj][localIndexStart[rank]+iii];
				Jump(inputs,1);
			}
			inputs.close();
		}
	}
	
	// ---> Print it all in a single file
	output.open(outputFilename.c_str(),ios::trunc);
	for(int iii = 0; iii < OrderBunch; ++iii)
	{
		output << buffer[iii];
		for(int jjj =0; jjj < VtN; ++jjj)
		{
			 output << " " << dataBuffer[iii][jjj];
		}
		output << endl;
	}
	
	output.close();
	
	if(erase)
	{
		// ---> Read each file
		for(int rank = 0; rank < nodes; ++rank)
		{
			inputFilenameBase = commonFile + "r" + ToString(rank) + "_" + filenameBase + "_";
			for(int iii = 0; iii < localVtN[rank]; ++iii)
			{
				Vt = localVtStart[rank] + iii*dVt;
				if(abs(Vt)>1E-4)
				{
					inputFilename = inputFilenameBase + ToString(Vt) + ".dat";
				}
				else
				{
					inputFilename = inputFilenameBase + "0.dat";
				}
				remove(inputFilename.c_str());
			}
		}
	}
	
};