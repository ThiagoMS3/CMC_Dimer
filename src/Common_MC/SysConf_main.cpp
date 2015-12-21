#include "SysConf.h"

/*		Main instructions of the SysConf class
 *
 * 		-- Constructors
 * 		-- System setup
 * 		-- RNG setup
 *
 * 		-- Common MC instructions
 * 		-- Configuration randomizers
 * 		-- Output and debug
 *
 * 		TODO : properly break it in parts
 *
 */

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
// MC methods - Common instructions
// *******************************************************************

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
