#include "SysConf.h"

// *******************************************************************
// General cluster CMC methods
// ---> The descriptions are in the *.h file
// *******************************************************************

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

void SysConf::GetSublattice(vector<double>& output, vector<double> & NSiteCount)
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

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount,vector<double> & Ndimer,vector<double> & MeanSublattice,vector<complex<double> >& complexPhaseVector,vector<double>& realPhaseVector,complex<double>& complexPhase,double& realPhase)
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
	{
		Ncount[iii] = 0;
	}

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

	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
	{
		GetSublattice(MeanSublattice,NSiteCount);
	}

	// Calculate the plaquette and columnar indexes
	int Sub0Value = 0;
	int Sub1Value = 0;
	int spinNeigh1 = 0;
	int spinNeigh0 = 0;
	int index0 = 0;
	int index1 = 0;

	for(int iii = 0; iii < N; ++iii)
	{
		complexPhaseVector[iii] = 0;
	}
	int chosenSublattice = -1;

	int lll = 0;
	int nnn = 0;
	for(int iii = 0; iii < N3Count; ++iii)
	{
		spinPos = N3_Map[iii];
		lll = spinPos/N;
		nnn = spinPos%N;

		// Complex phase
		chosenSublattice = SublatticePositions[lll];
		complexPhaseVector[nnn] += Jindex[chosenSublattice];
	}

	complexPhase = complexPhase/((double)N*L);
	realPhase = cos(3*arg(complexPhase));

	complexPhase = 0.;
	for(int iii = 0; iii < N; ++iii)
	{
		complexPhaseVector[iii] = complexPhaseVector[iii]/((double)L);
		realPhaseVector[iii] = cos(3*arg(complexPhaseVector[iii]));
		complexPhase += complexPhaseVector[iii];
	}
	complexPhase = complexPhase/((double)N);
	realPhase = cos(3*arg(complexPhase));
}

void SysConf::GetSiteNf(vector<double> & NSiteCount,vector<double> & Ncount,vector<double> & Ndimer,vector<double> & MeanSublattice)
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
	{
		Ncount[iii] = 0;
	}

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

	if(SimType.compare("Moessner")==0||SimType.compare("Manual")==0)
	{
		GetSublattice(MeanSublattice,NSiteCount);
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

void 	SysConf::GetN3Correlation(double& corr, vector<double>& Star3Network)
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

//void	SysConf::GetSzSzCorrelation(vector<double>& corr)
//{
//	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
//	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] = 0;
//	}
//
//	double firstSpin = 0;
//
//	for(int iii = 0; iii < L; ++iii)
//	{
//		firstSpin = spinConf[idxConv(N,iii,N/2)];
//		for(int nnn = 0; nnn < N/2; ++nnn)
//		{
//			corr[N/2+nnn] += firstSpin*spinConf[idxConv(N,iii,nnn)];
//		}
//
//		for(int nnn = N/2; nnn < N; ++nnn)
//		{
//			corr[nnn-N/2] 	+= firstSpin*spinConf[idxConv(N,iii,nnn)];
//		}
//	}
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] 		= corr[nnn]/L;
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
//		for(int iii = 0; iii < L; ++iii)
//		{
//			firstSpin = spinConf[idxConv(N,iii,N/2)];
//
//			for(int jjj = 0; jjj < NbOfNeights; ++jjj)
//			{
//				index = idxConv(NbOfNeights,iii,jjj);
//				posNeigh  = neighboursTable[index];
//
//
//			#if BORDER_TYPE == 4
//			// ANTI-PERIODIC BORDERS
//				firstDimer = (firstSpin*weightTable[index]*spinConf[idxConv(N,posNeigh,N/2)] + 1)/2.;
//			#else
//				firstDimer = (firstSpin*spinConf[idxConv(N,posNeigh,N/2)] + 1)/2.;
//			#endif
//
//
//				for(int nnn = 0; nnn < N/2; ++nnn)
//				{
//				#if BORDER_TYPE == 4
//				// ANTI-PERIODIC BORDERS
//					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#else
//					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#endif
//					corr[N/2+nnn] += firstDimer*dimer;
//				}
//
//				for(int nnn = N/2; nnn < N; ++nnn)
//				{
//				#if BORDER_TYPE == 4
//				// ANTI-PERIODIC BORDERS
//					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#else
//					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
//				#endif
//					corr[nnn-N/2] 	+= firstDimer*dimer;
//				}
//			}
//		}
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] 		= corr[nnn]/(NbOfNeights*L);
//	}
//}
//
//void	SysConf::GetSxSxCorrelation(vector<double>& corr)
//{
//	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
//	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] = 0;
//	}
//
//	char firstDummy = 0;
//	char secondDummy = 0;
//
//	int firstSx = 0;
//	int layerSx = 0;
//
//	for(int iii = 0; iii < L; ++iii)
//	{
//		firstDummy = spinConf[idxConv(N,iii,N/2)];
//		secondDummy = -spinConf[idxConv(N,iii,N/2+1)];
//		if( firstDummy == secondDummy )
//		{
//			firstSx = 1;
//		}
//		else
//		{
//			firstSx = 0;
//		}
//
//		// Cover *most* of the layers
//		for(int nnn = 0; nnn < N/2; ++nnn)
//		{
//			firstDummy = spinConf[idxConv(N,iii,nnn)];
//			secondDummy = -spinConf[idxConv(N,iii,nnn+1)];
//			if( firstDummy == secondDummy )
//			{
//				layerSx = 1;
//			}
//			else
//			{
//				layerSx = 0;
//			}
//			corr[N/2+nnn] += firstSx*layerSx;
//		}
//
//		for(int nnn = N/2; nnn < N - 1; ++nnn)
//		{
//			firstDummy = spinConf[idxConv(N,iii,nnn)];
//			secondDummy = -spinConf[idxConv(N,iii,nnn+1)];
//			if( firstDummy == secondDummy )
//			{
//				layerSx = 1;
//			}
//			else
//			{
//				layerSx = 0;
//			}
//			corr[nnn-N/2] += firstSx*layerSx;
//		}
//
//		// Deal with boundary case : nnn = N - 1
//		firstDummy = spinConf[idxConv(N,iii,N - 1)];
//		secondDummy = -spinConf[idxConv(N,iii,0)];
//		if( firstDummy == secondDummy )
//		{
//			layerSx = 1;
//		}
//		else
//		{
//			layerSx = 0;
//		}
//		corr[N/2 - 1] += firstSx*layerSx;
//	}
//
//	for(int nnn = 0; nnn < N; ++nnn)
//	{
//		corr[nnn] 		= corr[nnn]/(L*Kz*Kz); // Kz = DBeta
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

	int nStep = 60;
	int deltaN = N/nStep;
	for(int iii = 0; iii < L; ++iii)
	{
		for(int mmm = 0; mmm < N; mmm += deltaN)
		{
			firstSpin = spinConf[idxConv(N,iii,mmm)];
			for(int nnn = 0; nnn < mmm; ++nnn)
			{
				corr[N - mmm+nnn] += firstSpin*spinConf[idxConv(N,iii,nnn)];
			}

			for(int nnn = mmm; nnn < N; ++nnn)
			{
				corr[nnn-mmm] 	+= firstSpin*spinConf[idxConv(N,iii,nnn)];
			}
		}

	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] 		= corr[nnn]/(nStep*L);
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

	int nStep = 60;
	int deltaN = N/nStep;
	for(int iii = 0; iii < L; ++iii)
	{
		for(int mmm = 0; mmm < N; mmm += deltaN)
		{
			firstSpin = spinConf[idxConv(N,iii,mmm)];

			for(int jjj = 0; jjj < NbOfNeights; ++jjj)
			{
				index = idxConv(NbOfNeights,iii,jjj);
				posNeigh  = neighboursTable[index];


			#if BORDER_TYPE == 4
			// ANTI-PERIODIC BORDERS
				firstDimer = (firstSpin*weightTable[index]*spinConf[idxConv(N,posNeigh,mmm)] + 1)/2.;
			#else
				firstDimer = (firstSpin*spinConf[idxConv(N,posNeigh,mmm)] + 1)/2.;
			#endif


				for(int nnn = 0; nnn < mmm; ++nnn)
				{
				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#else
					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#endif
					corr[N - mmm+nnn] += firstDimer*dimer;
				}

				for(int nnn = mmm; nnn < N; ++nnn)
				{
				#if BORDER_TYPE == 4
				// ANTI-PERIODIC BORDERS
					dimer = (spinConf[idxConv(N,iii,nnn)]*weightTable[index]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#else
					dimer = (spinConf[idxConv(N,iii,nnn)]*spinConf[idxConv(N,posNeigh,nnn)] + 1)/2.;
				#endif
					corr[nnn-mmm] 	+= firstDimer*dimer;
				}
			}
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] 		= corr[nnn]/(nStep*NbOfNeights*L);
	}
}

void	SysConf::GetSxSxCorrelation(vector<double>& corr)
{
	//   corr[nnn] 		= Sum_iii [ <A_(iii,0) A_(iii,nnn*dB)>_iii ]
	//   siteMean[iii]	= Sum_nnn [ <A_(iii,nnn)>_nnn ]

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] = 0;
	}

	char firstDummy = 0;
	char secondDummy = 0;

	int firstSx = 0;
	int layerSx = 0;

	int nStep = 60;
	int deltaN = N/nStep;
	for(int iii = 0; iii < L; ++iii)
	{
		for(int mmm = 0; mmm < N; mmm += deltaN)
		{
			firstDummy = spinConf[idxConv(N,iii,mmm)];
			if(mmm != N - 1)
			{
				secondDummy = -spinConf[idxConv(N,iii,mmm+1)];
			}
			else
			{
				secondDummy = -spinConf[idxConv(N,iii,0)];
			}
			firstSx = ( firstDummy == secondDummy );

			// Cover *most* of the layers
			for(int nnn = 0; nnn < mmm; ++nnn)
			{
				firstDummy = spinConf[idxConv(N,iii,nnn)];
				secondDummy = -spinConf[idxConv(N,iii,nnn+1)];
				layerSx = ( firstDummy == secondDummy );
				corr[N - mmm+nnn] += firstSx*layerSx;
			}

			for(int nnn = mmm; nnn < N - 1; ++nnn)
			{
				firstDummy = spinConf[idxConv(N,iii,nnn)];
				secondDummy = -spinConf[idxConv(N,iii,nnn+1)];
				layerSx = ( firstDummy == secondDummy );
				corr[nnn-mmm] += firstSx*layerSx;
			}

			// Deal with boundary case : nnn = N - 1
			firstDummy = spinConf[idxConv(N,iii,N - 1)];
			secondDummy = -spinConf[idxConv(N,iii,0)];
			layerSx = ( firstDummy == secondDummy );
			corr[N - 1 - mmm] += firstSx*layerSx;
		}
	}

	for(int nnn = 0; nnn < N; ++nnn)
	{
		corr[nnn] 		= corr[nnn]/(nStep*L*Kz*Kz); // Kz = DBeta
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
