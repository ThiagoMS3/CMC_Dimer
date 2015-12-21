#include "SysConf.h"

/*		Monte-Carlo instructions of the SysConf class
 *
 * 		-- Diagonal N3 (legacy)
 * 		-- Diagonal N0 (legacy)
 * 		-- Mixed
 *
 * 		TODO : properly break it in parts
 *
 */

// *******************************************************************
// MC methods - Diagonal = N3 --- LEGACY
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
// MC methods - Diagonal = N0 --- LEGACY
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
