//
//  MCSim.cpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCSim.hpp"


template<class lattice, class polymer>
MCSim<lattice, polymer>::MCSim()
{
	lat = new lattice;
	pol = new polymer(lat);
}

template<class lattice, class polymer>
MCSim<lattice, polymer>::~MCSim()
{
	delete lat;
	delete pol;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Init()
{
	Nfinal = Nrelax + Nmeas;

	cycle = 0;
	acceptAveLiq = 0.;
	acceptAvePoly = 0.;

	InitRNG();
		
	lat->Init();
	pol->Init();
		
	tStart = std::chrono::high_resolution_clock::now();
	tCycle = std::chrono::high_resolution_clock::now();
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::InitRNG()
{
	int seed;
	
	FILE* tmp = fopen("/dev/urandom", "rb");
	
	if ( (tmp != NULL) && (fread((void*) &seed, sizeof(seed), 1, tmp) != 0) )
		std::cout << "Using entropy-harvested random seed: " << seed << std::endl;
	
	else
	{
		seed = (int) time(NULL);
		std::cout << "Using system time as RNG seed: " << seed << std::endl;
	}
	
	fclose(tmp);
	
	lat->rngEngine.seed(seed);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Run()
{
	acceptCountPoly = 0;
	
	for ( int i = 0; i < Nchain; i++ )
		UpdateTAD<>(lat, pol, &acceptCountPoly);
	
	if ( latticeType == "MCLiqLattice" )
	{
		acceptCountLiq = 0;
		
		int NliqMoves = std::ceil(NliqMC * Ldens * Ntot);
		
		for ( int i = 0; i < NliqMoves; i++ )
			UpdateSpin<>(lat, pol, &acceptCountLiq);
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);
	}
	
	acceptAvePoly += acceptCountPoly / ((double) Nchain);
	
	cycle++;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	std::cout << "************" << std::endl;
	std::cout << "Performed " << cycle << " out of " << Nfinal*Ninter << " MC cycles" << std::endl;

	double polyRate = acceptAvePoly / ((long double) Ninter);
	std::cout << "Polymer acceptance rate: " << 100*polyRate << "%" << std::endl;
	
	acceptAvePoly = 0;
	
	if ( latticeType == "MCLiqLattice" )
	{
		double liqRate = acceptAveLiq / ((long double) Ninter);
		std::cout << "Liquid acceptance rate: " << 100*liqRate << "%" << std::endl;
		
		acceptAveLiq = 0;
	}
	
	std::chrono::high_resolution_clock::time_point tInter = tCycle;
	tCycle = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double, std::ratio<60,1>> dTotal = tCycle - tStart;
	std::chrono::duration<double, std::ratio<1,1>>  dCycle = tCycle - tInter;

	std::cout << "Total runtime: " << dTotal.count() << " mins (" << Ninter/dCycle.count() << " cycles/s)" << std::endl;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::DumpVTK(int idx)
{
	lat->ToVTK(idx);
	pol->ToVTK(idx);
}


template class MCSim<MCLattice, MCPoly>;
template class MCSim<MCLattice, MCHeteroPoly>;

template class MCSim<MCLiqLattice, MCPoly>;
template class MCSim<MCLiqLattice, MCHeteroPoly>;
