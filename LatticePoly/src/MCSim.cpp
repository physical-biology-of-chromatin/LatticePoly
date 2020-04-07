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
	cycle = 0;
	
	acceptAveLiq = 0.;
	acceptAvePoly = 0.;
	
	tStart = std::chrono::high_resolution_clock::now();
	
	InitRNG();
		
	lat->Init(rngEngine);
	pol->Init(rngEngine);
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

	rngEngine.seed(0);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	double polyRate = acceptAvePoly / ((long double) cycle);
	std::cout << "Polymer acceptance rate: " << 100*polyRate << "%" << std::endl;
	
	if ( latticeType == "MCLiqLattice" || latticeType == "MCCGLattice" )
	{
		double liqRate = acceptAveLiq / ((long double) cycle);
		std::cout << "Liquid acceptance rate: " << 100*liqRate << "%" << std::endl;
	}
	
	tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::ratio<60,1>> tElapsed = tEnd - tStart;
	
	std::cout << "Total runtime: " << tElapsed.count() << " mins (" << cycle/tElapsed.count() << " MC cycles/min)" << std::endl;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Run()
{
	acceptCountLiq = 0;
	acceptCountPoly = 0;
	
	if ( latticeType == "MCLiqLattice" || latticeType == "MCCGLattice" )
	{
		int NliqMoves = std::ceil(NliqMC * Ldens * Ntot);
		
		MCCGLattice* liqlat = static_cast<MCCGLattice*>(lat);
		
		for ( int i = 0; i < Nchain; i++ )
			UpdateTAD<>(liqlat, pol, rngEngine, rngDistrib, &acceptCountPoly);
		
		for ( int i = 0; i < NliqMoves; i++ )
			UpdateSpin<>(liqlat, pol, rngEngine, rngDistrib, &acceptCountLiq);
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);
	}

	else
	{
		for ( int i = 0; i < Nchain; i++ )
			UpdateTAD<>(lat, pol, rngEngine, rngDistrib, &acceptCountPoly);
	}
	
	acceptAvePoly += acceptCountPoly / ((double) Nchain);
	cycle++;
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

template class MCSim<MCCGLattice, MCPoly>;
template class MCSim<MCCGLattice, MCHeteroPoly>;
