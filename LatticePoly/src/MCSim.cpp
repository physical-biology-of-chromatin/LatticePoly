//
//  MCSim.cpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <cstring>
#include <dirent.h>
#include <algorithm>

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

	InitSimRange();
	InitRNG();

	lat->Init(Ninit);
	pol->Init(Ninit);
		
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
		
	lat->rngEngine.seed(seed);
	
	fclose(tmp);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::InitSimRange()
{
	int liqId = 0;
	int polyId = 0;
	
	if ( RestartFromFile )
	{
		DIR* dir;
		dirent* pdir;
		
		std::vector<std::string> files;
		dir = opendir(outputDir.c_str());
		
		while ( (pdir = readdir(dir)) )
		{
			char* tmp = std::strtok(pdir->d_name, ".");
			tmp = std::strtok(NULL, ".");
			
			if ( tmp != NULL )
			 {
				 if ( std::strcmp(tmp, "vtp") == 0 )
					 files.push_back(pdir->d_name);
			 }
		}
			
		std::sort(files.rbegin(), files.rend());
		
		auto polyFound = std::find_if(files.begin(), files.end(),
									  [](const std::string& s){return s.find("poly") != std::string::npos;});
		auto liqFound = std::find_if(files.begin(), files.end(),
									 [](const std::string& s){return s.find("liq") != std::string::npos;});
		
		if ( polyFound != files.end() )
			polyId = std::atoi(polyFound->c_str() + std::strlen("poly"));
		
		else
		{
			RestartFromFile = false;

			std::cout << "Could not locate any polymer configuration files in directory " << outputDir << " - starting fresh" << std::endl;
		}
		
		if ( (RestartFromFile) && (latticeType == "MCLiqLattice") )
		{
			if ( liqFound != files.end() )
				liqId = std::atoi(liqFound->c_str() + std::strlen("liq"));
			
			else
			{
				RestartFromFile = false;

				std::cout << "Could not locate any liquid configuration files in directory " << outputDir << " - starting fresh" << std::endl;
			}
		}
	}

	Ninit = (latticeType == "MCLiqLattice") ? std::min(polyId, liqId) : polyId;
	Nfinal = Nrelax + Nmeas;
	
	if ( Ninit >= Nfinal )
		throw std::runtime_error("MCSim: Found configuration file with index " + std::to_string(Ninit) + " higher than Nfinal");
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Run()
{
	acceptCountPoly = 0;
	
	for ( int i = 0; i < Nchain; ++i )
		UpdateTAD<>(lat, pol, &acceptCountPoly);
	
	if ( latticeType == "MCLiqLattice" )
	{
		acceptCountLiq = 0;
		
		int NliqMoves = std::ceil(NliqMC * Ldens * Ntot);
		
		for ( int i = 0; i < NliqMoves; ++i )
			UpdateSpin<>(lat, pol, &acceptCountLiq);
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);
	}
	
	acceptAvePoly += acceptCountPoly / ((double) Nchain);
	
	++cycle;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	std::cout << "************" << std::endl;
	std::cout << "Performed " << cycle << " out of " << (Nfinal-Ninit)*Ninter << " MC cycles" << std::endl;

	double polyRate = acceptAvePoly / ((long double) Ninter);
	acceptAvePoly = 0;

	std::cout << "Polymer acceptance rate: " << 100*polyRate << "%" << std::endl;
		
	if ( latticeType == "MCLiqLattice" )
	{
		double liqRate = acceptAveLiq / ((long double) Ninter);
		acceptAveLiq = 0;

		std::cout << "Liquid acceptance rate: " << 100*liqRate << "%" << std::endl;
	}
	
	std::chrono::high_resolution_clock::time_point tInter = tCycle;
	
	tCycle = std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double, std::ratio<60,1>> dTotal = tCycle - tStart;
	std::chrono::duration<double, std::ratio<1,1>>  dCycle = tCycle - tInter;

	std::cout << "Total runtime: " << dTotal.count() << " mins (" << Ninter/dCycle.count() << " cycles/s)" << std::endl;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::DumpVTK(int frame)
{
	lat->ToVTK(frame);
	pol->ToVTK(frame);
}


template class MCSim<MCLattice, MCPoly>;
template class MCSim<MCLattice, MCHeteroPoly>;

template class MCSim<MCLiqLattice, MCPoly>;
template class MCSim<MCLiqLattice, MCHeteroPoly>;
