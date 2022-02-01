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
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

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
	InitRNG();
	InitSimRange();

	lat->Init(Ninit);
	pol->Init(Ninit);
		
	NliqMoves = (latticeType == "MCLattice") ? 0 : NliqMC * static_cast<MCLiqLattice*>(lat)->nLiq;
	
	cycle = 0;
	acceptAveLiq = 0.;
	acceptAvePoly = 0.;


	
	

	
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
void MCSim<lattice, polymer>::InitSimRange()
{
	int liqId = 0;
	int polyId = 0;
	
	if ( RestartFromFile )
	{
		dirent* pdir;
		std::vector<std::string> files;
		
		DIR* dir = opendir(outputDir.c_str());
		
		while ( (pdir = readdir(dir)) )
		{
			std::string fileName = pdir->d_name;
			size_t pos = fileName.find_last_of(".");
			
			if ( (pos != std::string::npos) && (fileName.substr(pos+1) == "vtp") )
				files.push_back(fileName.substr(0, pos));
		}
		
		closedir(dir);

		std::sort(files.rbegin(), files.rend());
		
		auto polyFind = std::find_if(files.begin(), files.end(),
									 [](const std::string& s){return s.find("poly") != std::string::npos;});
		auto liqFind = std::find_if(files.begin(), files.end(),
									[](const std::string& s){return s.find("liq") != std::string::npos;});
		
		if ( polyFind != files.end() )
			polyId = std::atoi(polyFind->c_str() + std::strlen("poly"));
		else
			RestartFromFile = false;

		if ( liqFind != files.end() )
			liqId = std::atoi(liqFind->c_str() + std::strlen("liq"));
		else
			RestartFromFile = RestartFromFile && (latticeType == "MCLattice");

		if ( !RestartFromFile )
			std::cout << "Could not locate required configuration files in directory " << outputDir << " - starting fresh" << std::endl;
	}

	Ninit = (latticeType == "MCLattice") ? polyId : std::min(polyId, liqId);
	Nfinal = Nrelax + Nmeas;
	
	if ( Ninit >= Nfinal )
		throw std::runtime_error("MCSim: Found configuration file with index " + std::to_string(Ninit) + " higher than Nfinal");
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Run()
{
	acceptCountPoly = 0;

	for ( int i = 0; i < (int)(0*pol->activeForks.size()+pol->Ntad); ++i )
		UpdateTAD<>(lat, pol, &acceptCountPoly);



	
		
	
	acceptAvePoly += acceptCountPoly / ((double) pol->Ntad);
	
	if ( latticeType != "MCLattice" )
	{
		acceptCountLiq = 0;
		
		for ( int i = 0; i < NliqMoves; ++i )
			UpdateSpin<>(lat, pol, &acceptCountLiq);
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);
		
		
	}
	MCTad* empty_tad = nullptr;
	if ( latticeType == "MCLattice" )
		pol->OriginMove(empty_tad);

		

	pol->ForkMove();
	
	++cycle;
	
	if(cycle==900000)
	{
		double rnd = lat->rngDistrib(lat->rngEngine);		std::ostringstream streamObj;
		streamObj << rnd;
		std::setprecision(9);
		std::string strObj = streamObj.str();
		
		std::ofstream outfile(strObj+"fork74.res", std::ios_base::app | std::ios_base::out);
		
		int reptation74[800];
		for ( int i = 0; i < (int) 800 ; ++i )
			reptation74[i]=0;
		
		for ( int i = 500000; i < (int) pol->nreptation74.size() ; ++i )
			++reptation74[pol->nreptation74.at(i)];
		for ( int i = 0; i < (int) 800 ; ++i )
			outfile << reptation74[i] << std::endl;
		
		rnd = lat->rngDistrib(lat->rngEngine);
		std::ostringstream streamObj1;
		streamObj1 << rnd;
		std::setprecision(9);
		std::string strObj1 = streamObj1.str();
		std::ofstream outfile1(strObj1+"fork73.res", std::ios_base::app | std::ios_base::out);

		int reptation73[800];
		for ( int i = 0; i < (int) 800 ; ++i )
			reptation73[i]=0;
		
		for ( int i = 500000; i < (int) pol->nreptation73.size() ; ++i )
			++reptation73[pol->nreptation73.at(i)];
		for ( int i = 0; i < (int) 800 ; ++i )
			outfile1 << reptation73[i] << std::endl;

		rnd = lat->rngDistrib(lat->rngEngine);
		std::ostringstream streamObj2;
		streamObj2 << rnd;
		std::setprecision(9);
		std::string strObj2 = streamObj2.str();
		std::ofstream outfile2(strObj2+"fork64.res", std::ios_base::app | std::ios_base::out);
		
		int reptation64[800];
		for ( int i = 0; i < (int) 800 ; ++i )
			reptation64[i]=0;
		
		for ( int i = 500000; i < (int) pol->nreptation64.size() ; ++i )
			++reptation64[pol->nreptation64.at(i)];
		for ( int i = 0; i < (int) 800 ; ++i )
			outfile2 << reptation64[i] << std::endl;
		
		rnd = lat->rngDistrib(lat->rngEngine);
		std::ostringstream streamObj3;
		streamObj3 << rnd;
		std::setprecision(9);
		std::string strObj3 = streamObj3.str();
		std::ofstream outfile3(strObj3+"fork44.res", std::ios_base::app | std::ios_base::out);
		
		int reptation44[800];
		for ( int i = 0; i < (int) 800 ; ++i )
			reptation44[i]=0;

		for ( int i = 500000; i < (int) pol->nreptation44.size() ; ++i )
			++reptation44[pol->nreptation44.at(i)];
		for ( int i = 0; i < (int) 800 ; ++i )
			outfile3 << reptation44[i] << std::endl;
		
		rnd = lat->rngDistrib(lat->rngEngine);
		std::ostringstream streamObj4;
		streamObj4 << rnd;
		std::setprecision(9);
		std::string strObj4 = streamObj4.str();
		std::ofstream outfile4(strObj4+"fork14.res", std::ios_base::app | std::ios_base::out);
		
		int reptation14[800];
		for ( int i = 0; i < (int) 800 ; ++i )
			reptation14[i]=0;
		
		for ( int i = 500000; i < (int) pol->nreptation14.size() ; ++i )
			++reptation14[pol->nreptation14.at(i)];
		for ( int i = 0; i < (int) 800 ; ++i )
			outfile4 << reptation14[i] << std::endl;
	}
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	std::cout << "************" << std::endl;
	std::cout << "Performed " << cycle << " out of " << (Nfinal-Ninit)*Ninter << " MC cycles" << std::endl;

	std::cout << "Polymer acceptance rate: " << 100*acceptAvePoly / ((long double) Ninter) << "%" << std::endl;
		
	acceptAvePoly = 0;

	if ( latticeType != "MCLattice" )
	{
		std::cout << "Liquid acceptance rate: " << 100*acceptAveLiq / ((long double) Ninter) << "%" << std::endl;
		
		acceptAveLiq = 0;
	}
	
	auto tInter = tCycle;
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
template class MCSim<MCLiqLattice, MCHeteroPoly>;

template class MCSim<MCLattice, MCReplicPoly>;
template class MCSim<MCLiqLattice, MCReplicPoly>;
