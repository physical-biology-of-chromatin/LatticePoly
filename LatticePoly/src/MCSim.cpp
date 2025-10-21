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
	lat  = new lattice;
	pol  = new polymer(lat);
	pol1 = new polymer(lat); //TWO CHAIN
	pol2 = new polymer(lat); 
	pol3 = new polymer(lat); 
	pol4 = new polymer(lat);
	pol5 = new polymer(lat);
	pol6 = new polymer(lat);
	pol7 = new polymer(lat);
	pol8 = new polymer(lat);
	pol9 = new polymer(lat);



}

template<class lattice, class polymer>
MCSim<lattice, polymer>::~MCSim()
{
	delete lat;
	delete pol;
	delete pol1;
	delete pol2;
	delete pol3;
	delete pol4;
	delete pol5;
	delete pol6;
	delete pol7;
	delete pol8;
	delete pol9;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Init()
{
	InitRNG();
	InitSimRange();

	lat->Init(Ninit);
	pol->Init(Ninit);
	pol1->Init(Ninit); 
	pol2->Init(Ninit);
	pol3->Init(Ninit);
	pol4->Init(Ninit);
	pol5->Init(Ninit);
	pol6->Init(Ninit);
	pol7->Init(Ninit);
	pol8->Init(Ninit);
	pol9->Init(Ninit);

	pol_list={pol, pol1, pol2, pol3, pol4, pol5, pol6, pol7, pol8, pol9};


	int sphcount = 0;
	int polycount = 0;

	for ( int vi = 0; vi < Ntot; ++vi )
	{
	 	if ( lat->bitTable[0][vi] >= 0 )
	 		sphcount++;
		if ( lat->bitTable[0][vi] > 0 )
	 		polycount++;
	}
	std::cout << "Polymer vol fraction " <<  double(polycount)/double (sphcount) << std::endl;
	std::cout << "Sphere " <<  double (sphcount) << std::endl;

		
	NliqMoves = (latticeType == "MCLattice") ? 0 : NliqMC * static_cast<MCLiqLattice*>(lat)->nLiq;
	
	cycle = 0;
	acceptAveLiq = 0.;
	acceptAvePoly = 0.;
	acceptAveTopo = 0.;
	
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
void MCSim<lattice, polymer>::Run(int frame)
{
	acceptCountPoly = 0;
	acceptCountPolyTopo = 0;
	/*if ( frame < Nrelax and J_ext > 0.)
	{
		for ( int i = 0; i < ( int(Nchain/5) - (int) pol->activeExtruders.size() - (int) pol1->activeExtruders.size() - (int) pol2->activeExtruders.size()  ); ++i ) //TWO CHAIN
		{
			double rndL = lat->rngDistrib(lat->rngEngine); //TWO CHAIN
			if ( rndL < 0.33 )
				pol->LoadExtruders();
			else if ( 0.33 <= rndL && rndL < 0.66 )
				pol1->LoadExtruders();	 
			else if ( 0.66 <= rndL && rndL < 0.99 )
				pol2->LoadExtruders();	
		}	
		double rnd = lat->rngDistrib(lat->rngEngine);
		if( rnd < extrusion )	
			pol->Extrusion();
		pol->UnloadExtruders();	

		double rnd1 = lat->rngDistrib(lat->rngEngine); //TWO CHAIN
		if( rnd1 < extrusion )	
			pol1->Extrusion();
		pol1->UnloadExtruders();	

		double rnd2 = lat->rngDistrib(lat->rngEngine); 
		if ( rnd2 < extrusion )
			pol2->Extrusion();
		pol2->UnloadExtruders();	
	}
	if ( frame == Nrelax - 1 and J_ext > 0.)
	{
		if( pol->activeExtruders.size()>0 )
		{
			for ( int i = 0; i < (int) pol->activeExtruders.size(); ++i )
			{
				MCTad* Barrier = pol->activeExtruders.at(i)->loops;
				pol->activeExtruders.at(i) -> isCohesin = 0;
				pol->activeExtruders.at(i) -> loops = 0;
				pol->activeExtruders.at(i) -> loopDir = -1;
				Barrier -> isBarrier = 0;
				Barrier -> loops = 0;	
			}
			pol->activeExtruders.erase(std::remove_if(pol->activeExtruders.begin(), pol->activeExtruders.end(), [](const MCTad* tadMono){return tadMono->isCohesin == 0;}), pol->activeExtruders.end());		
		}	

		if( pol1->activeExtruders.size()>0 ) //TWO CHAIN
		{
			for ( int i = 0; i < (int) pol1->activeExtruders.size(); ++i )
			{
				MCTad* Barrier = pol1->activeExtruders.at(i)->loops;
				pol1->activeExtruders.at(i) -> isCohesin = 0;
				pol1->activeExtruders.at(i) -> loops = 0;
				pol1->activeExtruders.at(i) -> loopDir = -1;
				Barrier -> isBarrier = 0;
				Barrier -> loops = 0;	
			}
			pol1->activeExtruders.erase(std::remove_if(pol1->activeExtruders.begin(), pol1->activeExtruders.end(), [](const MCTad* tadMono){return tadMono->isCohesin == 0;}), pol1->activeExtruders.end());		
		}

		if( pol2->activeExtruders.size()>0 ) //TWO CHAIN
		{
			for ( int i = 0; i < (int) pol2->activeExtruders.size(); ++i )
			{
				MCTad* Barrier = pol2->activeExtruders.at(i)->loops;
				pol2->activeExtruders.at(i) -> isCohesin = 0;
				pol2->activeExtruders.at(i) -> loops = 0;
				pol2->activeExtruders.at(i) -> loopDir = -1;
				Barrier -> isBarrier = 0;
				Barrier -> loops = 0;	
			}
			pol2->activeExtruders.erase(std::remove_if(pol2->activeExtruders.begin(), pol2->activeExtruders.end(), [](const MCTad* tadMono){return tadMono->isCohesin == 0;}), pol2->activeExtruders.end());		
		}			
	}	
	else if ( frame >= Nrelax  and J_ext > 0. )
	{
		for ( int i = 0; i < ( NExtruders - (int) pol->activeExtruders.size() -  (int) pol1->activeExtruders.size() ); ++i ) //TWO CHAIN
		{
			double rndL = lat->rngDistrib(lat->rngEngine); //TWO CHAIN
			if ( rndL < 0.33 )
				pol->LoadExtruders();
			else if ( 0.33 <= rndL && rndL < 0.66 )
				pol1->LoadExtruders();	 
			else if ( 0.66 <= rndL && rndL < 0.99 )
				pol2->LoadExtruders();	 
		}	
		double rnd = lat->rngDistrib(lat->rngEngine);
		if( rnd < extrusion )	
			pol->Extrusion();
		pol->UnloadExtruders();	

		double rnd1 = lat->rngDistrib(lat->rngEngine); //TWO CHAIN
		if( rnd1 < extrusion )	
			pol1->Extrusion();
		pol1->UnloadExtruders();

		double rnd2 = lat->rngDistrib(lat->rngEngine); 
		if ( rnd2 < extrusion )
			pol2->Extrusion();
		pol2->UnloadExtruders();		
					
	}*/
	for ( auto current_p : pol_list )
	{
	
	for ( int i = 0; i < current_p->Ntad; ++i )
	{
		if ( frame < Nrelax )
			UpdateNoTopo<>(lat, current_p, &acceptCountPoly);
		else
			UpdateTAD<>(lat, current_p, &acceptCountPoly, &acceptCountPolyTopo);
		}	
	}
	

	if ( latticeType != "MCLattice" )
	{
		acceptCountLiq = 0;
		
		for ( int i = 0; i < NliqMoves; ++i )
		{
			if ( frame < Nrelax )
				UpdateSpin<>(static_cast<MCLattice*>(lat), static_cast<MCPoly*>(pol), &acceptCountLiq);
			else
				UpdateSpin<>(lat, pol, &acceptCountLiq);
		}
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);
	}
		
	++cycle;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	std::cout << "************" << std::endl;
	std::cout << "Performed " << cycle << " out of " << (Nfinal-Ninit)*Ninter << " MC cycles" << std::endl;

	std::cout << "Polymer acceptance rate: " << 100*acceptAvePoly / ((long double) Ninter) << "%" << std::endl;

	std::cout << "Topo move acceptance rate: " << 100*acceptAveTopo  / ((long double) Ninter) << "%" << std::endl;
		
	acceptAvePoly = 0;
	acceptAveTopo = 0;

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
	pol->ToVTK(frame, "A");
	pol1->ToVTK(frame, "B"); 
	pol2->ToVTK(frame, "C");//TWO CHAIN
	pol3->ToVTK(frame, "D");
	pol4->ToVTK(frame, "E");
	pol5->ToVTK(frame, "F");
	pol6->ToVTK(frame, "G");
	pol7->ToVTK(frame, "H");
	pol8->ToVTK(frame, "I");
	pol9->ToVTK(frame, "J");
	
}


template class MCSim<MCLattice, MCPoly>;

template class MCSim<MCLattice, MCHeteroPoly>;
template class MCSim<MCLattice, MCLivingPoly>;
template class MCSim<MCLattice, MCReplicPoly>;

template class MCSim<MCLiqLattice, MCHeteroPoly>;
template class MCSim<MCLiqLattice, MCLivingPoly>;
template class MCSim<MCLiqLattice, MCReplicPoly>;
