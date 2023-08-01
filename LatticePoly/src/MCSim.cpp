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
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataWriter.h>


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
	
	NbindedCohesin =  0;
	active_forks =  0;
	binded_forks =  0;
	
	
	
	
	
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

	if ( (frame == 0) && (polyType != "MCPoly") )
		static_cast<MCHeteroPoly*>(pol)->BuildHetTable();
	
	acceptCountPoly = 0;
	NbindedCohesin = (polyType == "MCReplicPoly") ?  static_cast<MCReplicPoly*>(pol)->NbindedCohesin : 0;
	active_forks = (polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol)->activeForks.size() : 0;
	binded_forks = (polyType == "MCReplicPoly") and Jf_sister!=0 ?  static_cast<MCReplicPoly*>(pol)->NbindedForks : 0;

	//two different enhancement according to the topology
	for ( int i = 0; i < pol->Ntad + enhancement_cohesin*NbindedCohesin + enhancement_fork* (active_forks- binded_forks) + enhancement_sister*binded_forks ; ++i )
	{

		if ( frame < Nrelax + NG1 or 0==0)
			UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCPoly*>(pol), &acceptCountPoly);
		
		else
			UpdateTAD<>(lat, pol, &acceptCountPoly);
	}
	
	
	
	acceptAvePoly += acceptCountPoly / ((double) pol->Ntad);
	
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

	if ( frame >= Nrelax + NG1 )
	{
		if ( latticeType == "MCLattice" )
			UpdateRepl(static_cast<MCLattice*>(lat), pol);
		else
			UpdateRepl(static_cast<MCLiqLattice*>(lat), pol);

	}
		
	++cycle;

	
	

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
		
	if ( frame == Nrelax + Nmeas)
		pol->PrintCohesins();

	lat->ToVTK(frame);

	pol->ToVTK(frame);
	
	
	/*if ( frame == Nrelax + Nmeas)
	{
		vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();

		// Set the file name
		std::string path = outputDir + "output.vtm";
		
		writer->SetFileName(path.c_str());
		// Set the input of the writer to the vtkMultiBlockDataSet
		// configure the writer to not create vtp files
		writer->SetCompressorTypeToNone();
		writer->SetDataModeToBinary();
		writer->SetWriteMetaFile(0);
		writer->SetEncodeAppendedData(0);

		writer->SetInputData(pol->mbds);

		writer->Write();

		

		
	}*/

}


template class MCSim<MCLattice, MCPoly>;

template class MCSim<MCLattice, MCHeteroPoly>;
template class MCSim<MCLattice, MCLivingPoly>;
template class MCSim<MCLattice, MCReplicPoly>;

template class MCSim<MCLiqLattice, MCHeteroPoly>;
template class MCSim<MCLiqLattice, MCLivingPoly>;
template class MCSim<MCLiqLattice, MCReplicPoly>;
