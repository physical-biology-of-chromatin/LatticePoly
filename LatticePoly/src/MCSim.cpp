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
	
	pol0 = new polymer(lat);
	pol1 = new polymer(lat);
	pol2 = new polymer(lat);
	pol3 = new polymer(lat);
	pol4 = new polymer(lat);
	pol5 = new polymer(lat);
	pol6 = new polymer(lat);
	pol7 = new polymer(lat);
	pol8 = new polymer(lat);
	pol9 = new polymer(lat);
	pol10 = new polymer(lat);
	pol11 = new polymer(lat);
	pol12 = new polymer(lat);
	pol13 = new polymer(lat);
	pol14 = new polymer(lat);
	pol15 = new polymer(lat);
	pol16 = new polymer(lat);



	pol_yeast={
		pol0 ,
		pol1 ,
		pol2 ,
		pol3,
		pol4 ,
		pol5 ,
		pol6,
		pol7 ,
		pol8 ,
		pol9 ,
		pol10 ,
		pol11 ,
		pol12 ,
		pol13,
		pol14,
		pol15,
		pol16,

	};
}

template<class lattice, class polymer>
MCSim<lattice, polymer>::~MCSim()
{
	delete lat;
	delete pol1;
	delete pol2;
	delete pol3;
	delete pol4;


}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Init()
{
	InitRNG();
	InitSimRange();
	lat->Init(Ninit);
	
	std::vector<double> x_pos_chrom={-0.14747378,  0.03028516,  0.27210213, -0.5210614 ,  0.50625318,
		-0.17220201, -0.33236478,  0.72759515, -0.76223487,  0.423846  ,
		0.29928386, -0.86521121,  0.97667577, -0.57512943, -0.12851069,
		0.764649  , -0.99914605};
	
	std::vector<double> y_pos_chrom={0.13509806, -0.34508377,  0.35490905, -0.09216842, -0.32203683,
		0.64058291, -0.63994816,  0.26571658,  0.31463947, -0.90573427,
		0.95416412, -0.50140758, -0.21471943,  0.81806243, -0.99170812,
		0.64444698,  0.04131783};
	
	
	int chrom_pos[3]={0,0,0};

	std::vector<int> indexes;
	for( int i = 0; i < 17; i++ )
		indexes.push_back( i );
	std::shuffle (indexes.begin(), indexes.end(), lat->rngEngine);
	
	for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
	{

		chrom_pos[0]=int((0.6*(0.5*L)*x_pos_chrom[i])+ 0.5*L );
		chrom_pos[1]=int((0.6*(0.5*L)*y_pos_chrom[i])+ 0.5*L );
		chrom_pos[2]=(L/2);

		
		pol_yeast.at(indexes.at(i))->Init(Ninit,indexes.at(i),chrom_pos);

	}
	


		
	




	NliqMoves = (latticeType == "MCLattice") ? 0 : NliqMC * static_cast<MCLiqLattice*>(lat)->nLiq;
	
	cycle = 0;
	acceptAveLiq = 0.;
	acceptAvePoly = 0.;
	
	NbindedCohesin =  0;
	active_forks =  0;
	binded_forks =  0;
	NbindedCohesin_loops=0;
	
	
	
	
	
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
	
	if(polyType=="MCReplicPoly")
	{
		bool fullyRepl=true;
		for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
			if(pol_yeast.at(i)->Ntad!=2* static_cast<MCReplicPoly*>(pol_yeast.at(i))->individual_Nchain)
				fullyRepl=false;
		if(fullyRepl)
			throw std::runtime_error("finshed repli");
	}
		
		


		

	/*if ( (cycle == (unsigned long long) Nrelax*Ninter) && (polyType != "MCPoly") )
		for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
			static_cast<MCHeteroPoly*>(pol_yeast.at(i))->BuildHetTable();
	 */
	
	NbindedCohesin=0;
	active_forks=0;
	NbindedCohesin_loops=0;
	binded_forks=0;
	

	int N_moves=0;
	for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
		N_moves=N_moves+pol_yeast.at(i)->Ntad;
	for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
	{
		NbindedCohesin = NbindedCohesin + ((polyType == "MCReplicPoly") ?  static_cast<MCReplicPoly*>(pol_yeast.at(i))->NbindedCohesin : 0);
		active_forks = active_forks + ((polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol_yeast.at(i))->activeForks.size() : 0);
		binded_forks = binded_forks + ((polyType == "MCReplicPoly") and Jf_sister!=0 ?  static_cast<MCReplicPoly*>(pol_yeast.at(i))->NbindedForks : 0);
		NbindedCohesin_loops = NbindedCohesin_loops+((polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol_yeast.at(i))->active_extruders.size() : 0);
	}
		
	

	//two different enhancement according to the topology
	
	for ( int i = 0; i < N_moves + enhancement_cohesin*(NbindedCohesin+2*NbindedCohesin_loops) + enhancement_fork* (active_forks- binded_forks) + enhancement_sister*binded_forks ; ++i )
	{
		int t = lat->rngEngine() % (int) pol_yeast.size();


		if ( frame < Nrelax + NG1)
			UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCPoly*>(pol_yeast.at(t)), &acceptCountPoly);

		else
			UpdateTAD<>(lat, (pol_yeast.at(t)), &acceptCountPoly);

	}
	 
	
	
	acceptAvePoly += acceptCountPoly / ((double) N_moves);
	
	if ( latticeType != "MCLattice" )
	{
		acceptCountLiq = 0;
		

		for ( int i = 0; i < NliqMoves; ++i )
		{
			if ( frame < Nrelax )
				UpdateSpin<>(static_cast<MCLattice*>(lat), static_cast<MCPoly*>(pol0), &acceptCountLiq);
			else
				UpdateSpin<>(lat, pol0, &acceptCountLiq);
		}
		
		acceptAveLiq += acceptCountLiq / ((double) NliqMoves);

		
		
	}

	if ( frame >= Nrelax + NG1 )
	{
		if ( latticeType == "MCLattice" )
		{
			
			active_forks=0;
			std::vector<MCTad*> activeOrigins;
			std::vector<int> respective_chain;


			for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
			{
				active_forks = active_forks + ((polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol_yeast.at(i))->activeForks.size() : 0);
				auto origins =static_cast<MCReplicPoly*>(pol_yeast.at(i))->activeOrigins;
				for ( int j = 0; j < (int) origins.size() ; ++j )
				{
					activeOrigins.push_back(origins.at(j));
					respective_chain.push_back(i);
				}
			}

			if ( (int) activeOrigins.size() > 0 )
			{

				//auto originsCopy =activeOrigins;
				//auto respective_chainCopy =respective_chain;

				/*auto shuffle_seed=lat->rngEngine;
				std::shuffle (originsCopy.begin(), originsCopy.end(),shuffle_seed);
				std::shuffle (respective_chainCopy.begin(), respective_chainCopy.end(),shuffle_seed);
				*/
				
				for ( int i=0 ; i < (int)activeOrigins.size(); i++) //for every element in indexes
				{
					
					MCTad* origin = activeOrigins[i]; //select origin taf
					double rndReplic = lat->rngDistrib(lat->rngEngine);
					
					int Nocc = active_forks % 2 == 0 ? int(active_forks) : int(active_forks)+ 1;
					// -1 since origin firing implicate 2 new monomer in the system
					if ( rndReplic < exp(-double(cycle)/(5*60/0.0003)) * double(2*Ndf- Nocc) * originRate and origin->status==0)
					{
						auto chrom=respective_chain.at(i);
						//std::cout << "assign chrom" <<  chrom << std::endl;
						
						
						static_cast<MCReplicPoly*>(pol_yeast.at(chrom))->Replicate(origin);
						//std::cout << "end repli" << std::endl;

						
						active_forks=0;
						for ( int k = 0; k < (int) pol_yeast.size()  ; ++k )
							active_forks = active_forks + ((polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol_yeast.at(k))->activeForks.size() : 0);
						
					}
				}

				
			}
			auto shuffled_pol_yeast =pol_yeast;
			std::shuffle (shuffled_pol_yeast.begin(), shuffled_pol_yeast.end(), lat->rngEngine);

			for ( int i = 0; i < (int) shuffled_pol_yeast.size()  ; ++i )
					UpdateRepl(static_cast<MCLattice*>(lat), shuffled_pol_yeast.at(i));
		}
		else
		{
			auto shuffled_pol_yeast =pol_yeast;
			std::shuffle (shuffled_pol_yeast.begin(), shuffled_pol_yeast.end(), lat->rngEngine);
			
			for ( int i = 0; i < (int) shuffled_pol_yeast.size()  ; ++i )
				UpdateRepl(static_cast<MCLiqLattice*>(lat), shuffled_pol_yeast.at(i));
		}
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

	//if ( frame == Nrelax + Nmeas)
		//pol->PrintCohesins();
	
	lat->ToVTK(frame);

	for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
		pol_yeast.at(i)->ToVTK(frame,std::to_string(i));


	
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
template class MCSim<MCLattice, MCHeteroPoly_looped>;
