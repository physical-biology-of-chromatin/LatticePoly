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
	
	std::vector<double> x_pos_chrom={-0.15375203955642344,
		0.0315744613550785,
		0.28368608169152465,
		-0.5432440614587442,
		0.5278054227160841,
		-0.17953301093540425,
		-0.3465142390318366,
		0.7585703792380805,
		-0.9243455561378048,
		0.4238459950479107,
		0.2992838644448729,
		-0.86521120975323,
		0.976675773628176,
		-0.5751294291397393,
		-0.12851068979899324,
		0.7646489954560431};
	
	std::vector<double> y_pos_chrom={0.14084946290918265,
		-0.3597747017215528,
		0.37001825820131595,
		-0.0960922253710707,
		-0.3357466062041174,
		0.6678538519389444,
		-0.6671920813772709,
		0.277028685854142,
		0.38155640847493627,
		-0.9057342725556136,
		0.954164120307897,
		-0.5014075812324265,
		-0.21471942904125782,
		0.8180624302199665,
		-0.9917081236973845,
		0.6444469828838243};
	
	
	int chrom_pos[3]={0,0,0};

	std::vector<int> indexes;
	for( int i = 0; i < 16; i++ )
		indexes.push_back( i );
	std::shuffle (indexes.begin(), indexes.end(), lat->rngEngine);
	
	for ( int i = 0; i < (int) indexes.size()  ; ++i )
	{
		chrom_pos[0]=int((0.6*(0.5*L)*x_pos_chrom[i])+ 0.5*L );
		chrom_pos[1]=int((0.6*(0.5*L)*y_pos_chrom[i])+ 0.5*L );
		chrom_pos[2]=(L/2);
		pol_yeast.at(indexes.at(i))->Init(Ninit,indexes.at(i),chrom_pos);
	}
	
	
	
	/*
	pol0->Init(Ninit,0);
	pol1->Init(Ninit,1);
	pol2->Init(Ninit,2);
	pol3->Init(Ninit,3);
	pol4->Init(Ninit,5);
	pol5->Init(Ninit,5);
	pol6->Init(Ninit,6);
	pol7->Init(Ninit,7);
	pol8->Init(Ninit,8);
	pol9->Init(Ninit,9);
	pol10->Init(Ninit,10);
	pol11->Init(Ninit,11);
	pol12->Init(Ninit,12);
	pol13->Init(Ninit,13);
	pol14->Init(Ninit,14);
	pol15->Init(Ninit,15);
	pol16->Init(Ninit,16);

*/




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
	
	
		

	if ( (frame == Nrelax) && (polyType != "MCPoly") )
		for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
			static_cast<MCHeteroPoly*>(pol_yeast.at(i))->BuildHetTable();
	

	
		/*static_cast<MCHeteroPoly*>(pol0)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol1)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol2)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol3)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol4)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol5)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol6)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol7)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol8)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol9)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol10)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol11)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol12)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol13)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol14)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol15)->BuildHetTable();
		static_cast<MCHeteroPoly*>(pol16)->BuildHetTable();*/
	



	//to modify
	acceptCountPoly = 0;
	NbindedCohesin = (polyType == "MCReplicPoly") ?  static_cast<MCReplicPoly*>(pol0)->NbindedCohesin : 0;
	active_forks = (polyType == "MCReplicPoly") ?  (int) static_cast<MCReplicPoly*>(pol0)->activeForks.size() : 0;
	binded_forks = (polyType == "MCReplicPoly") and Jf_sister!=0 ?  static_cast<MCReplicPoly*>(pol0)->NbindedForks : 0;

	//std::cout << pol->Ntad + enhancement_cohesin*NbindedCohesin + enhancement_fork* (active_forks- binded_forks) + enhancement_sister*binded_forks<< std::endl;

	/*int N_moves=pol0->Ntad+pol1->Ntad+pol2->Ntad+pol3->Ntad+pol4->Ntad+pol5->Ntad+pol6->Ntad+pol7->Ntad+pol8->Ntad+pol9->Ntad+pol10->Ntad+pol11->Ntad+pol2->Ntad+pol13->Ntad+pol14->Ntad+pol15->Ntad+pol16->Ntad;*/
	int N_moves=0;
	for ( int i = 0; i < (int) pol_yeast.size()  ; ++i )
		N_moves=N_moves+pol_yeast.at(i)->Ntad;
	

	//two different enhancement according to the topology
	for ( int i = 0; i < N_moves  ; ++i )
	{
		int t = lat->rngEngine() % (int) pol_yeast.size();
		


		if ( frame < Nrelax + NG1 or 0==1)
		{

			UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCHeteroPoly*>(pol_yeast.at(t)), &acceptCountPoly);
			//UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCHeteroPoly*>(pol2), &acceptCountPoly);
			//UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCHeteroPoly*>(pol3), &acceptCountPoly);
			//UpdateTAD<>(static_cast<MCLattice*>(lat), static_cast<MCHeteroPoly*>(pol4), &acceptCountPoly);


		}

		
		else
		{
			UpdateTAD<>(lat, (pol_yeast.at(t)), &acceptCountPoly);

			//UpdateTAD<>(lat, pol2, &acceptCountPoly);
			//UpdateTAD<>(lat, pol3, &acceptCountPoly);
			//UpdateTAD<>(lat, pol4, &acceptCountPoly);


		}
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
