//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <random>

#include "MCReplicPoly.hpp"


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCLivingPoly(_lat) {}

void MCReplicPoly::Init(int Ninit,int chrom, int chrom_pos[3])
{
	MCLivingPoly::Init(Ninit,chrom,chrom_pos);
	total_activated_cars=0;
	NbindedForks=0;
	NbindedCohesin=0;
	std::vector<int> lattice_neigh_load1={0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, 1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4, 4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9, 10, 10, 11, 12};
	std::vector<int> lattice_neigh_load2={0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  1,  3,  5,  7, 9, 11, 12,  2,  4,  6,  8, 10, 11, 12,  3,  5,  8, 10, 11,  4,  6, 7,  9, 12,  5, 10, 12,  6,  9, 11,  7,  9, 12,  8, 10, 11,  9, 11, 10, 12, 11, 12};
	
	for(int n=0; n< 55 ; ++n)
	{
		lattice_neigh1[n]=lattice_neigh_load1[n];
		lattice_neigh2[n]=lattice_neigh_load2[n];

	}
	
	//for ( int vi = 0; vi < Ntot; ++vi )
	//	ReplTable[0][vi] = 0;
	activeForks.reserve(individual_Nchain);
	
	
	

	binded_particles.reserve(Ndf);
	for (int i = 0; i < (int) Ndf; ++i)
		binded_particles.push_back({});

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	
	
	if ( !RestartFromFile or RestartFromFile )
	{
		if(n_barriers!=0)
		{
			std::ifstream Carsfile(CARpath);
			
			if ( !Carsfile.good() )
			throw std::runtime_error("MCReplicPoly: Couldn't open file " + CARpath);
			
			std::string line_cars;
			
			while ( std::getline(Carsfile, line_cars) )
			{
				std::istringstream ss(line_cars);
				
				float d1;
				float d2;

				if ( ss >> d1 >>d2)
				{

					if(int(d1)==chrom)
					{
						ChIP.push_back(d2);
					}
				}
			}
			std::cout <<(int) chrom <<  std::endl;

			std::cout <<(int) ChIP.size()<<  std::endl;
			std::cout <<Ntad<<  std::endl;

			if (Ntad != (int) ChIP.size() )
				throw std::runtime_error("Nchain and CARs size do not match");
			
			
			Carsfile.close();
		}
		
		if(StartFromPODLS==true)
		{
			std::ifstream PODLSfile(PODLSpath);
			

			if ( !PODLSfile.good() )
			throw std::runtime_error("MCReplicPoly: Couldn't open file " + PODLSpath);
			
			std::string line_podls;
			while ( std::getline(PODLSfile, line_podls) )
			{
				std::istringstream ss(line_podls);
				
				float d1;
				float d2;

				
				if ( ss >> d1 >>d2)
				{
					if(int(d1)==chrom)
					{
						PODLS.push_back(d2);
					}

				}
			}
			
			std::cout <<"Podls"<<PODLS.size() <<  std::endl;
			std::cout <<Ntad<<  std::endl;



			if (Ntad != (int) PODLS.size() )
			throw std::runtime_error("Nchain and PODLS size do not match");
			
			individual_Nchain=(int) PODLS.size();
			individual_Ndf=int(individual_Nchain*Ndf/1531);

			//load_rfd
			RFD.reserve(individual_Nchain);
			for (int i = 0; i < (int) individual_Nchain; ++i)
				RFD.push_back(0);


			PODLSfile.close();
		}
	}
	
	
	if(StartFromPODLS==true)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<> d(PODLS.begin(), PODLS.end());
		origins={};
		for(int n=0; n<int(Ntad/5); ++n)
		{
			
			int origin=d(gen);
			origins.push_back(origin);
		}
	}
	else
	{
		origins={};
		std::ifstream OriginFile(OriginsPath);
		
		std::cout <<"Uniform Origins"<<  std::endl;
		
		if ( !OriginFile.good() )
		throw std::runtime_error("MCReplicPoly: Couldn't open file " + OriginsPath);
		
		std::string line;
		while ( std::getline(OriginFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			
			if ( ss >> d1 )
			{
				if (d1 >Ntad )
					throw std::runtime_error("Nchain and origin ID do not match");
				origins.push_back(d1);
				
			}
		}
		OriginFile.close();

	}
	

	//load Origins in activeOrigins vector
	for (int i=0 ; i < (int) origins.size();++i)
		activeOrigins.push_back( &tadConf[origins[i]]);


	if(n_barriers>0)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<> d(ChIP.begin(), ChIP.end());

		std::ofstream outfile_car(outputDir+"/car.res", std::ios_base::app | std::ios_base::out);
		
		//If I saturate the CARs all CARs are boundaries
		int nonzerobin=0;
		for (int i=0 ; i < (int) ChIP.size();++i)
			if(ChIP.at(i)!=0)
				++nonzerobin;
		int chromosome_n_barries = int(n_barriers*Ntad/1531);
		std::cout <<"n_barriers"<< chromosome_n_barries <<std::endl;

		if(chromosome_n_barries>nonzerobin)
		{
			std::cout << "Saturated_CARS" <<std::endl;

			chromosome_n_barries=nonzerobin;
		}
			
			
		active_cars={};
		while((int) active_cars.size() < int(chromosome_n_barries) )
		{
			int car=d(gen);
			if(std::find(active_cars.begin(),active_cars.end(),car) == active_cars.end())
			   active_cars.push_back(car);
		}

		
		
		for (int i=0 ; i < (int) active_cars.size();++i)
			tadConf[active_cars[i]].isCAR=true;
		
		
		//std::cout << "Extruder before " << N_extruders <<std::endl;
		
		//individual_N_extruders= int(N_extruders*Ntad);
		
		//std::cout << "Extruder after " << individual_N_extruders <<std::endl;

	}
	loaded_mcms={};
	
	std::cout <<"Finished Repli initi"<<  std::endl;


}


void MCReplicPoly::TrialMove(double* dE)
{

	MCHeteroPoly::TrialMove(dE);

}

void  MCReplicPoly::OriginMove_explicit(const int spinTable[Ntot])
{
	if ( (int) activeOrigins.size() > 0 )
	{
		auto originsCopy =activeOrigins;
		std::shuffle (originsCopy.begin(), originsCopy.end(), lat->rngEngine);
		
		
		
		for ( int i=0 ; i < (int)originsCopy.size(); i++) //for every element in indexes
		{
			
			MCTad* origin = originsCopy[i]; //select origin taf
			
			for ( int v = 0; v < 13 ; ++v )
			{
				int pos =(v == 0) ?  origin->pos : lat->bitTable[v][origin->pos];
				if(spinTable[pos]>0)
				{
					double rndReplic = lat->rngDistrib(lat->rngEngine); //is it correct?
					if(rndReplic < originRate  and origin->status==0 and !origin->isFork())
					{
						//previously to match with implicit: if(rndReplic < Ntot*originRate/12  and origin->status==0 and !origin->isFork())
						//a particle can activate only one origin
						if(std::find(Spin_pos_toDelete.begin(),Spin_pos_toDelete.end(),pos) == Spin_pos_toDelete.end())
						{

							//Replicate(origin);
							//check if replication occurs
							if(origin->status!=0)
								Spin_pos_toDelete.push_back(pos);
							
							//in case of licensing
							Spin_pos_toDelete.push_back(pos);
							std::cout << origin->SisterID  << std::endl;
							loaded_mcms.push_back(origin->SisterID);

						}
						if ( (int) loaded_mcms.size() == Ndf)
						{
							std::ofstream mcms(outputDir+"/loaded_mcm.res", std::ios_base::app | std::ios_base::out);
							for ( int mcm = 0; mcm < (int) loaded_mcms.size() ; ++mcm )
								mcms << loaded_mcms[mcm] << std::endl;
							
							throw std::runtime_error("End of Licencing");
							
						}
					}
				}
			}
		}
		
	}
}
void  MCReplicPoly::OriginMove_implicit() //This function is not used anymore in the Yeast_full_genome branch since origin firing is governed in MCSim comsidering all chromosomes 
{
	if ( (int) activeOrigins.size() > 0 )
	{
		std::cout <<"Origin move"<<  std::endl;
		auto originsCopy =activeOrigins;
		std::shuffle (originsCopy.begin(), originsCopy.end(), lat->rngEngine);
		
		
		for ( int i=0 ; i < (int)originsCopy.size(); i++) //for every element in indexes
		{
			
			MCTad* origin = originsCopy[i]; //select origin taf
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			
			int Nocc = activeForks.size() % 2 == 0 ? int(activeForks.size()) : int(activeForks.size())+ 1;
			// -1 since origin firing implicate 2 new monomer in the system
			if ( rndReplic < double(2*individual_Ndf- Nocc) * originRate and origin->status==0 and  Ntad < individual_Nchain + int(stop_replication))//Ntad < Nchain -2 + int(Nchain * stop_replication))
			{
				//int origin_pos=std::find(tadConf.begin(),tadConf.end(),origin) == tadConf.end();	
				auto origin_pos = std::distance(tadConf.begin(), std::find_if(tadConf.begin(), tadConf.end(),[origin](const MCTad& t) { return &t == origin; }));
				Replicate(origin);
				//std::ofstream mcms(outputDir+"/fired_origins.res", std::ios_base::app | std::ios_base::out);
                                //mcms<<origin_pos <<individual_Nchain<< std::endl;
				std::cout <<"FIRED ORIGIN "<< origin_pos<<"CHAIN SIZE "<<individual_Nchain<<  std::endl;

			}
		}
	}
}


void MCReplicPoly::ForkMove()
{
	/*int tot=0;
	for ( int v = 0; v < Ntot ; ++v )
		tot=tot+ReplTable[0][v];
	if(tot!=int(activeForks.size())*55)
		std::cout <<tot<<  std::endl;
	*/
	if ( activeForks.size() > 0   )
	{
		auto activeForksCopy =activeForks;
		for ( int i=0 ; i < (int)activeForksCopy.size(); i++)
		{
			MCTad* fork = activeForks[i];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( fork->status==0 and rndReplic < replicRate and Ntad < individual_Nchain - 1  + int(stop_replication))//Nchain + int(Nchain*stop_replication) )
				Replicate(fork);
		}

	}
	
}


void MCReplicPoly::Replicate(MCTad* tad)
{

	if (tad->isRightEnd() || tad->isLeftEnd())
		return;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
		
	double rnd = lat->rngDistrib(lat->rngEngine);
	
	//origin replication
	if ( !tad->isFork() )
	{	
		//int origin_pos=std::find(tadConf.begin(),tadConf.end(),origin) == tadConf.end();      
        /*        auto origin_pos = std::distance(tadConf.begin(), std::find_if(tadConf.begin(), tadConf.end(),[tad](const MCTad& t) { return &t == tad; }));
                std::ofstream mcms(outputDir+"/fired_origins.res", std::ios_base::app | std::ios_base::out);
                mcms<<origin_pos <<" "<<individual_Nchain<< std::endl;
                std::cout <<"FIRED ORIGIN "<< origin_pos<<"CHAIN SIZE "<<individual_Nchain<<  std::endl;*/
		if(std::find(activeOrigins.begin(),activeOrigins.end(),tad) == activeOrigins.end())
			throw std::runtime_error("origin from another chromosome");

			


		// Can't replicate tad if it's already adjacent to a fork
		if ( nb1->isFork() || nb2->isFork() )
			return;
		
		else
		{
			// Replicate extremities at half the normal rate
			if ( nb1->isLeftEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}

			else
			{

				activeForks.push_back(nb1);
				UpdateReplTable(nb1);//increase energy around new fork
				nb1->binding_site = nb2;

			}
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0. )
					return;
			}
			
			else
			{
				
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
				nb2->binding_site = nb1;




			}
		}
		
	}
	
	//Update the center of mass if I am at the end of replication
	
	else //fork displacement
	{
		// Replicating left fork means displacing it to its left neighbor, so we need to check if it's already a fork or chain end
		if ( tad->isLeftFork() )
		{
			if ( nb1->isLeftFork() || nb1->isRightEnd() )
				// Probably should never happen, do nothing
				return;
			
			if ( nb1->isRightFork() || nb1->isLeftEnd() ) // SPECIAL CASE
			{
				// Merge forks/replicate extremities at half the normal rate
				if ( rnd < 0.5 )
					return;
				
				//MERGING
				//if(tad->binding_site->isFork())
				//	tad->binding_site->binding_site=nb2;
				
				//if(nb1->isCAR)
				//	TurnCohesive(nb1);

			}
		
			else
			{

				activeForks.push_back(nb1); // STANDARD CASE
				UpdateReplTable(nb1);//increase energy around new fork
				
				nb1->binding_site=tad->binding_site;
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb1;

				//if(nb2->isCAR)
				//	TurnCohesive(nb2);

			}
		}
		
		// Same for right forks
		else if ( tad->isRightFork() )
		{
			if ( nb2->isRightFork() || nb2->isLeftEnd() )
				return;
			
			if  ( nb2->isLeftFork() || nb2->isRightEnd() ) // SPECIAL CASE
			{
				if ( rnd < 0.5 )
					return;
				//MERGING
				//if(tad->binding_site->isFork())
				//	tad->binding_site->binding_site=nb2;
				
				//if(nb2->isCAR)
				//	TurnCohesive(nb2);

				
			}
		
			else
			{
				
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
				nb2->binding_site=tad->binding_site;
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb2;
				
				//if(nb1->isCAR)
				//	TurnCohesive(nb1);

			}
		}
		
		// Delete old forks
		auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
		activeForks.erase(fork);
		UpdateReplTable(tad);

		
		if ( nb1->isFork() || nb2->isFork() )
		{
			auto fork2 = std::find(activeForks.begin(), activeForks.end(), nb1->isFork() ? nb1 : nb2);
			activeForks.erase(fork2);
			if(nb1->isFork()){
				UpdateReplTable(nb1);//since it's a fork about to be replicated, energy decreases
			}
			else{
				UpdateReplTable(nb2);//since it's a fork about to be replicated, energy decreases
			}
			
		}
	}
	
	//if(Ntad==2*Nchain-2)
		//Update_rcms_before_separation();
		


	ReplicateTADs(tad);
	ReplicateBonds(tad);
	Update();


	
}

void MCReplicPoly::ReplicateTADs(MCTad* tad)
{
	
	MCTad tadReplic;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];

	// Replicate left end/fork, if applicable
	if ( nb1->isLeftEnd() || nb1->isRightFork() )
	{
		
		
		tadReplic = *nb1;
		tadConf.push_back(tadReplic);
		nb1->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb1);
		
		if(latticeType=="MCLiqLattice")
		{
		   //when merging two forks or replicating end (iff the opposite one is replicated), I create a partcle at fork pos
			if(nb1->isLeftEnd() and tadConf.at(individual_Nchain-1).status!=0)
				Spin_pos_toCreate.push_back(tad->pos);
			else if(nb1->isRightFork())
				Spin_pos_toCreate.push_back(tad->pos);
		}
		
		if(nb1->isCAR)
		{

			tadConf.back().isCAR=true;
			TurnCohesive(nb1);
			//std::cout <<  " TURN COHESIVE nb1 " << std::endl;


			
		}
		if(nb1->domain!=-1)
			tadConf.back().domain=nb1->domain;
		if(nb1->type!=-1)
			tadConf.back().type=nb1->type;
		if( nb1->insulator_type.size()!=0)
				tadConf.back().insulator_type=nb1->insulator_type;

			


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);
	
	
	if(tad->isCAR)
	{

		tadConf.back().isCAR=true;
		//std::cout <<  " TURN COHESIVE tad " << std::endl;

		TurnCohesive(tad);

	}
	
	if(tad->domain!=-1)
		tadConf.back().domain=tad->domain;
	if(tad->type!=-1)
		tadConf.back().type=tad->type;
	if(tad->insulator_type.size()!=0)
		tadConf.back().insulator_type=tad->insulator_type;

	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		

		tadReplic = *nb2;
		tadConf.push_back(tadReplic);
		nb2->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb2);
		
		if(latticeType=="MCLiqLattice")
		{
			//when merging two forks or replicating end (iff the opposite one is replicated), I create a partcle at fork pos

			if(nb2->isRightEnd() and tadConf.at(0).status!=0)
				Spin_pos_toCreate.push_back(tad->pos);
			else if(nb2->isLeftFork())
				Spin_pos_toCreate.push_back(tad->pos);

		}
		
		if(nb2->isCAR)
		{
			tadConf.back().isCAR=true;
			//std::cout <<  " TURN COHESIVE nb2 " << std::endl;
			TurnCohesive(nb2);

		}
		
		if(nb2->domain!=-1)
			tadConf.back().domain=nb2->domain;
		if(nb2->type!=-1)
			tadConf.back().type=nb2->type;
		if(nb2->insulator_type.size()!=0)
			tadConf.back().insulator_type=nb2->insulator_type;


	}
}

void MCReplicPoly::ReplicateBonds(MCTad* tad)
{
	MCBond bondReplic1;
	MCBond bondReplic2;

	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
	
	// Create/modify bonds between relevant neighbors and replicated tad
	MCBond* bond1 = tad->isRightFork() ? tad->bonds[2] : tad->bonds[0];
	MCBond* bond2 = tad->isLeftFork() ? tad->bonds[2] : tad->bonds[1];
	
	bondReplic1.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad : bond1->id1;
	bondReplic1.id2 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
				
	bondReplic1.dir = bond1->dir;

	bondReplic2.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
	bondReplic2.id2 = (nb2->isRightEnd() || nb2->isLeftFork()) ? Ntad+1 : bond2->id2;
		
	bondReplic2.dir = bond2->dir;
	
	// For right forks, update bond1 to link left (replicated) neighbor to new tad
	if ( tad->isRightFork() )
	{
		bond1->id2 = bondReplic1.id2;
		
		SetBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			SetBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		SetBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			SetBond(*bond3);
			UnsetFork(nb1);
		}
	}
	
	else
		tadTopo.push_back(bondReplic2);
}

void MCReplicPoly::UnsetFork(MCTad* tad)
{
	if ( tad->isFork() )
	{
		tad->bonds[2] = nullptr;
		tad->neighbors[2] = nullptr;
		
		--tad->links;
	}
}

void MCReplicPoly::Update()
{
	// Update bonds
	if ( (int) tadTopo.size() > Nbond )
	{
		for ( auto bond = tadTopo.begin()+Nbond; bond != tadTopo.end(); ++bond )
			SetBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	// Update tads
	if ( (int) tadConf.size() > Ntad )
	{
		for ( auto tad = tadConf.begin()+Ntad; tad != tadConf.end(); ++tad )
		{
			
			if(tadConf.at(tad->SisterID).isCentromere)
			{
				MCTad* binding_centromere = &tadConf.at(tad->SisterID);
				tad->isCentromere = true;
				tad->isCohesin = true;
				binding_centromere->isCohesin=true;

				tad->binding_site = binding_centromere;
				binding_centromere->binding_site = &(*tad);
			}
			
			if(tadConf.at(tad->SisterID).isrDNA)
				tad->isrDNA=true;
			
			if ( tad->type != -1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable_tads[tad->type][vi];
				}
			}
			if ( tad->domain != -1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable_domain[tad->domain][vi];
				}
			}
			if ( tad->insulator_type.size() != 0 )
			{
				for ( int type_id = 0; type_id< tad->insulator_type.size();  ++type_id )
				{
					for ( int v = 0; v < 13; ++v )
					{
						int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
						
						++hetTable_insulator[tad->insulator_type.at(type_id)][vi];
					}
				}
			}
			
			++lat->bitTable[0][tad->pos];
		}
		
		Ntad = (int) tadConf.size();
	}
	
	// Update origins
	activeOrigins.erase(std::remove_if(activeOrigins.begin(), activeOrigins.end(), [](const MCTad* tad){return tad->status != 0;}),
						  activeOrigins.end());

	
	// Update fork/origin counters
	Nfork = (int) activeForks.size();
	
	//check how many forks are binded to their sister
	NbindedForks = 0;
	
	 for (int i=0; i < (int) activeForks.size();++i)
	 {
		 if (activeForks.at(i)->binding_site->isFork())
			 ++NbindedForks;
	 }

	//Update RFD vector

	if((int)activeForks.size()>0)
	{

		for (int i = 0; i < (int) activeForks.size(); ++i)
		{	
			int fork_pos = (int) std::distance(tadConf.data(), activeForks.at(i));
			if(RFD.at(fork_pos)==0)
			{
				int direction = activeForks.at(i)->isLeftFork() ? -1 :1 ;
				RFD.at(fork_pos)= direction;
			}
		}

	}



	 //enlarge_box
	/*
	if(Ntad==Nchain+int(0.2*Nchain) or Ntad==Nchain+1+int(0.2*Nchain))
		for ( int vi = 0; vi < Ntot; ++vi )
			if(lat->bitTable[0][vi]==10)
				lat->bitTable[0][vi]=0;
	if(Ntad==Nchain+int(0.4*Nchain) or Ntad==Nchain+1+int(0.4*Nchain))
		for ( int vi = 0; vi < Ntot; ++vi )
			if(lat->bitTable[0][vi]==20)
				lat->bitTable[0][vi]=0;
	if(Ntad==Nchain+int(0.6*Nchain) or Ntad==Nchain+1+int(0.6*Nchain))
		for ( int vi = 0; vi < Ntot; ++vi )
			if(lat->bitTable[0][vi]==30)
				lat->bitTable[0][vi]=0;
	if(Ntad==Nchain+int(0.8*Nchain) or Ntad==Nchain+1+int(0.8*Nchain))
		for ( int vi = 0; vi < Ntot; ++vi )
			if(lat->bitTable[0][vi]==40)
				lat->bitTable[0][vi]=0;
	if(Ntad==Nchain+int(1*Nchain) or Ntad==Nchain+1+int(1*Nchain))
		for ( int vi = 0; vi < Ntot; ++vi )
			if(lat->bitTable[0][vi]==50)
				lat->bitTable[0][vi]=0;
				
	*/

		
	//After every replication update I search for cohesive cohesins partner
	if(cohesionMode!=1)
		Find_cohesive_CAR();
}

double MCReplicPoly::GetEffectiveEnergy() //chiedere Maxime
{
	if (tadTrial->isFork() )
	{
		
		double Etot = 0.;

		if ( Jf > 0.  )
			Etot=Etot+Jf*(lat->ReplTable[0][tadUpdater->vo]-lat->ReplTable[0][tadUpdater->vn]);


		if ( Jf_sister > 0.  and tadTrial->binding_site->isFork())
		{

			double Jsister_replisome1=0.0;
			double Jsister_replisome2=0.0;

			double old_dist=0.0;
			double new_dist=0.0;
			for ( int dir = 0; dir < 3; ++dir )
			{
				//Here two sister forks are created among two NN, I just need to put the two in the same box when they are at box boundaries
				double distance=lat->xyzTable[dir][tadUpdater->vo]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance) > L/2. )
				{
					double pbcShift = std::copysign(L, distance);
					distance -= pbcShift;
				}
					
				old_dist=old_dist+SQR(distance);
					
				double distance1=lat->xyzTable[dir][tadUpdater->vn]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance1) > L/2. )
				{
					double pbcShift = std::copysign(L, distance1);
					distance1 -= pbcShift;
				}
				new_dist=new_dist+SQR(distance1);
			}
			double thr_distance = (neigh==1) ? 2 : 0.5;

			Jsister_replisome1= old_dist<=thr_distance ? 1 : old_dist/2;
			Jsister_replisome2= new_dist<=thr_distance ? 1 : new_dist/2;
			
			Etot=Etot-Jf_sister*(Jsister_replisome1-Jsister_replisome2);
			

	
		}

		return MCHeteroPoly::GetEffectiveEnergy()+Etot;
	}

	if (tadTrial->isCohesin  or (!tadTrial->isCohesin and tadTrial->isCAR and tadTrial->binding_site== &tadConf.at(tadTrial->SisterID)))
	{

		double Etot = 0.;
		
		
		
		if ( Jpair > 0.  )
		{
			
			double Jpair_anchors1=0.0;
			double Jpair_anchors2=0.0;
			
			double old_dist=0.0;
			double new_dist=0.0;
			double distance_old_new[3];
			for ( int dir = 0; dir < 3; ++dir )
			{
				double distance=lat->xyzTable[dir][tadUpdater->vo]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance) > L/2. )
				{
					double pbcShift = std::copysign(L, distance);
					distance -= pbcShift;
				}
				
				old_dist=old_dist+SQR(distance);
				
				double distance1=lat->xyzTable[dir][tadUpdater->vn]-lat->xyzTable[dir][tadTrial->binding_site->pos];
				while ( std::abs(distance1) > L/2. )
				{
					double pbcShift = std::copysign(L, distance1);
					distance1 -= pbcShift;
				}
				new_dist=new_dist+SQR(distance1);
				
				double distance_old_new_dir=lat->xyzTable[dir][tadUpdater->vn]-lat->xyzTable[dir][tadUpdater->vo];
				while ( std::abs(distance_old_new_dir) > L/2. )
				{
					double pbcShift = std::copysign(L, distance_old_new_dir);
					distance_old_new_dir -= pbcShift;
				}
				distance_old_new[dir]=distance_old_new_dir;

			}
			// if they are already binded compute Etot
			double thr_distance = (neigh==1) ? 2 : 0.5;
			if(old_dist<=thr_distance or 1==1)
			{
				Jpair_anchors1= old_dist<=thr_distance ? 1 : old_dist/2;
				Jpair_anchors2= new_dist<=thr_distance ? 1 : new_dist/2;
				
				Etot=Etot-Jpair*(Jpair_anchors1-Jpair_anchors2);
			}
			else// if old distance is greater than thr_distance, meaning they did not make any binindg yet
			{


				//compute distance vector between
				auto conf = BuildUnfoldedConf();
				int id1=(int) std::distance(tadConf.data(), tadTrial);
				int id2=(int) std::distance(tadConf.data(), tadTrial->binding_site);
				
				old_dist = 0.0;
				new_dist = 0.0;

				
				for ( int dir = 0; dir < 3; ++dir )
				{
					old_dist=old_dist+SQR(conf[id1][dir]-conf[id2][dir]);
					new_dist=new_dist+SQR(conf[id1][dir]+distance_old_new[dir]-conf[id2][dir]);
				}

				Jpair_anchors1= old_dist<=thr_distance ? 1 : old_dist/2;
				Jpair_anchors2= new_dist<=thr_distance ? 1 : new_dist/2;
				
				Etot=Etot-Jpair*(Jpair_anchors1-Jpair_anchors2);
				
			}
			
			
			
			
		}
		return MCHeteroPoly::GetEffectiveEnergy()+Etot;
	}
	
	return 	MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::TurnCohesive(MCTad* tad)
{

	
	if( std::find(cohesive_CARs.begin(),cohesive_CARs.end(),tad) == cohesive_CARs.end())
	{
		//std::cout <<  " TURN COHESIVE " << std::endl;

		//std::cout <<  "status = "<< tad->status << std::endl;

		double rnd = lat->rngDistrib(lat->rngEngine);
		double activation_rate = ForkTableMode==0? keco1 : keco1*lat->ReplTable[0][tad->pos];
		//with rate keco, I turn the replicated active CAR cohesive
		if(rnd<activation_rate)
		{
			//I turn cohesive only one of the two sisters
			double rnd2 = lat->rngDistrib(lat->rngEngine);
			if(rnd2<0.5)
			{
				
				cohesive_CARs.push_back(tad);
				++total_activated_cars;
			}
			else
			{
				cohesive_CARs.push_back(&tadConf.at(tad->SisterID));
				++total_activated_cars;
			}
		}	
	}
}
void MCReplicPoly::Find_cohesive_CAR()
{
	//std::cout <<  "FIND COHESIVE"<< std::endl;

	//find a  non-symmetric partner of cohesive CAR
	if(cohesionMode!=1  and cohesionMode!=2)
	{
		if(cohesive_CARs.size()>1 )
		{
			//copy and shuffle vector as cohesive CARs are removed
			auto cohesive_CARs_copy=cohesive_CARs;
			std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);


			for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )//loop over all cohesive CARs
			{

				
				if(!cohesive_CARs_copy.at(i)->isCohesin)
				{

					auto Sister_CAR=&tadConf.at( cohesive_CARs_copy.at(i)->SisterID);
					if(Sister_CAR->isCohesin) //if the homologous is already a cohesive coesin avoid crossing
						return;
					
					auto tad_shifter= Sister_CAR;
					//Make a symmetrical binding
					double rnd_symm = lat->rngDistrib(lat->rngEngine);

					if(cohesionMode!=4) //turn to -1 when cohesionmode=0 or 3 to have the case of strictly non-symmetrical binding 
						rnd_symm = -1; 
					if(rnd_symm > 0.0) //homologous binding, this can be modified to have both non-symmetrical and symmetrical binding (non used in pubblication)
					{
						cohesive_CARs_copy.at(i)->isCohesin=true;
						Sister_CAR->isCohesin=true;
						Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
						cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
						cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
						NbindedCohesin+=2;
						//PrintCohesins();
					}
					else
					{
						//when cohesionmode=0 I can search for an active CAR in both directions 
						double rnd = lat->rngDistrib(lat->rngEngine);
						bool same_direction=false;
						if(cohesionMode==3)
							same_direction=true;
							
						if (same_direction==true)
							rnd = Sister_CAR->status==1 ? 0.0 : 0.6; // impose same direction
						if(rnd>0.5)
						{
							// random, if rnd >0.5 go right
							while( !tad_shifter->isFork() and !tad_shifter->isRightEnd() and !tad_shifter->isLeftEnd())
							{
								//the starting point for my search (symmetrical position of cohesive CAR)
								tad_shifter=tad_shifter->neighbors[0];


								if( tad_shifter->isCohesin and !tad_shifter->isFork()) //I am stopped if already a cohesin or is a not a fork
								{
									//I go to previous position as the car is occupied
									tad_shifter=tad_shifter->neighbors[1];
									
									Sister_CAR=tad_shifter;

									//establishement of the binding and cohesin status
									cohesive_CARs_copy.at(i)->isCohesin=true;
									Sister_CAR->isCohesin=true;
									Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
									cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;

									//I delete the cohesive CAR from the vector since has now a binding partner 
									cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
									NbindedCohesin+=2;
									
									//PrintCohesins();
									
									return;
								}
								
								if(tad_shifter->isCAR and !tad_shifter->isCohesin  and !tad_shifter->isFork()) // I found an active CAR to use as an anchor
								{	
				
									Sister_CAR=tad_shifter;
									
									cohesive_CARs_copy.at(i)->isCohesin=true;
									Sister_CAR->isCohesin=true;
									Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
									cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
									cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
									NbindedCohesin+=2;
									
									
									//PrintCohesins();
									
									return;
								}
							}
						}
						else
						{	//same as above but I shift in other direction
							while( !tad_shifter->isFork() and !tad_shifter->isRightEnd() and !tad_shifter->isLeftEnd())
							{
								tad_shifter=tad_shifter->neighbors[1];

								if(tad_shifter->isCohesin and !tad_shifter->isFork())
								{
									tad_shifter=tad_shifter->neighbors[0];
									
									Sister_CAR=tad_shifter;									
									cohesive_CARs_copy.at(i)->isCohesin=true;
									Sister_CAR->isCohesin=true;
									Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
									cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
									cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
									NbindedCohesin+=2;
									//std::cout <<  "FOUND + cohesin " << std::endl;
									//std::cout <<  "anchor found with status  " << Sister_CAR->status<< " and binding of status " << Sister_CAR->binding_site->status<< std::endl;
									
									//PrintCohesins();
									
									return;
								}
								
								if(tad_shifter->isCAR and !tad_shifter->isCohesin  and !tad_shifter->isFork())
								{	//check if CAR is cohesive
									
									Sister_CAR=tad_shifter;
									//cohesive_CARs_copy.at(i)->binding_site->binding_site=nullptr;
									
									cohesive_CARs_copy.at(i)->isCohesin=true;
									Sister_CAR->isCohesin=true;
									Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
									cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
									//auto del = std::find(cohesive_CARs.begin(), cohesive_CARs.end(), cohesive_CARs_copy.at(i));
									cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
									//cohesive_CARs.erase(del);
									NbindedCohesin+=2;
									
									
									//std::cout <<  "FOUND + " << std::endl;
									//std::cout <<  "anchor found with status  " << Sister_CAR->status<< " and binding of status " << Sister_CAR->binding_site->status<< std::endl;
									
									//PrintCohesins();
									
									return;
								}
							}
						}
					}
				}
			}
		}
	}

	if(cohesionMode==1) //handcuff model for cohesion: not used in the paper
	{
		if(cohesive_CARs.size()>1 and std::find(cohesive_CARs.begin(),cohesive_CARs.end(),tadTrial) != cohesive_CARs.end())
		{


			auto cohesive_CARs_copy=cohesive_CARs;
			std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);
			for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )
			{
				double rnd = lat->rngDistrib(lat->rngEngine);

				if(rnd<1 and cohesive_CARs_copy.at(i)!=tadTrial and !cohesive_CARs_copy.at(i)->isCohesin and cohesive_CARs_copy.at(i)->status!=tadTrial->status)
				{
					auto conf=BuildUnfoldedConf();
					double distance=0.0;
					int id1=(int) std::distance(tadConf.data(), cohesive_CARs_copy.at(i));
					int id2=(int) std::distance(tadConf.data(), tadTrial);

					for ( int dir = 0; dir < 3; ++dir )
						distance=distance+SQR(conf[id1][dir]-conf[id2][dir]);
					
					//std::cout <<  "Try to find interaction between  "<< id2<<" and "<<id1<< std::endl;

					if(distance<= 2.0/sqrt(2))
					{

						std::cout <<  "Distance = "<< distance<<std::endl;

						cohesive_CARs_copy.at(i)->isCohesin=true;
						tadTrial->isCohesin=true;
						tadTrial->binding_site=cohesive_CARs_copy.at(i);
						cohesive_CARs_copy.at(i)->binding_site=tadTrial;
						cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
						NbindedCohesin+=2;
						std::cout <<  "car found partner between "<< id2<<" and "<<id1<< std::endl;
						//PrintCohesins();
					}
				}
			}
		}
	}
	if(cohesionMode==2) //symmetrical cohesion
	{

		if(cohesive_CARs.size()>1 )
		{
			
			auto cohesive_CARs_copy=cohesive_CARs;
			std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);
			for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )
			{
				if(!cohesive_CARs_copy.at(i)->isCohesin and cohesive_CARs_copy.at(i)->status!=0) //check that is not already a cohesin and it has been replicated
				{
					auto Sister_CAR=&tadConf.at(cohesive_CARs_copy.at(i)->SisterID);
					if(!Sister_CAR->isCohesin)
					{
						cohesive_CARs_copy.at(i)->isCohesin=true;
						Sister_CAR->isCohesin=true;
						Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
						cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
						cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
						NbindedCohesin+=2;
					}
				}
			}
		}
	}


	//check that all the cohesive cohesins are connected with other SC
	for ( int i = 0; i < (int) tadConf.size(); ++i ) 
		if(tadConf.at(i).isCohesin)
			if(tadConf.at(i).status==tadConf.at(i).binding_site->status)
				throw std::runtime_error("cohesive cohesin illegal after loading");


}
void MCReplicPoly::LoadExtruders()
{



	if(Ntad!=2*individual_Nchain and originRate!=0)
		return;


	//load with a certain rate, if rnd is greater exit function
	double rnd = lat->rngDistrib(lat->rngEngine);
	if(rnd>loading_rate)
		return;
	
	//select a random monomer in the chain
	int t = lat->rngEngine() % Ntad;

	auto Loader_starting_monomer = &tadConf[t];
	//Loading is not permitted in fork/cohesin/end sites
	while(Loader_starting_monomer->isFork() or Loader_starting_monomer->isLeftEnd() or Loader_starting_monomer->isRightEnd() or Loader_starting_monomer->isCohesin or Loader_starting_monomer->isCAR )
	{
		t = lat->rngEngine() % Ntad;
		Loader_starting_monomer = &tadConf[t];
	}
	
	
	//find two terminal
	auto LeftAnchor = Loader_starting_monomer->neighbors[0];
	auto RightAnchor = Loader_starting_monomer->neighbors[1];


	
	if(RightAnchor->isFork() or RightAnchor->isLeftEnd() or RightAnchor->isRightEnd() or  LeftAnchor->isFork() or LeftAnchor->isLeftEnd() or LeftAnchor->isRightEnd() or RightAnchor->isCohesin or LeftAnchor->isCohesin  or RightAnchor->isCAR or LeftAnchor->isCAR )
		return;
	


	RightAnchor->isCohesin=true;
	LeftAnchor->isCohesin=true;
	RightAnchor->binding_site=LeftAnchor;
	LeftAnchor->binding_site=RightAnchor;
	//I load in the extruders vector the two
	active_extruders.push_back(LeftAnchor);

	for (int i=0 ; i < active_extruders.size() ; ++i)
		if(active_extruders.at(i)->status!=active_extruders.at(i)->binding_site->status)
			throw std::runtime_error("extruder illegal after loading");

}

void MCReplicPoly::Move_Last_Extruders()//Function to move only last monomer 
{
	

	//we set here speed of extrusion equal 1
	double extruder_speed=1;
	//always pick the last element of active extruders
	int t = (int) active_extruders.size() -1;
	auto LeftAnchor = active_extruders.at(t);
	auto RightAnchor = LeftAnchor->binding_site;


	double rnd2 = lat->rngDistrib(lat->rngEngine);
	if(rnd2>0.5) //I move left leg
	{
		int check_size= active_extruders.size();

		//move left anchor to the left and check it is not stalled
		LeftAnchor = LeftAnchor->neighbors[0];
		if(LeftAnchor->isFork() or LeftAnchor->isLeftEnd() or LeftAnchor->isRightEnd() or LeftAnchor->isCAR or LeftAnchor->isCohesin )
			return;
		


		//delete info of old extruder
		LeftAnchor->neighbors[1]->isCohesin=false;
		LeftAnchor->neighbors[1]->binding_site=nullptr;


		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
		

		//new info of the new left anchor
		LeftAnchor->isCohesin=true;
		RightAnchor->binding_site=LeftAnchor;
		LeftAnchor->binding_site=RightAnchor;
		//I load in the extruders vector the two
		active_extruders.push_back(LeftAnchor);

		if (check_size != (int) active_extruders.size() )
			throw std::runtime_error("size do not match");


		



	}
	else // I move right leg
	{

		RightAnchor = RightAnchor->neighbors[1];
		if(RightAnchor->isFork() or RightAnchor->isLeftEnd() or RightAnchor->isRightEnd() or RightAnchor->isCAR or RightAnchor->isCohesin )
			return;
		//delete info of old extruder
		RightAnchor->neighbors[0]->isCohesin=false;
		RightAnchor->neighbors[0]->binding_site=nullptr;

		//active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
		//new info of the new left anchor
		RightAnchor->isCohesin=true;
		LeftAnchor->binding_site=RightAnchor;
		RightAnchor->binding_site=LeftAnchor;

		//no nead to change active_extruders that contain only left legs
		for (int i=0 ; i < active_extruders.size() ; ++i)
			if(active_extruders.at(i)->status!=active_extruders.at(i)->binding_site->status)
				throw std::runtime_error("extruder illegal after right move");

	}
	


	//At the end I verify if I am stacking loops
	//Reload updated anchors of extruders
	LeftAnchor = active_extruders.at(t);

	RightAnchor = active_extruders.at(t)->binding_site;
	

	//my two unchors are followed by other anchors of bigger loop
	if(LeftAnchor->neighbors[0]->isCohesin and RightAnchor->neighbors[1]->isCohesin)
	{
		//check for an important special case: they two anchors are followed by cohesive cohesin
		if (LeftAnchor->neighbors[0]->binding_site->status!=LeftAnchor->neighbors[0]->binding_site->status)//check the status of the anchor neighbour on the left. its biding site must not be in other SC
			return;
		if (RightAnchor->neighbors[1]->binding_site->status!=RightAnchor->neighbors[1]->status)//check the status of the anchor neighbour on the left. its biding site must not be in other SC
			return;



		RightAnchor = RightAnchor->neighbors[1];
		LeftAnchor = LeftAnchor->neighbors[0];

		//delete info of old extruder
		LeftAnchor->neighbors[1]->isCohesin=false;
		LeftAnchor->neighbors[1]->binding_site=nullptr;
		RightAnchor->neighbors[0]->isCohesin=false;
		RightAnchor->neighbors[0]->binding_site=nullptr;

		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());

		//load a new left anchor in the the active_extruder vector (it will be a duplicate od the existing bigger loop)
		active_extruders.push_back(LeftAnchor);
		for (int i=0 ; i < active_extruders.size() ; ++i)
			if(active_extruders.at(i)->status!=active_extruders.at(i)->binding_site->status)
				throw std::runtime_error("extruder illegal after stacking");

	}

}

void MCReplicPoly::Move_Extruders()
{
	double extruder_speed=1;
	//load with a certain rate, if rnd is greater exit function
	double rnd = lat->rngDistrib(lat->rngEngine);
	if(rnd>extruder_speed)
		return;
	
	int t = lat->rngEngine() % (int) active_extruders.size();
	auto LeftAnchor = active_extruders.at(t);
	auto RightAnchor = active_extruders.at(t)->binding_site;
	double rnd2 = lat->rngDistrib(lat->rngEngine);
	if(rnd2>0.5) //I move left leg
	{
		LeftAnchor = LeftAnchor->neighbors[0];
		if(LeftAnchor->isFork() or LeftAnchor->isLeftEnd() or LeftAnchor->isRightEnd() or LeftAnchor->isCAR or LeftAnchor->isCohesin )
			return;
		
		//delete info of old extruder
		LeftAnchor->neighbors[1]->isCohesin=false;
		LeftAnchor->neighbors[1]->binding_site=nullptr;
		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
		//new info of the new left anchor
		LeftAnchor->isCohesin=true;
		RightAnchor->binding_site=LeftAnchor;
		LeftAnchor->binding_site=RightAnchor;
		//I load in the extruders vector the two
		active_extruders.push_back(LeftAnchor);
	}
	else // I move right leg
	{
		RightAnchor = RightAnchor->neighbors[1];
		if(RightAnchor->isFork() or RightAnchor->isLeftEnd() or RightAnchor->isRightEnd() or RightAnchor->isCAR or RightAnchor->isCohesin )
			return;
		//delete info of old extruder
		RightAnchor->neighbors[0]->isCohesin=false;
		RightAnchor->neighbors[0]->binding_site=nullptr;
		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
		//new info of the new left anchor
		RightAnchor->isCohesin=true;
		LeftAnchor->binding_site=RightAnchor;
		RightAnchor->binding_site=LeftAnchor;
		//no nead to change active_extruders that contain only left legs
	}
	
	//At the end I verify if I am stacking loops
	//Reload updated anchors of extruders
	LeftAnchor = active_extruders.at(t);
	RightAnchor = active_extruders.at(t)->binding_site;
	//my two unchors are followed by other anchors of bigger loop

	if(LeftAnchor->neighbors[0]->isCohesin and RightAnchor->neighbors[1]->isCohesin)
	{

		RightAnchor = RightAnchor->neighbors[1];
		

		LeftAnchor = LeftAnchor->neighbors[0];
		//delete info of old extruder
		LeftAnchor->neighbors[1]->isCohesin=false;
		LeftAnchor->neighbors[1]->binding_site=nullptr;
		RightAnchor->neighbors[0]->isCohesin=false;
		RightAnchor->neighbors[0]->binding_site=nullptr;
		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());

		//load a new left anchor in the the active_extruder vector (it will be a duplicate od the existing bigger loop)
		active_extruders.push_back(LeftAnchor);
		
	}
	if(LeftAnchor->neighbors[0]->isLeftEnd() and RightAnchor->neighbors[1]->isRightEnd())
	{
		//delete info of old extruder
		LeftAnchor->isCohesin=false;
		LeftAnchor->binding_site=nullptr;
		RightAnchor->isCohesin=false;
		RightAnchor->binding_site=nullptr;
		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
	}

}


void MCReplicPoly::unLoadExtruders()
{
	/*if( (int) active_extruders.size()>0)
	{
		std::cout <<  "CHECK valency " << std::endl;

		for ( int i = 0; i < (int) active_extruders.size(); ++i )
		{
			std::cout <<  "valency of " <<i+1<<" over "<<active_extruders.size()<<" = "<<  active_extruders.at(i)->N_loaded_extruders<< std::endl;

		}
	 }*/
	

	if((int) active_extruders.size()>0)
	{
		for ( int i = 0; i < (int) active_extruders.size(); ++i )
		{
			for ( int j = 0; j < (int) active_extruders.at(i)->N_loaded_extruders ; ++j ) //I iterate to have each cohesin with constant unloading time
			{
				double rnd = lat->rngDistrib(lat->rngEngine);
				if(rnd < unloading_rate)
				{
					if(active_extruders.at(i)->N_loaded_extruders==1)
					{
						active_extruders.at(i)->N_loaded_extruders=0;
						active_extruders.at(i)->isCohesin=false;
						active_extruders.at(i)->binding_site->isCohesin=false;
					}
					else
						--active_extruders.at(i)->N_loaded_extruders;
				}

			}
		}
		active_extruders.erase(std::remove_if(active_extruders.begin(), active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), active_extruders.end());
	}
	
}
void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->isFork()) //increase energy at fork site
	{
			if( neigh==true)
			{
			for ( int v = 0; v < 55 ; ++v )
			{
				int vo =(lattice_neigh1[v] == 0) ? tadUpdater->vo: lat->bitTable[lattice_neigh1[v]][tadUpdater->vo];
				int vn =(lattice_neigh1[v] == 0) ? tadUpdater->vn: lat->bitTable[lattice_neigh1[v]][tadUpdater->vn];
		
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];
				int vi2 = (lattice_neigh2[v] == 0) ? vn: lat->bitTable[lattice_neigh2[v]][vn];

				--lat->ReplTable[0][vi1];
				++lat->ReplTable[0][vi2];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
				
				--lat->ReplTable[0][vo];
				++lat->ReplTable[0][vn];
			}
		}
	}
	//if CAR is replicated can be turned into cohesive (only for handcuff model)
	if(Jpair>0. and tadTrial->isCAR and tadTrial->status!=0 and cohesionMode==1)
		Find_cohesive_CAR();

}

void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
	
	if(tad->isFork())
	{
		if(neigh==true)
		{
			for ( int v = 0; v < 55 ; ++v )
			{
				
				int vo =(lattice_neigh1[v] == 0) ?  tad->pos : lat->bitTable[lattice_neigh1[v]][tad->pos];
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];

				--lat->ReplTable[0][vi1];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				--lat->ReplTable[0][vo];
			}
		}
	}
	else
	{
		if(neigh==true)
		{
			for ( int v = 0; v < 55 ; ++v )
			{
				
				
				int vo =(lattice_neigh1[v] == 0) ?  tad->pos : lat->bitTable[lattice_neigh1[v]][tad->pos];
				int vi1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];

				
				++lat->ReplTable[0][vi1];


			}
		}
		else{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				++lat->ReplTable[0][vo];
			}
		}
	}
}


vtkSmartPointer<vtkPolyData> MCReplicPoly::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCHeteroPoly_looped::GetVTKData();
	
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto status = vtkSmartPointer<vtkIntArray>::New();
	auto sisterID = vtkSmartPointer<vtkIntArray>::New();
	auto cohesin = vtkSmartPointer<vtkIntArray>::New();
	auto cars = vtkSmartPointer<vtkIntArray>::New();



	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);
	
	sisterID->SetName("Sister ID");
	sisterID->SetNumberOfComponents(1);
	
	cohesin->SetName("Cohesin");
	cohesin->SetNumberOfComponents(1);
	
	cars->SetName("CAR");
	cars->SetNumberOfComponents(1);


	

	
	for ( int t = 0; t < Ntad; ++t )
	{
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
		
		forks->InsertNextValue(fork);
		status->InsertNextValue(tadConf[t].status);
		sisterID->InsertNextValue(tadConf[t].SisterID);
		//MODIFY THE OUTPUT
		cohesin->InsertNextValue(tadConf[t].isCohesin);
		cars->InsertNextValue(tadConf[t].isCAR);
	}
	
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(status);
	polyData->GetPointData()->AddArray(sisterID);
	polyData->GetPointData()->AddArray(cohesin);
	polyData->GetPointData()->AddArray(cars);

	


	
	return polyData;
}

void MCReplicPoly::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	MCHeteroPoly::SetVTKData(polyData);

	vtkDataArray* status = polyData->GetPointData()->GetArray("Replication status");
	vtkDataArray* sisterID = polyData->GetPointData()->GetArray("Sister ID");



	for ( int t = 0; t < Ntad; ++t )
	{
		tadConf[t].status = (int) status->GetComponent(t, 0);
		tadConf[t].SisterID = (int) sisterID->GetComponent(t, 0);

	}
}

void MCReplicPoly::PrintRFD()
{
	//Note: this function is only for statistics purpose and not adapted for partially replicated chromosomes

	std::ofstream outfile_rfd(outputDir+"/"+std::to_string(individual_Nchain)+"_RFD.res", std::ios_base::app | std::ios_base::out);
	RFD.at(0)=-1;
	RFD.back()=1;
	for ( int i = 0; i < individual_Nchain ; ++i )		
		outfile_rfd << RFD.at(i) << std::endl;

	outfile_rfd << -100 << std::endl;


}


void MCReplicPoly::PrintCohesins()
{
    //Note: this function is only for statistics purpose and not adapted for partially replicated chromosomes
	int Sister_chromatid1 = originRate!=0 ? Ntad/2 : Ntad; 

	std::ofstream outfile_trans(outputDir+"/"+std::to_string(Sister_chromatid1)+"_cohesion_pattern_trans.res", std::ios_base::app | std::ios_base::out);
	std::ofstream outfile_cis1(outputDir+"/"+std::to_string(Sister_chromatid1)+"_cohesion_pattern_cis1.res", std::ios_base::app | std::ios_base::out);
	std::ofstream outfile_cis2(outputDir+"/"+std::to_string(Sister_chromatid1)+"_cohesion_pattern_cis2.res", std::ios_base::app | std::ios_base::out);
	std::ofstream outfile_cars(outputDir+"/"+std::to_string(Sister_chromatid1)+"_cars.res", std::ios_base::app | std::ios_base::out);

	std::vector<int> check;
	std::cout << "PRINTING COHESINS" << std::endl;


	for ( int i = 0; i < Sister_chromatid1 ; ++i )
	{
		if(tadConf.at(i).isCAR)
		{
			outfile_cars << i << std::endl;
		}
		if(tadConf.at(i).isCohesin)
		{
			if(tadConf.at(i).binding_site->status!=tadConf.at(i).status)
			{
				//std::cout << "Cohesion: SC1 bound at " << i<< "with SC2 at "<<tadConf.at(i).binding_site->SisterID << std::endl;
				int binding_mon= (int) tadConf.at(i).binding_site->SisterID;

				outfile_trans << i << " " <<binding_mon<< std::endl;

			}
			else
			{
				//std::cout << "Looping: anchor at " << i<< " binding with anchor at "<<(int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;
				 int binding_mon= (int) std::distance(tadConf.data(), tadConf.at(i).binding_site);
				if(i<binding_mon)
					outfile_cis1 << i << " " <<binding_mon<< std::endl;

			}
			
		}
	}	
	if(Sister_chromatid1!=Ntad)
	{
		for ( int i = 0; i < Sister_chromatid1 ; ++i )
			{
			int mon_sister=tadConf.at(i).SisterID;
			if(tadConf.at(mon_sister).isCohesin)
			{
				if(tadConf.at(mon_sister).binding_site->status!=tadConf.at(mon_sister).status)
				{
					//std::cout << "Cohesion: SC2 bound at " << (int) tadConf.at(i).SisterID << "with SC1 at "<< (int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;

				}
				else
				{
					int anch1=mon_sister;
					int anch2= (int) tadConf.at((int) std::distance(tadConf.data(), tadConf.at(mon_sister).binding_site)).SisterID;
					if(anch1<anch2)
						outfile_cis2 << anch1 << " " << anch2<< std::endl;

					//std::cout << "Looping: anchor at " << i<< " binding with anchor at "<<(int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;

				}
				check.push_back((int) std::distance(tadConf.data(), tadConf.at(i).binding_site));
				
				
			}
		}	
	}
	outfile_cis1 << -1 << " " << -1 << std::endl;
	outfile_cis2 << -1 << " " << -1 << std::endl;
	outfile_trans << -1 << " " << -1 << std::endl;
	

	/*std::set<int> setOfNumbers(check.begin(), check.end());
	if (setOfNumbers.size() == check.size())
		std::cout<<"Vector has only unique values" <<std::endl;
	else
		std::cout<<"Vector is not unique" <<std::endl;*/
}
