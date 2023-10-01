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

void MCReplicPoly::Init(int Ninit)
{
	MCLivingPoly::Init(Ninit);
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
	
	for ( int vi = 0; vi < Ntot; ++vi )
		ReplTable[0][vi] = 0;

	
	activeForks.reserve(Nchain);
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
				
				double d1;
				
				if ( ss >> d1 )
				{

				ChIP.push_back(d1);
					
				}
			}
			if (Nchain != (int) ChIP.size() )
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
				
				if ( ss >> d1 )
				{
					PODLS.push_back(d1);

				}
			}
			if (Nchain != (int) PODLS.size() )
			throw std::runtime_error("Nchain and PODLS size do not match");
			
			PODLSfile.close();
		}
	}
	
	
	if(StartFromPODLS==true)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<> d(PODLS.begin(), PODLS.end());
		origins={};
		for(int n=0; n<int(Nchain/5); ++n)
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
				if (d1 >Nchain )
					throw std::runtime_error("Nchain and origin ID do not match");
				origins.push_back(d1-1);
				
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

		active_cars={};
		while((int) active_cars.size() < n_barriers )
		{
			int car=d(gen);
			if(std::find(active_cars.begin(),active_cars.end(),car) == active_cars.end())
			   active_cars.push_back(car);

			
		}

		
		
		for (int i=0 ; i < (int) active_cars.size();++i)
			tadConf[active_cars[i]].isCAR=true;
		
		
		std::cout << "Extruder before " << N_extruders <<std::endl;
		
		N_extruders = N_extruders*active_cars.size();
		std::cout << "Extruder after " << N_extruders <<std::endl;

	}
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
					if(rndReplic < Ntot*originRate/12  and origin->status==0 and !origin->isFork())
					{
						
						//a particle can activate only one origin
						if(std::find(Spin_pos_toDelete.begin(),Spin_pos_toDelete.end(),pos) == Spin_pos_toDelete.end())
						{
							
							Replicate(origin);
							//check if replication occurs
							if(origin->status!=0)
							Spin_pos_toDelete.push_back(pos);
						}
					}
				}
			}
		}
		
	}
}
void  MCReplicPoly::OriginMove_implicit()
{
	if ( (int) activeOrigins.size() > 0 )
	{
		auto originsCopy =activeOrigins;
		std::shuffle (originsCopy.begin(), originsCopy.end(), lat->rngEngine);
		
		
		for ( int i=0 ; i < (int)originsCopy.size(); i++) //for every element in indexes
		{
			
			MCTad* origin = originsCopy[i]; //select origin taf
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			
			int Nocc = activeForks.size() % 2 == 0 ? int(activeForks.size()) : int(activeForks.size())+ 1;
			
			// -1 since origin firing implicate 2 new monomer in the system
			if ( rndReplic < double(2*Ndf- Nocc) * originRate and origin->status==0 and  Ntad < Nchain + int(stop_replication))//Ntad < Nchain -2 + int(Nchain * stop_replication))
			{
				
				Replicate(origin);
				
			}
		}
	}
}
	
	/*
	if(Ntad>=int(.95*Nchain+Nchain))
	{
		
		std::ostringstream streamObj;
		streamObj << originRate;
		std::string strObj = streamObj.str();

		std::ofstream outfile(outputDir+"repltime"+std::to_string(Ndf)+ "_" + strObj+".res", std::ios_base::app | std::ios_base::out);

		outfile << MCsteps << std::endl;


		
		
		exit(0);
		
	}*/

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
			if ( fork->status==0 and rndReplic < replicRate and Ntad < Nchain - 1  + int(stop_replication))//Nchain + int(Nchain*stop_replication) )
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
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb2;
				
				if(nb1->isCAR)
					TurnCohesive(nb1);

			}
		
			else
			{

				activeForks.push_back(nb1); // STANDARD CASE
				UpdateReplTable(nb1);//increase energy around new fork
				
				nb1->binding_site=tad->binding_site;
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb1;

				if(nb2->isCAR)
					TurnCohesive(nb2);

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
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb2;
				
				if(nb2->isCAR)
					TurnCohesive(nb2);

				
			}
		
			else
			{
				
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
				nb2->binding_site=tad->binding_site;
				if(tad->binding_site->isFork())
					tad->binding_site->binding_site=nb2;
				
				if(nb1->isCAR)
					TurnCohesive(nb1);

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
	
	if(Ntad==2*Nchain-2)
		Update_rcms_before_separation();
		

	
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
			if(nb1->isLeftEnd() and tadConf.at(Nchain-1).status!=0)
				Spin_pos_toCreate.push_back(tad->pos);
			else if(nb1->isRightFork())
				Spin_pos_toCreate.push_back(tad->pos);
		}
		
		if(nb1->isCAR)
		{

			tadConf.back().isCAR=true;
			//MODIFY
			//tadConf.back().isChoesin=true;
			//nb1->isChoesin=true;
			//nb1->binding_site = &tadConf.back();
			//tadConf.back().binding_site=nb1;
		}


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);
	
	
	if(tad->isCAR)
	{

		tadConf.back().isCAR=true;
		//MODIFY
		//tadConf.back().isChoesin=true;
		//tad->isChoesin=true;
		//tad->choesin_binding_site = &tadConf.back();
		//tadConf.back().choesin_binding_site=tad;
	}

	
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
			//MODIFY

			//tadConf.back().isChoesin=true;
			//nb2->isCohesin=true;
			//nb2->binding_site = &tadConf.back();
			//tadConf.back().binding_site=nb2;
		}

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
				tad->isCentromere=true;
			
			if ( tad->type == 1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable[vi];
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
	/*NbindedForks = 0;
	for (int i=0; i < (int) activeForks.size();++i)
	{
		if (activeForks.at(i)->binding_site->isFork())
		{
			int pos_binded=activeForks.at(i)->binding_site->pos;
			if ( Jf_sister > 0.  and neigh==1)
			{
				
				for ( int v = 0; v < 55 ; ++v )
				{
					
					int vo =(lattice_neigh1[v] == 0) ? activeForks.at(i)->pos : lat->bitTable[lattice_neigh1[v]][activeForks.at(i)->pos];
					int v1 = (lattice_neigh2[v] == 0) ? vo: lat->bitTable[lattice_neigh2[v]][vo];
					if(v1==pos_binded)
					{
						++NbindedForks;
						break;
					}
				}
			}
			
			if (Jf_sister > 0. and neigh==0 )
			{

				for ( int v = 0; v < 13 ; ++v )
				{
					
					
					int vo =(v == 0) ?  activeForks.at(i)->pos : lat->bitTable[v][activeForks.at(i)->pos];
					if(vo==pos_binded)
					{
						++NbindedForks;
						break;
					}
				}
			}
		}

	}*/
	//check how many forks are binded to their sister
	NbindedForks = 0;
	 for (int i=0; i < (int) activeForks.size();++i)
		 if (activeForks.at(i)->binding_site->isFork())
			 ++NbindedForks;


	 
	 
	
	if(cohesionMode!=1)
		Find_cohesive_CAR();
}

double MCReplicPoly::GetEffectiveEnergy() //chiedere Maxime
{
	if (tadTrial->isFork() )
	{
		
		double Etot = 0.;

		if ( Jf > 0.  )
			Etot=Etot+Jf*(ReplTable[0][tadUpdater->vo]-ReplTable[0][tadUpdater->vn]);


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

	//if (tadTrial->isCohesin  or (!tadTrial->isCohesin and tadTrial->isCAR and tadTrial->binding_site== &tadConf.at(tadTrial->SisterID))) //need to ask maxime and to test it
	if (tadTrial->isCohesin) //need to ask maxime and to test it

	{
		//if((!tadTrial->isCohesin and tadTrial->isCAR and tadTrial->binding_site == &tadConf.at(tadTrial->SisterID)))
			//std::cout <<  "CAR NOT COHESIN"<< std::endl;

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

	//std::cout <<  "Turn COHESIVE"<< std::endl;

	
	if( std::find(cohesive_CARs.begin(),cohesive_CARs.end(),tad) == cohesive_CARs.end())
	{

		double rnd = lat->rngDistrib(lat->rngEngine);
		//int original_total_activated_cars=total_activated_cars;
		double activation_rate = ForkTableMode==0? keco1 : keco1*ReplTable[0][tad->pos];
		if(rnd<activation_rate)
		{
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
		
		
		
		/*if(cohesionMode!=2) //activate both ends only when it is not homologous
		{
			if(rnd<activation_rate)
			{

				cohesive_CARs.push_back(&tadConf.at(tad->SisterID));
				++total_activated_cars;
			}
		}
		if((cohesionMode==0 or cohesionMode==3 or cohesionMode==4)  and (total_activated_cars-original_total_activated_cars)==1) //in case of non homologous cohesin stack mantain transitory binding between two Sc if I activate only one CAR
		{
			tad->binding_site=&tadConf.at(tad->SisterID);
			tadConf.at(tad->SisterID).binding_site=tad;

		}*/
		
	}
}

void MCReplicPoly::Find_cohesive_CAR()
{
	//std::cout <<  "FIND COHESIVE"<< std::endl;

	if(cohesionMode!=1  and cohesionMode!=2)
	{
		if(cohesive_CARs.size()>1 )
		{

			auto cohesive_CARs_copy=cohesive_CARs;
			std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);


			for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )//loop over all cohesive CARs
			{

				
				if(!cohesive_CARs_copy.at(i)->isCohesin)
				{
					//std::cout <<  "CAR STATUS  " << cohesive_CARs_copy.at(i)->status<<std::endl;

					auto Sister_CAR=&tadConf.at( cohesive_CARs_copy.at(i)->SisterID);
					if(Sister_CAR->isCohesin) //if the homologous is already a cohesive coesin avoid crossing
						return;
					
					auto tad_shifter= Sister_CAR;
					//Make a symmetrical binding
					double rnd_symm = lat->rngDistrib(lat->rngEngine);

					if(cohesionMode!=4)
						rnd_symm = -1;
					if(rnd_symm > 0.0) //homologous binding
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
								//std::cout <<  "shifter status - s " << tad_shifter->status<<std::endl;
								tad_shifter=tad_shifter->neighbors[0];
								//std::cout <<  "shifter is fork " << tad_shifter->isFork()<<std::endl;


								if( tad_shifter->isCohesin and !tad_shifter->isFork())
								{
									
									
									tad_shifter=tad_shifter->neighbors[1];
									
									//cohesive_CARs_copy.at(i)->binding_site->binding_site=nullptr;
									
									cohesive_CARs_copy.at(i)->isCohesin=true;
									Sister_CAR->isCohesin=true;
									Sister_CAR->binding_site=cohesive_CARs_copy.at(i);
									cohesive_CARs_copy.at(i)->binding_site=Sister_CAR;
									//auto del = std::find(cohesive_CARs.begin(), cohesive_CARs.end(), cohesive_CARs_copy.at(i));
									cohesive_CARs.erase(std::remove_if(cohesive_CARs.begin(), cohesive_CARs.end(), [](const MCTad* tad){return tad->isCohesin;}), cohesive_CARs.end());
									//cohesive_CARs.erase(del);
									NbindedCohesin+=2;
									//std::cout <<  "FOUND - cohesin " << std::endl;
									//std::cout <<  "anchor found with status  " << Sister_CAR->status<< " and binding of status " << Sister_CAR->binding_site->status<< std::endl;
									
									PrintCohesins();
									
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
									
									
									//std::cout <<  "FOUND - " << std::endl;
									//std::cout <<  "anchor found with status  " << Sister_CAR->status<< " and binding of status " << Sister_CAR->binding_site->status<< std::endl;
									
									PrintCohesins();
									
									return;
								}
							}
						}
						else
						{
							while( !tad_shifter->isFork() and !tad_shifter->isRightEnd() and !tad_shifter->isLeftEnd())
							{
								//std::cout <<  "shifter status +  s " << tad_shifter->status<<std::endl;
								tad_shifter=tad_shifter->neighbors[1];

								if(tad_shifter->isCohesin and !tad_shifter->isFork())
								{
									tad_shifter=tad_shifter->neighbors[0];
									std::cout <<  "shifter status cohesin   " << tad_shifter->status<<std::endl;
									
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
									//std::cout <<  "FOUND + cohesin " << std::endl;
									//std::cout <<  "anchor found with status  " << Sister_CAR->status<< " and binding of status " << Sister_CAR->binding_site->status<< std::endl;
									
									PrintCohesins();
									
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
									
									PrintCohesins();
									
									return;
								}
							}
						}
					}
				}
			}
		}
	}

	if(cohesionMode==1)
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
	if(cohesionMode==2)
	{

		if(cohesive_CARs.size()>1 )
		{
			

				
			auto cohesive_CARs_copy=cohesive_CARs;
			std::shuffle (cohesive_CARs_copy.begin(), cohesive_CARs_copy.end(), lat->rngEngine);
			for ( int i = 0; i < (int) cohesive_CARs_copy.size(); ++i )
			{
				if(!cohesive_CARs_copy.at(i)->isCohesin and cohesive_CARs_copy.at(i)->status!=0)
				{
					auto Sister_CAR=&tadConf.at( cohesive_CARs_copy.at(i)->SisterID);
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
}

void MCReplicPoly::LoadExtruders()
{
	//load with a certain rate, if rnd is greater exit function
	double rnd = lat->rngDistrib(lat->rngEngine);
	if(rnd>loading_rate)
		return;
	
	//select a random monomer in the chain
	int t = lat->rngEngine() % Ntad;
	//std::cout <<  "loading  " << t << std::endl;

	auto Loader_starting_monomer = &tadConf[t];
	//Loading is not permitted in fork/cohesin/end sites
	while(Loader_starting_monomer->isFork() or Loader_starting_monomer->isLeftEnd() or Loader_starting_monomer->isRightEnd() or Loader_starting_monomer->isCohesin or Loader_starting_monomer->isCAR )
	{
		t = lat->rngEngine() % Ntad;
		Loader_starting_monomer = &tadConf[t];
	}
	
	//if I land in a CAR which is not a cohesin I go either right or left of it
	/*if(Loader_starting_monomer->isCAR)
	{
		rnd = lat->rngDistrib(lat->rngEngine);
		if(rnd > 0.5)
			Loader_starting_monomer=Loader_starting_monomer->neighbors[0];
		else
			Loader_starting_monomer=Loader_starting_monomer->neighbors[1];
		
		if(Loader_starting_monomer->isFork() or Loader_starting_monomer->isLeftEnd() or Loader_starting_monomer->isRightEnd()) //I need to check if this move created inconsistency
			return;
	}*/

	
	
	//find two terminal
	auto LeftAnchor = Loader_starting_monomer->neighbors[0];
	auto RightAnchor = Loader_starting_monomer->neighbors[1];
	
	//if after the movements I'm at a fork or end I return
	if(RightAnchor->isFork() or RightAnchor->isLeftEnd() or RightAnchor->isRightEnd() or  LeftAnchor->isFork() or LeftAnchor->isLeftEnd() or LeftAnchor->isRightEnd())
		return;
	
	//move left anchor to the left
	while(!LeftAnchor->isCohesin)
	{
		
		//if after the movements I'm at a fork or end I return
		if(LeftAnchor->isFork() or LeftAnchor->isLeftEnd() or LeftAnchor->isRightEnd() )
			return;
		if(LeftAnchor->isCAR ) //if I find a CAR check for permeability (I know it cannot be a cohesin)
		{
			rnd = lat->rngDistrib(lat->rngEngine);
			if(rnd<permeability)
				break;
		}
		LeftAnchor=LeftAnchor->neighbors[0];

	}
	//move right anchor to the right
	while(!RightAnchor->isCohesin)
	{
		if(RightAnchor->isFork() or RightAnchor->isLeftEnd() or RightAnchor->isRightEnd() )
		{
			std::cout <<  "  END  " << std::endl;
			return;
		}
		if(RightAnchor->isCAR )//if I find a CAR check for permeability
		{
			rnd = lat->rngDistrib(lat->rngEngine);
			if(rnd<permeability)
				break;
		}
		RightAnchor=RightAnchor->neighbors[1];
	}

	
	//if both are cohesin
	if(RightAnchor->isCohesin and LeftAnchor->isCohesin)
	{
		//if are cohesive cohesin move back of one step
		if(RightAnchor->binding_site->status != RightAnchor->status)
			RightAnchor=RightAnchor->neighbors[0];
		if(LeftAnchor->binding_site->status != LeftAnchor->status)
			LeftAnchor=LeftAnchor->neighbors[1];
		
		//if are still both cohesin meand that are occupied by extruders so I check if the two anchor correspond to existing extruder
		if(RightAnchor->isCohesin and LeftAnchor->isCohesin)
		{
			if(RightAnchor->binding_site==LeftAnchor)//existing extruder, update  N_loaded_extruders without new binding
			{
				++LeftAnchor->N_loaded_extruders;
				std::cout <<  "  REINFORCEMENT OF EXISTING MONOMER " << std::endl;
				int active_extruders_count=0;
				for (int i=0 ; i < (int) active_extruders.size() ; ++i)
					active_extruders_count=active_extruders_count+ (int) active_extruders.at(i)->N_loaded_extruders;
				std::cout <<  "  NEW CONNECTION MADE " <<active_extruders_count<< std::endl;

				return;

			}
			else //I am between two extruders so I shift one step back
			{
				RightAnchor=RightAnchor->neighbors[0];
				LeftAnchor=LeftAnchor->neighbors[1];
			}
			

		}
	}
	
	//if only one is cohesin I move back one step, valid for all cases
	if(RightAnchor->isCohesin)
		RightAnchor=RightAnchor->neighbors[0];
	if(LeftAnchor->isCohesin)
		LeftAnchor=LeftAnchor->neighbors[1];
	
	//both should not be cohesin now, I make the connection
	if(RightAnchor->isCohesin or LeftAnchor->isCohesin)
		std::cout <<  "ERROR " << std::endl;

	
	RightAnchor->isCohesin=true;
	LeftAnchor->isCohesin=true;
	RightAnchor->binding_site=LeftAnchor;
	LeftAnchor->binding_site=RightAnchor;
	++LeftAnchor->N_loaded_extruders;
	//I load in the extruders vector the two
	active_extruders.push_back(LeftAnchor);
	int active_extruders_count=0;
	for (int i=0 ; i < (int) active_extruders.size() ; ++i)
		active_extruders_count=active_extruders_count+ (int) active_extruders.at(i)->N_loaded_extruders;
	std::cout <<  "  NEW CONNECTION MADE " <<active_extruders_count<< std::endl;


	
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

				--ReplTable[0][vi1];
				++ReplTable[0][vi2];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
				
				--ReplTable[0][vo];
				++ReplTable[0][vn];
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

				--ReplTable[0][vi1];
			}
		}
		else
		{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				--ReplTable[0][vo];
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

				
				++ReplTable[0][vi1];


			}
		}
		else{
			for ( int v = 0; v < 13 ; ++v )
			{
				int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
				++ReplTable[0][vo];
			}
		}
	}
}


vtkSmartPointer<vtkPolyData> MCReplicPoly::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCHeteroPoly::GetVTKData();
	
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
