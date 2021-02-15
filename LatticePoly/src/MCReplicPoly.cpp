//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCReplicPoly.hpp"

#include <iterator>
#include <algorithm>


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);
	
	MCsteps=0;
	MCrepl=0;
	activeForks.reserve(Nchain);

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	
	Nfork = (int) activeForks.size();

	
	std::vector<int> originsvector;
	for (int i=0; i<Nchain; ++i) originsvector.push_back(i);
	std::random_shuffle ( originsvector.begin(), originsvector.end() );
	for ( int i=0 ; i < 93; i++)
	{
		origins.push_back(originsvector[i]);
	}
	/*
	origins={0,7,12,17,36,40,68,76,80,99,110,126,135,170,176,185,188,203,205,253,263,281,286,307,326,348,355,370,381,387,404,444,454,488,503,511,515,551,562,576,591,598,602,611,622,632,644,656,666,676,687,703,719,723,731,737,754,763,780,792,813,817,826,846,888,904,908,927,932,963,992,1020,1021,1042,1044,1046,1082,1084,1103,1108,1123,1143,1158,1163,1169,1175,1176,1189,1199,1202,1204,1219,1225};
	*/
	std::cout << "Origins :";
	for (std::vector<int>::iterator it=origins.begin(); it!=origins.end(); ++it)
	std::cout << ' ' << *it;


}

void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	
	// Nucleate replication bubble at random position with rate originRate (or pick position from prescribed origin locations)
	/*
	double rndOrigin = lat->rngDistrib(lat->rngEngine);
	
	if (MCsteps>Nrelax*Ninter and rndOrigin < originRate / (double) Ntad )
	{
		int t = 50;
		MCTad* tad = &tadConf[t];
		
		// Replicate chosen tad, if not yet replicated
		if ( tad->status == 0 )
			Replicate(tad);
	}
	*/
	if ( origins.size() > 0 and MCsteps> (Nrelax+1)*Ninter*Nchain )
	{


		/*Copy origins vector
		auto originsCopy =origins;

		//Pick random origin and move it (i.e. replicate it) with rate Originrate
		for ( int i=0 ; i < (int)originsCopy.size(); i++)
		{
			MCTad* origin = &tadConf[origins[i]];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < originRate/(double) 1 and origin->status==0)
			{
				Replicate(origin);
				origins.erase(origins.begin()+i);
			}
		}
	}*/
		
		int randorigin=(int) lat->rngEngine() % origins.size();
		MCTad* origin = &tadConf[origins[randorigin]];
		double rndReplic = lat->rngDistrib(lat->rngEngine);
		if ( rndReplic < ((10-Nfork/2)*origins.size()*originRate)/(double) Ntad and origin->status==0)
		{
			Replicate(origin);
			origins.erase(origins.begin()+randorigin);
		}
	}
		
	
	if ( Nfork > 0)
	{
		/*Copy activeFork vector
		auto activeForksCopy =activeForks;
		// Pick random fork and move it (i.e. replicate it) with rate replicRate
		for ( int i=0 ; i < (int)activeForksCopy.size(); i++)
		{
			MCTad* fork = activeForks[i];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < replicRate/(double) 1 and fork->isFork())
			{
				Replicate(fork);
			}
		}
	}*/


		
		int randfork=(int) lat->rngEngine() % activeForks.size();
		MCTad* fork = activeForks[randfork];
		double rndReplic = lat->rngDistrib(lat->rngEngine);
		if ( rndReplic < activeForks.size()*replicRate/(double) Ntad and fork->isFork())
		{
			Replicate(fork);
		}
	}
	
	MCsteps+=1;
	
}

void MCReplicPoly::Replicate(MCTad* tad)
{
	if (tad->isRightEnd() || tad->isLeftEnd())
		return;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
		
	double rnd = lat->rngDistrib(lat->rngEngine);
	
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
				UpdateReplTable(nb1);
			}
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
			{
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);
			}
		}
	}
	
	else
	{
		// Replicating left fork means displacing it to its left neighbor, so we need to check if it's already a fork or chain end
		if ( tad->isLeftFork() )
		{
			if ( nb1->isLeftFork() || nb1->isRightEnd() )
				// Probably should never happen, do nothing
				return;
			
			if ( nb1->isRightFork() || nb1->isLeftEnd() )
			{
				// Merge forks/replicate extremities at half the normal rate
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb1);
			}
		}
		
		// Same for right forks
		else if ( tad->isRightFork() )
		{
			if ( nb2->isRightFork() || nb2->isLeftEnd() )
				return;
			
			if  ( nb2->isLeftFork() || nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb2);
			}
		}
		
		// Delete old forks
		auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
		activeForks.erase(fork);
		
		if ( nb1->isFork() || nb2->isFork() )
		{
			auto fork2 = std::find(activeForks.begin(), activeForks.end(), nb1->isFork() ? nb1 : nb2);
			activeForks.erase(fork2);
		}
	}
	
	ReplicateTADs(tad);
	ReplicateBonds(tad);

	Update();
	UpdateReplTable(tad);

	
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


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);

	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		tadConf.push_back(tadReplic);
		nb2->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb2);


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
		
		CreateBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			CreateBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		CreateBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			CreateBond(*bond3);
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
			CreateBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	// Update tads
	if ( (int) tadConf.size() > Ntad )
	{
		for ( auto tad = tadConf.begin()+Ntad; tad != tadConf.end(); ++tad )
		{
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
	
	Nfork = (int) activeForks.size();
}

double MCReplicPoly::GetEffectiveEnergy() const
{
	if ( Jf > 0.  )
	{
		
		if ( tadTrial->isFork()){
			return 	MCHeteroPoly::GetEffectiveEnergy() +Jf * (ReplTable[0][tadUpdater->vo]-ReplTable[0][tadUpdater->vn]);
		}
	}
	if ( Jf > 0.  )
	{
		if (tadTrial->status == 1)
		{
			double Jbott1=0.0;
			double Jbott2=0.0;

			int SistPos = tadConf.at(tadTrial->SisterID).pos;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				if(SistPos==vi2 and tadTrial->SisterID%1==0){
					Jbott2=Jpair;
				}
				if(SistPos==vi1 and tadTrial->SisterID%1==0){
					Jbott1=Jpair;
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		if (tadTrial->status == -1)
		{
			double Jbott1=0.0;
			double Jbott2=0.0;

			int SistPos = tadConf.at(tadTrial->SisterID).pos;
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				if(SistPos==vi2 and tadConf.at(tadTrial->SisterID).SisterID%10==0){
					Jbott2=Jpair;
				}
				if(SistPos==vi1 and tadConf.at(tadTrial->SisterID).SisterID%10==0){
					Jbott1=Jpair;
				}
			}


			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
	}
	return MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->isFork())
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			--ReplTable[0][vi1];
			++ReplTable[0][vi2];
		}
	}
	if ( tadTrial->status == -1)
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];

			--ReplTable[1][vi1];
			++ReplTable[1][vi2];
			
		}
	}
	if ( tadTrial->status == 1)
	{
	
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			
			--ReplTable[2][vi1];
			++ReplTable[2][vi2];
		}
	}
}

void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
	if(tad->status == -1)
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[1][vi];
		}
	}
	if(tad->status == 1)
	{
	for ( int v = 0; v < 13; ++v )
		{

			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[2][vi];
		}
	}
	
	if(tad->isRightEnd() || tad->isLeftEnd())
	{
	for ( int v = 0; v < 13; ++v )
	{
		int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
	
		--ReplTable[0][vi];
		}
	}
	else{
		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[0][vi];
		}
	}
}
std::vector<double3> MCReplicPoly::GetPBCConf()
{
	std::vector<MCTad*> leftEnds;
	std::vector<MCTad*> builtTads;
	
	std::vector<double3> conf(Ntad);
	
	builtTads.reserve(Ntad);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			conf[t][i] = lat->xyzTable[i][tadConf[t].pos];
		
		if ( tadConf[t].isLeftEnd() )
			leftEnds.push_back(&tadConf[t]);
	}
	
	// Grow chains recursively, starting from their respective left extremities
	auto leftEnd = leftEnds.begin();
	
	while ( (int) builtTads.size() < Ntad )
	{
		MCTad *tad1, *tad2;
		tad1 = *leftEnd;
		
		bool builtTad1 = (std::find(builtTads.begin(), builtTads.end(), tad1) != builtTads.end());

		if ( !builtTad1 )
		{
			builtTads.push_back(tad1);

			// Traverse main branch
			while ( (tad2 = tad1->neighbors[1]) )
			{
				bool builtTad2 = (std::find(builtTads.begin(), builtTads.end(), tad2) != builtTads.end());

				if ( !builtTad2 )
				{
					BuildPBCPair(builtTads, conf, tad1, tad2);
					
					// Traverse side branches
					if ( tad2->isFork() )
					{
						MCTad *tad3, *tad4;
						tad3 = tad2->neighbors[2];
						
						BuildPBCPair(builtTads, conf, tad2, tad3);
					
						while ( (tad4 = (tad2->isLeftFork() ? tad3->neighbors[1] : tad3->neighbors[0])) )
						{
							BuildPBCPair(builtTads, conf, tad3, tad4);
						
							if ( tad4->isFork() )
								break;
							
							tad3 = tad4;
						}
					}
				}
				
				tad1 = tad2;
			}
		}
		
		++leftEnd;
	}
	
	int chainNum = Ntad / Nchain;
	int chainLength = (chainNum == 1) ? Ntad : Nchain;
	
	std::vector<double3> centers(chainNum);
	
	for ( int c = 0; c < chainNum; ++c )
	{
		auto end1 = conf.begin() + c*chainLength;
		auto end2 = conf.begin() + (c+1)*chainLength;

		centers[c] = GetPBCCenterMass(end1, end2);
	}
	
	centerMass = {0., 0., 0.};
	
	for ( int c = 0; c < chainNum; ++c )
	{
		for ( int i = 0; i < 3; ++i )
			centerMass[i] += centers[c][i] / ((double) chainNum);
	}
	
	return conf;
}

void MCReplicPoly::BuildPBCPair(std::vector<MCTad*>& builtTads, std::vector<double3>& conf, MCTad* tad1, MCTad* tad2)
{
	int id1 = (int) std::distance(tadConf.data(), tad1);
	int id2 = (int) std::distance(tadConf.data(), tad2);
	
	FixPBCPair(conf, id1, id2);
	
	builtTads.push_back(tad2);
}
