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
	
	activeForks.reserve(Nchain);

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	
	Nfork = (int) activeForks.size();
	
	// Deterministic origin locations can also be set here (or read from file) in a new array
}

void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);

	// Nucleate replication bubble at random position with rate originRate (or pick position from prescribed origin locations)
	double rndOrigin = lat->rngDistrib(lat->rngEngine);
	
	if ( rndOrigin < originRate / (double) Ntad )
	{
		int t = lat->rngEngine() % (Nchain-2) + 1;
		MCTad* tad = &tadConf[t];
		
		// Replicate chosen tad, if not yet replicated
		if ( tad->status == 0 )
			Replicate(tad);
	}
	
	if ( Nfork > 0 )
	{
		// Pick random fork and move it (i.e. replicate it) with rate replicRate
		double rndReplic = lat->rngDistrib(lat->rngEngine);

		if ( rndReplic < replicRate * Nfork / (double) Ntad )
		{
			int f = lat->rngEngine() % Nfork;
			MCTad* fork = activeForks[f];

			Replicate(fork);
		}
	}
}

void MCReplicPoly::Replicate(MCTad* tad)
{
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
				activeForks.push_back(nb1);
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
				activeForks.push_back(nb2);
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
				activeForks.push_back(nb1);
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
				activeForks.push_back(nb2);
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
	}
	
	// Replicate TAD
	tadReplic = *tad;
		
	tadConf.push_back(tadReplic);
	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		
		tadConf.push_back(tadReplic);
	}
}

void MCReplicPoly::ReplicateBonds(MCTad* tad)
{
	MCBond bondReplic1;
	MCBond bondReplic2;

	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
	
	// If replicated tad is a fork, replicated bonds should link new fork to both the replicated and main branches
	MCBond* bond1 = tad->isRightFork() ? tad->bonds[2] : tad->bonds[0];
	MCBond* bond2 = tad->isLeftFork() ? tad->bonds[2] : tad->bonds[1];
	
	// Create/modify bonds between relevant neighbors and replicated tad
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
	
	MCTad *tad1, *tad2, *tad3, *tad4;

	while ( (int) builtTads.size() < Ntad )
	{
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
	
	std::vector<double3> centreMassPBC(chainNum);
	std::fill(centreMassPBC.begin(), centreMassPBC.end(), (double3) {0., 0., 0.});
	
	for ( int c = 0; c < chainNum; ++c )
	{
		for ( int i = 0; i < 3; ++i )
		{
			for ( int t = c*chainLength; t < (c+1)*chainLength; ++t )
				centreMassPBC[c][i] += conf[t][i] / ((double) chainLength);
				
			double deltaCentreMass = centreMassPBC[c][i] - centreMass[i];
				
			while ( std::abs(deltaCentreMass) > L/2. )
			{
				double pbcShift = std::copysign(L, deltaCentreMass);
					
				for ( int t = c*chainLength; t < (c+1)*chainLength; ++t )
					conf[t][i] -= pbcShift;
					
				deltaCentreMass -= pbcShift;
				centreMassPBC[c][i] -= pbcShift;
			}
		}
	}
	
	centreMass = {0., 0., 0.};
	
	for ( int c = 0; c < chainNum; ++c )
	{
		for ( int i = 0; i < 3; ++i )
			centreMass[i] += centreMassPBC[c][i] / ((double) chainNum);
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
