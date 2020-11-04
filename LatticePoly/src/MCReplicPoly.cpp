//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCReplicPoly.hpp"

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
		
	double rndOrigin = lat->rngDistrib(lat->rngEngine);
	
	// Nucleate replication bubble at random position with rate originRate (or pick position from prescribed origin locations)
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
		// Pick random fork and move it (along the chain) with rate replicRate
		double rndReplic = lat->rngDistrib(lat->rngEngine);

		if ( rndReplic < replicRate * Nfork / (double) Ntad )
		{
			int f = lat->rngEngine() % Nfork;
			MCTad* fork = activeForks[f];

			// Moving an existing fork implies replicating it
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
			// Replicate exremities at half the normal rate
			if ( !nb1->isLeftEnd() )
				activeForks.push_back(nb1);
			else if ( rnd < 0.5 )
				return;
			
			if ( !nb2->isRightEnd() )
				activeForks.push_back(nb2);
			else if ( rnd < 0.5 )
				return;
		}
	}
	
	else
	{
		// Replicating left fork means displacing it to its left neighbor, so we need to check if it's already a fork
		if ( tad->isLeftFork() )
		{
			if ( nb1->isRightFork() )
				// Coalesce!
				return;
		
			else if ( nb1->isLeftFork() )
				// Probably should never happen, do nothing
				return;
			
			else
			{
				auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
				activeForks.erase(fork);
				
				if ( !nb1->isLeftEnd() )
					activeForks.push_back(nb1);
			}
		}
		
		// Same for right fork
		else if ( tad->isRightFork() )
		{
			if  ( nb2->isLeftFork() )
				// Coalesce!
				return;
		
			else if ( nb2->isRightFork() )
				// Probably should never happen, do nothing
				return;
			
			else
			{
				auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
				activeForks.erase(fork);
				
				if ( !nb2->isRightEnd() )
					activeForks.push_back(nb2);
			}
		}
	}
	
	ReplicateTAD(tad);
	ReplicateBonds(tad);

	Update();
}

void MCReplicPoly::ReplicateTAD(MCTad* tad)
{
	MCTad tadReplic;
	
	// Replicate left end (if it is a neighbor) to avoid formation of end loops
	if ( tad->neighbors[0]->isLeftEnd() )
	{
		tadReplic = *tad->neighbors[0];
		
		tadConf.push_back(tadReplic);
	}
	
	// Replicate TAD
	tadReplic = *tad;
		
	tadConf.push_back(tadReplic);
	
	// Same for right end
	if ( tad->neighbors[1]->isRightEnd() )
	{
		tadReplic = *tad->neighbors[1];
		
		tadConf.push_back(tadReplic);
	}
}

void MCReplicPoly::ReplicateBonds(MCTad* tad)
{
	MCBond bondReplic1;
	MCBond bondReplic2;

	// If replicated tad is a fork, replicated bonds should link new fork to both the replicated and main branches
	MCBond* bond1 = tad->isRightFork() ? tad->bonds[2] : tad->bonds[0];
	MCBond* bond2 = tad->isLeftFork() ? tad->bonds[2] : tad->bonds[1];
	
	// Create/modify bonds between relevant neighbors and replicated tad
	bondReplic1.id1 = tad->neighbors[0]->isLeftEnd() ? Ntad : bond1->id1;
	bondReplic1.id2 = tad->neighbors[0]->isLeftEnd() ? Ntad+1 : Ntad;
				
	bondReplic1.dir = bond1->dir;

	bondReplic2.id1 = tad->neighbors[0]->isLeftEnd() ? Ntad+1 : Ntad;
	bondReplic2.id2 = tad->neighbors[1]->isRightEnd() ? Ntad+1 : bond2->id2;
		
	bondReplic2.dir = bond2->dir;
	
	// For right fork, update bond1 to link left (replicated) neighbor to new tad
	if ( tad->isRightFork() )
	{
		bond1->id2 = bondReplic1.id2;
		
		CreateBond(*bond1);
		UnsetFork(tad);
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for right neighbor of left fork
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		CreateBond(*bond2);
		UnsetFork(tad);
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
