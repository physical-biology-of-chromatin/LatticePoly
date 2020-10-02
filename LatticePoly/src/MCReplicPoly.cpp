//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>

#include "MCReplicPoly.hpp"


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);
	
	int t1 = lat->rngEngine() % Nchain;
	int t2 = lat->rngEngine() % Nchain;

	int tOrig = std::min(t1, t2);
	int tEnd = std::max(t1, t2);
			
	Replicate(tOrig, tEnd);
	
	Update();
}

void MCReplicPoly::Replicate(int tOrig, int tEnd)
{
	auto origin = tadConf.begin() + tOrig;
	auto end = tadConf.begin() + tEnd;
	
	if ( !origin->isFork() )
	{
		auto nextFork = std::find_if(origin, end+1, [](const MCTad& t){return t.isFork();});
		end = nextFork-1;
		
		int dist = (int) std::distance(origin, end);
		
		if ( dist > 1 )
		{
			ReplicateTADs(origin, end);
			ReplicateBonds(origin, end);
		}
	}
}

void MCReplicPoly::ReplicateTADs(std::vector<MCTad>::iterator origin, std::vector<MCTad>::iterator end)
{
	MCTad tadReplic;
	
	if ( origin->isLeftEnd() )
	{
		tadReplic = *origin;
		
		tadConf.push_back(tadReplic);
	}
	
	for ( auto tad = origin+1; tad != end; ++tad )
	{
		tadReplic = *tad;
		
		tadConf.push_back(tadReplic);
	}
	
	if ( end->isRightEnd() )
	{
		tadReplic = *end;
		
		tadConf.push_back(tadReplic);
	}
}

void MCReplicPoly::ReplicateBonds(std::vector<MCTad>::iterator origin, std::vector<MCTad>::iterator end)
{
	MCBond bondReplic;

	int Nreplic = (int) tadConf.size() - Ntad;
	
	int tOrig = (int) std::distance(tadConf.begin(), origin);
	int tEnd = (int) std::distance(tadConf.begin(), end);
	
	MCBond* bond = origin->isLeftEnd() ? origin->bonds[0] : origin->bonds[1];

	if ( !origin->isLeftEnd() )
	{
		bondReplic.id1 = tOrig;
		bondReplic.id2 = Ntad;
				
		bondReplic.dir = bond->dir;
	
		tadTopo.push_back(bondReplic);
		
		++bond;
	}
				
	for ( int t = 0; t < Nreplic-1; ++t )
	{
		bondReplic.id1 = t+Ntad;
		bondReplic.id2 = t+Ntad+1;
		
		bondReplic.dir = bond->dir;
		
		tadTopo.push_back(bondReplic);

		++bond;
	}
	
	if ( !end->isRightEnd() )
	{
		bondReplic.id1 = Ntad+Nreplic-1;
		bondReplic.id2 = tEnd;
		
		bondReplic.dir = bond->dir;

		tadTopo.push_back(bondReplic);
		
		++bond;
	}
}

void MCReplicPoly::Update()
{
	if ( (int) tadTopo.size() > Nbond )
	{
		for ( auto bond = tadTopo.begin()+Nbond; bond != tadTopo.end(); ++bond )
			CreateBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
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
}
