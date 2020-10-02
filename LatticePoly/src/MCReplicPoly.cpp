//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCReplicPoly.hpp"


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);
	
	int t1 = lat->rngEngine() % Ntad;
	int t2 = lat->rngEngine() % Ntad;

	int tmin = std::min(t1, t2);
	int tmax = std::max(t1, t2);
	
	Replicate(tmin, tmax);
	
	Update();
}

void MCReplicPoly::Replicate(int tmin, int tmax)
{
	MCTad* origin = &tadConf[tmin];
	MCTad* end = &tadConf[tmax];
	
	for ( auto tad = origin; tad != end + 1; ++tad )
	{
		if ( tad->isFork() )
		{
			end = tad-1;
			break;
		}
		
		MCTad tadRepl(*tad);
		
		tadConf.push_back(tadRepl);
	}
		
	long Nbranch = end-origin;

	if ( Nbranch < 2 )
		tadConf.resize(Ntad);
	
	else
	{
		MCLink bondRepl;

		if ( !origin->isLeftEnd() )
		{
			bondRepl.id1 = tmin;
			bondRepl.id2 = Ntad;
		
			tadTopo.push_back(bondRepl);
		}
		
		MCLink* bond = origin->isLeftEnd() ? origin->bonds[0] : origin->bonds[1];
		
		for ( int t = Ntad + 1; t <= Ntad + (int) Nbranch; ++t )
		{
			bondRepl.id1 = t-1;
			bondRepl.id2 = t;
			
			bondRepl.dir = bond->dir;
			
			tadTopo.push_back(bondRepl);

			++bond;
		}
		
		if ( !end->isRightEnd() )
		{
			bondRepl.id1 = Ntad + (int) Nbranch;
			bondRepl.id2 = tmin + (int) Nbranch;
		
			tadTopo.push_back(bondRepl);
		}
	}
}

void MCReplicPoly::Update()
{
	if ( (int) tadTopo.size() > Nbond )
	{
		for ( auto bond = tadTopo.begin() + Nbond; bond != tadTopo.end(); ++bond )
			CreateBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	if ( (int) tadConf.size() > Ntad )
	{
		for ( auto tad = tadConf.begin() + Ntad; tad != tadConf.end(); ++tad )
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
