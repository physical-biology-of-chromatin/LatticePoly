//
//  MCHeteroPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCHeteroPoly.hpp"


MCHeteroPoly::MCHeteroPoly(MCLattice* _lat): MCPoly(_lat) {};

void MCHeteroPoly::Init(std::mt19937_64& rngEngine)
{
	MCPoly::Init(rngEngine);
		
	for ( int i = 0; i < Ntot; i++ ) 
	{
		tadTable[i] = 0;
	}
	
	for ( int i = 0; i < Ndom; i++ )
	{
		int idx = rngEngine() % (Nchain-Nloc+1);

		for ( int j = 0; j < Nloc; j++ )
		{
			tadType[idx+j] = 1;
			tadTable[tadConf[idx+j]] = 1;
		}
	}
}

void MCHeteroPoly::TrialMoveSpinTAD(std::mt19937_64& rngEngine, double* dE)
{
	MCPoly::TrialMoveTAD(rngEngine, dE);
	
	if ( tad->legal )
	{
		MCLiqLattice* liqLat = static_cast<MCLiqLattice*>(lat);

		liqLat->idxSpin1 = tad->en;
		liqLat->idxSpin2 = tad->v2;
		
		*dE += liqLat->GetSpinEnergy();
	}
}

void MCHeteroPoly::AcceptMoveSpinTAD()
{
	MCLiqLattice* liqLat = static_cast<MCLiqLattice*>(lat);
	
	liqLat->AcceptMoveSpin();
	
	MCHeteroPoly::AcceptMoveTAD();
}

void MCHeteroPoly::AcceptMoveTAD()
{
	MCPoly::AcceptMoveTAD();
	
	if ( tadType[tad->n] == 1 )
	{
		tadTable[tad->en] -= 1;
		tadTable[tad->v2] += 1;
	}
}

double MCHeteroPoly::GetSpecificEnergy()
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tad->legal )
	{
		if ( tadType[tad->n] == 1 )
		{
			for ( int i = 0; i < 13; i++ )
			{
				if ( i == 0 )
				{
					E1 -= tadTable[tad->en] - 1.;
					E2 -= tadTable[tad->v2];
				}
				
				else
				{
					int v1 = lat->bitTable[i][tad->en];
					int v2 = lat->bitTable[i][tad->v2];

					E1 -= tadTable[v1];
					E2 -= tadTable[v2];
				}
			}
		}
	}
	
	return Jpp * (E2-E1);
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot])
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int i = 0; i < 12; i++ )
		{
			int v1 = lat->bitTable[i+1][tad->en];
			int v2 = lat->bitTable[i+1][tad->v2];
		
			E1 -= spinTable[v1];
			E2 -= spinTable[v2];
		}
	}
	
	return Jlp * (E2-E1);
}
