//
//  MCHeteroPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCHeteroPoly.hpp"


MCHeteroPoly::MCHeteroPoly(MCLattice* _lat): MCPoly(_lat) {};

void MCHeteroPoly::Init()
{
	MCPoly::Init();
		
	for ( int i = 0; i < Ntot; i++ )
		tadHetTable[i] = 0;
	
	for ( int i = 0; i < Ndom; i++ )
	{
		int idx = 0;

		for ( int j = 0; j < Nloc; j++ )
		{
			tadType[idx+j] = 1;
			tadHetTable[tadConf[idx+j]] = 1;
		}
	}
}

void MCHeteroPoly::AcceptMove()
{
	MCPoly::AcceptMove();
	
	if ( tadType[tad->n] == 1 )
	{
		tadHetTable[tad->en] -= 1;
		tadHetTable[tad->v2] += 1;
	}
}

double MCHeteroPoly::GetSpecificEnergy() const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int i = 0; i < 13; i++ )
		{
			if ( i == 0 )
			{
				E1 -= tadHetTable[tad->en] - 1.;
				E2 -= tadHetTable[tad->v2];
			}
			
			else
			{
				E1 -= tadHetTable[lat->bitTable[i][tad->en]];
				E2 -= tadHetTable[lat->bitTable[i][tad->v2]];
			}
		}
	}
	
	return Jpp * (E2-E1);
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int i = 0; i < 13; i++ )
		{
			int v1 = (i == 0) ? tad->en : lat->bitTable[i][tad->en];
			int v2 = (i == 0) ? tad->v2 : lat->bitTable[i][tad->v2];
			
			E1 -= spinTable[v1];
			E2 -= spinTable[v2];
		}
	}
	
	return Jlp * (E2-E1);
}
