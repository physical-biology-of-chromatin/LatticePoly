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
		tadHetTable[i] = 0;
	
	for ( int i = 0; i < Ndom; i++ )
	{
		int idx = rngEngine() % (Nchain-Nloc+1);

		for ( int j = 0; j < Nloc; j++ )
		{
			tadType[idx+j] = 1;
			tadHetTable[tadConf[idx+j]] = 1;
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
	MCHeteroPoly::AcceptMoveTAD();
	
	static_cast<MCLiqLattice*>(lat)->AcceptMoveSpin();
}

void MCHeteroPoly::AcceptMoveTAD()
{
	MCPoly::AcceptMoveTAD();
	
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
	
	if ( tad->legal )
	{
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
					int v1 = lat->bitTable[i][tad->en];
					int v2 = lat->bitTable[i][tad->v2];

					E1 -= tadHetTable[v1];
					E2 -= tadHetTable[v2];
				}
			}
		}
	}
	
	return Jpp * (E2-E1);
}

double MCHeteroPoly::GetBindingEnergy(const int spinTable[Ntot]) const
{
	double E1 = 0.;

	for ( int i = 0; i < 13; i++ )
	{
		if ( i == 0 )
		{
			E1 -= spinTable[tad->v2]*tadHetTable[tad->v2];
			E1 -= spinTable[tad->en]*(tadHetTable[tad->en]-tadType[tad->n]);
		}
		
		else
		{
			int v1 = lat->bitTable[i][tad->en];
			int v2 = lat->bitTable[i][tad->v2];
		
			E1 -= spinTable[tad->en] * tadHetTable[v1];
			E1 -= spinTable[tad->v2] * ( (v2 == tad->en) ? tadHetTable[v2]-tadType[tad->n] : tadHetTable[v2] );

			if ( tadType[tad->n] == 1 )
				E1 -= (v1 == tad->v2) ? 0 : spinTable[v1];
		}
	}
	
	return Jlp * E1;
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	for ( int i = 0; i < 13; i++ )
	{
		if ( i == 0 )
		{
			if ( spinTable[tad->en] != spinTable[tad->v2] )
			{
				E1 -= spinTable[tad->v2]*tadHetTable[tad->v2];
				E2 -= spinTable[tad->en]*tadHetTable[tad->v2];

				E1 -= spinTable[tad->en]*(tadHetTable[tad->en]-tadType[tad->n]);
				E2 -= spinTable[tad->v2]*(tadHetTable[tad->en]-tadType[tad->n]);		
			}	
		}

		else
		{
			int v1 = lat->bitTable[i][tad->en];
			int v2 = lat->bitTable[i][tad->v2];
		
			if ( spinTable[tad->en] != spinTable[tad->v2] )
			{
				E1 -= spinTable[tad->en] * tadHetTable[v1];
				E2 -= spinTable[tad->v2] * tadHetTable[v1];

				E1 -= spinTable[tad->v2] * ( (v2 == tad->en) ? tadHetTable[v2]-tadType[tad->n] : tadHetTable[v2] );
				E2 -= spinTable[tad->en] * ( (v2 == tad->en) ? tadHetTable[v2]-tadType[tad->n] : tadHetTable[v2] );
			}
		
			if ( tadType[tad->n] == 1 )
			{		
				E1 -= (v1 == tad->v2) ? 0 : spinTable[v1];
				E2 -= (v2 == tad->en) ? 0 : spinTable[v2];
			}
		}
	}
	
	return Jlp * (E2-E1);
}
