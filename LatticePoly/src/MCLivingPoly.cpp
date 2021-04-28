//
//  MCLivingPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 27/04/2021.
//  Copyright © 2021 ENS Lyon. All rights reserved.
//

#include "MCLivingPoly.hpp"


MCLivingPoly::MCLivingPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCLivingPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);

	if ( !RestartFromFile )
	{
		for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
		{
			if ( tad->type == 1 )
			{
				double rnd = lat->rngDistrib(lat->rngEngine);
				
				if ( rnd < inactiveRatio )
				{
					tad->type = 2;

					for ( int v = 0; v < 13; ++v )
					{
						int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
						
						--hetTable[vi];
					}
				}
			}
		}
	}
}

void MCLivingPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	
	if ( tadTrial->type == 2 )
	{
		double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < propRate / ((double) Ninter*Nmeas) )
		{
			tadTrial->type = 1;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
				++hetTable[vi];
			}
		}
	}
}
