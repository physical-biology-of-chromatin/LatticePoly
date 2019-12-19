//
//  MCTad.cpp
//  LatticePoly
//
//  Created by mtortora on 02/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCTad.hpp"


MCTad::MCTad(MCLattice* _lat): lat(_lat) {}

void MCTad::Init()
{
	n   = -1;
	en  = -1;
	v2  = -1;
	iv  = -1;
	nv1 = -1;
	nv2 = -1;
	
	dE = 0.;
	legal = false;
}

void MCTad::RandomMove(std::mt19937_64& rngEngine, const int config[2][Nchain])
{
	n  = rngEngine() % Nchain;
	en = config[0][n];
	
	if ( n == 0 )
	{
		int en2 = config[0][1];
		int cn2 = lat->opp[config[1][0]];
		
		int cm2 = std::max(cn2, config[1][1]);
		cn2     = std::min(cn2, config[1][1]);
		
		iv = rngEngine() % 11;
		
		if ( iv >= cn2 ) iv += 1;
		if ( iv >= cm2 ) iv += 1;
		
		v2 = (iv == 0) ? en2 : lat->bitTable[iv][en2];
		
		int b = lat->bitTable[0][v2];

		legal = ( (b == 0) || ( (b == 1) && (v2 == en2) ) );
		
		if ( legal )
		{
			cn2 = config[1][0];
			
			double E1 = lat->cTheta[cn2][config[1][1]];
			double E2 = lat->cTheta[lat->opp[iv]][config[1][1]];
			
			dE = E2 - E1;
		}
	}
	
	else if ( n == Nchain-1 )
	{
		int en2 = config[0][Nchain-2];
		int cn2 = config[1][Nchain-2];
		
		int cm2 = std::max(cn2, lat->opp[config[1][Nchain-3]]);
		cn2     = std::min(cn2, lat->opp[config[1][Nchain-3]]);
		
		iv = rngEngine() % 11;
		
		if ( iv >= cn2 ) iv += 1;
		if ( iv >= cm2 ) iv += 1;
		
		v2 = (iv == 0) ? en2 : lat->bitTable[iv][en2];
		
		int b = lat->bitTable[0][v2];
		
		legal = ( (b == 0) || ( (b == 1) && (v2 == en2) ) );
		
		if ( legal )
		{
			cn2 = config[1][Nchain-2];
			
			double E1 = lat->cTheta[config[1][Nchain-3]][cn2];
			double E2 = lat->cTheta[config[1][Nchain-3]][iv];
			
			dE = E2 - E1;
		}
	}
	
	else
	{
		int cn2 = config[1][n];
		int cm2 = config[1][n-1];
		
		int en2 = config[0][n-1];
				
		if ( lat->nbNN[0][cm2][cn2] > 0 )
		{
			iv = rngEngine() % lat->nbNN[0][cm2][cn2];
			
			if ( lat->nbNN[2*iv+1][cm2][cn2] >= cm2 ) iv += 1;
			
			nv1 = lat->nbNN[2*iv+1][cm2][cn2];
			nv2 = lat->nbNN[2*(iv+1)][cm2][cn2];
			
			v2 = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
			
			int b = lat->bitTable[0][v2];

			legal = ( (b == 0) || ( (b == 1) && ( (v2 == en2) || (v2 == config[0][n+1]) ) ) );
			
			if ( legal )
			{
				double E1,E2;
				
				if ( n == 1 )
				{
					E1 = lat->cTheta[cm2][cn2] + lat->cTheta[cn2][config[1][n+1]];
					E2 = lat->cTheta[nv1][nv2] + lat->cTheta[nv2][config[1][n+1]];
				}
				
				else if ( n == Nchain-2 )
				{
					E1 = lat->cTheta[config[1][n-2]][cm2] + lat->cTheta[cm2][cn2];
					E2 = lat->cTheta[config[1][n-2]][nv1] + lat->cTheta[nv1][nv2];
				}
				
				else
				{
					E1 = lat->cTheta[config[1][n-2]][cm2] + lat->cTheta[cm2][cn2] + lat->cTheta[cn2][config[1][n+1]];
					E2 = lat->cTheta[config[1][n-2]][nv1] + lat->cTheta[nv1][nv2] + lat->cTheta[nv2][config[1][n+1]];
				}
				
				dE = E2 - E1;
			}
		}
	}
}
