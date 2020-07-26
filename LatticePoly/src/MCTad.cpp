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

void MCTad::RandomMove(const int tadConf[Nchain], const int tadNbId[Nchain])
{
	n  = lat->rngEngine() % Nchain;
	en = tadConf[n];
	
	if ( n == 0 )
	{
		int en2 = tadConf[1];
		int cn2 = lat->opp[tadNbId[0]];
		
		int cm2 = std::max(cn2, tadNbId[1]);
		cn2     = std::min(cn2, tadNbId[1]);
		
		iv = lat->rngEngine() % 11;
		
		if ( iv >= cn2 ) iv += 1;
		if ( iv >= cm2 ) iv += 1;
		
		v2 = (iv == 0) ? en2 : lat->bitTable[iv][en2];
		
		int b = lat->bitTable[0][v2];

		legal = ( (b == 0) || ( (b == 1) && (v2 == en2) ) );
		
		if ( legal )
		{
			cn2 = tadNbId[0];
			
			double E1 = lat->cTheta[cn2][tadNbId[1]];
			double E2 = lat->cTheta[lat->opp[iv]][tadNbId[1]];
			
			dE = E2 - E1;
		}
	}
	
	else if ( n == Nchain-1 )
	{
		int en2 = tadConf[Nchain-2];
		int cn2 = tadNbId[Nchain-2];
		int cm2 = std::max(cn2, lat->opp[tadNbId[Nchain-3]]);
		
		cn2     = std::min(cn2, lat->opp[tadNbId[Nchain-3]]);
		
		iv = lat->rngEngine() % 11;
		
		if ( iv >= cn2 ) iv += 1;
		if ( iv >= cm2 ) iv += 1;
		
		v2 = (iv == 0) ? en2 : lat->bitTable[iv][en2];
		
		int b = lat->bitTable[0][v2];
		
		legal = ( (b == 0) || ( (b == 1) && (v2 == en2) ) );
		
		if ( legal )
		{
			cn2 = tadNbId[Nchain-2];
			
			double E1 = lat->cTheta[tadNbId[Nchain-3]][cn2];
			double E2 = lat->cTheta[tadNbId[Nchain-3]][iv];
			
			dE = E2 - E1;
		}
	}
	
	else
	{
		int cn2 = tadNbId[n];
		int cm2 = tadNbId[n-1];
		int en2 = tadConf[n-1];
				
		if ( lat->nbNN[0][cm2][cn2] > 0 )
		{
			iv = lat->rngEngine() % lat->nbNN[0][cm2][cn2];
			
			if ( lat->nbNN[2*iv+1][cm2][cn2] >= cm2 ) iv += 1;
			
			nv1 = lat->nbNN[2*iv+1][cm2][cn2];
			nv2 = lat->nbNN[2*(iv+1)][cm2][cn2];
			
			v2 = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
			
			int b = lat->bitTable[0][v2];

			legal = ( (b == 0) || ( (b == 1) && ( (v2 == en2) || (v2 == tadConf[n+1]) ) ) );
			
			if ( legal )
			{
				double E1,E2;
				
				if ( n == 1 )
				{
					E1 = lat->cTheta[cm2][cn2] + lat->cTheta[cn2][tadNbId[n+1]];
					E2 = lat->cTheta[nv1][nv2] + lat->cTheta[nv2][tadNbId[n+1]];
				}
				
				else if ( n == Nchain-2 )
				{
					E1 = lat->cTheta[tadNbId[n-2]][cm2] + lat->cTheta[cm2][cn2];
					E2 = lat->cTheta[tadNbId[n-2]][nv1] + lat->cTheta[nv1][nv2];
				}
				
				else
				{
					E1 = lat->cTheta[tadNbId[n-2]][cm2] + lat->cTheta[cm2][cn2] + lat->cTheta[cn2][tadNbId[n+1]];
					E2 = lat->cTheta[tadNbId[n-2]][nv1] + lat->cTheta[nv1][nv2] + lat->cTheta[nv2][tadNbId[n+1]];
				}
				
				dE = E2 - E1;
			}
		}
	}
}
