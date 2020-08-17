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
	dE  = 0.;

	n   = -1;
	nv1 = -1;
	nv2 = -1;
	
	vo  = -1;
	vn  = -1;
		
	legal = false;
}

void MCTad::RandomMove(const int tadConf[Nchain], const int tadBond[Nchain])
{
	n  = lat->rngEngine() % Nchain;
	vo = tadConf[n];
	
	if ( n == 0 )
	{
		int en2 = tadConf[1];
		int cn2 = lat->opp[tadBond[0]];
		
		int cm2 = std::max(cn2, tadBond[1]);
		cn2     = std::min(cn2, tadBond[1]);
		
		nv2 = lat->rngEngine() % 11;
		
		if ( nv2 >= cn2 ) ++nv2;
		if ( nv2 >= cm2 ) ++nv2;
		
		vn = (nv2 == 0) ? en2 : lat->bitTable[nv2][en2];
		
		int b = lat->bitTable[0][vn];

		legal = ( (b == 0) || ( (b == 1) && (vn == en2) ) );
		
		if ( legal )
		{
			cn2 = tadBond[0];
			
			double E1 = lat->cTheta[cn2][tadBond[1]];
			double E2 = lat->cTheta[lat->opp[nv2]][tadBond[1]];
			
			dE = E2 - E1;
		}
	}
	
	else if ( n == Nchain-1 )
	{
		int en2 = tadConf[Nchain-2];
		int cn2 = tadBond[Nchain-2];
		
		int cm2 = std::max(cn2, lat->opp[tadBond[Nchain-3]]);
		cn2     = std::min(cn2, lat->opp[tadBond[Nchain-3]]);
		
		nv1 = lat->rngEngine() % 11;
		
		if ( nv1 >= cn2 ) ++nv1;
		if ( nv1 >= cm2 ) ++nv1;
		
		vn = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
		
		int b = lat->bitTable[0][vn];
		
		legal = ( (b == 0) || ( (b == 1) && (vn == en2) ) );
		
		if ( legal )
		{
			cn2 = tadBond[Nchain-2];
			
			double E1 = lat->cTheta[tadBond[Nchain-3]][cn2];
			double E2 = lat->cTheta[tadBond[Nchain-3]][nv1];
			
			dE = E2 - E1;
		}
	}
	
	else
	{
		int cn2 = tadBond[n];
		int cm2 = tadBond[n-1];
		int en2 = tadConf[n-1];
				
		if ( lat->nbNN[0][cm2][cn2] > 0 )
		{
			int iv = lat->rngEngine() % lat->nbNN[0][cm2][cn2];
			
			if ( lat->nbNN[2*iv+1][cm2][cn2] >= cm2 ) ++iv;
			
			nv1 = lat->nbNN[2*iv+1][cm2][cn2];
			nv2 = lat->nbNN[2*(iv+1)][cm2][cn2];
			
			vn = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
			
			int b = lat->bitTable[0][vn];

			legal = ( (b == 0) || ( (b == 1) && ( (vn == en2) || (vn == tadConf[n+1]) ) ) );
			
			if ( legal )
			{
				double E1,E2;
				
				if ( n == 1 )
				{
					E1 = lat->cTheta[cm2][cn2] + lat->cTheta[cn2][tadBond[n+1]];
					E2 = lat->cTheta[nv1][nv2] + lat->cTheta[nv2][tadBond[n+1]];
				}
				
				else if ( n == Nchain-2 )
				{
					E1 = lat->cTheta[tadBond[n-2]][cm2] + lat->cTheta[cm2][cn2];
					E2 = lat->cTheta[tadBond[n-2]][nv1] + lat->cTheta[nv1][nv2];
				}
				
				else
				{
					E1 = lat->cTheta[tadBond[n-2]][cm2] + lat->cTheta[cm2][cn2] + lat->cTheta[cn2][tadBond[n+1]];
					E2 = lat->cTheta[tadBond[n-2]][nv1] + lat->cTheta[nv1][nv2] + lat->cTheta[nv2][tadBond[n+1]];
				}
				
				dE = E2 - E1;
			}
		}
	}
}
