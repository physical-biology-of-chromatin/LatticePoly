//
//  MCTadUpdater.cpp
//  LatticePoly
//
//  Created by mtortora on 01/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCTadUpdater.hpp"


MCTadUpdater::MCTadUpdater(MCLattice* _lat): lat(_lat) {}

void MCTadUpdater::TrialMove(const MCTad* tad, double* dE)
{
	*dE = 0;
	legal = false;
	
	vo = tad->pos;

	if ( tad->isLeftEnd() )
		TrialMoveLeftEnd(tad, dE);
	
	else if ( tad->isRightEnd() )
		TrialMoveRightEnd(tad, dE);
	
	else
	{
		if ( !tad->isFork() )
			TrialMoveLinear(tad, dE);
	}
}

void MCTadUpdater::TrialMoveLeftEnd(const MCTad* tad, double* dE)
{
	MCTad* nb = tad->neighbors[0];
	MCBond* bond = tad->bonds[0];
	
	int en2 = nb->pos;
	int cn2 = lat->opp[bond->dir];
	
	int cm2 = std::max(cn2, nb->bonds[1]->dir);
	cn2     = std::min(cn2, nb->bonds[1]->dir);
	
	nv2 = lat->rngEngine() % 11;
	
	if ( nv2 >= cn2 ) ++nv2;
	if ( nv2 >= cm2 ) ++nv2;
	
	vn = (nv2 == 0) ? en2 : lat->bitTable[nv2][en2];
	
	int b = lat->bitTable[0][vn];

	legal = ( (b == 0) || ( (b == 1) && (vn == en2) ) );
	
	if ( legal )
	{
		cn2 = bond->dir;

		for ( int b1 = 0; b1 < nb->links; ++b1 )
		{
			if ( nb->bonds[b1] != bond )
			{
				double E1 = lat->cTheta[cn2][nb->bonds[b1]->dir];
				double E2 = lat->cTheta[lat->opp[nv2]][nb->bonds[b1]->dir];
				
				*dE += E2 - E1;
			}
		}
	}
}

void MCTadUpdater::TrialMoveRightEnd(const MCTad* tad, double* dE)
{
	MCTad* nb = tad->neighbors[0];
	MCBond* bond = tad->bonds[0];
	
	int en2 = nb->pos;
	int cn2 = bond->dir;
	
	int cm2 = std::max(cn2, lat->opp[nb->bonds[0]->dir]);
	cn2     = std::min(cn2, lat->opp[nb->bonds[0]->dir]);
	
	nv1 = lat->rngEngine() % 11;
	
	if ( nv1 >= cn2 ) ++nv1;
	if ( nv1 >= cm2 ) ++nv1;
	
	vn = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
	
	int b = lat->bitTable[0][vn];
	
	legal = ( (b == 0) || ( (b == 1) && (vn == en2) ) );
	
	if ( legal )
	{
		cn2 = bond->dir;

		for ( int b1 = 0; b1 < nb->links; ++b1 )
		{
			if ( nb->bonds[b1] != bond )
			{
				double E1 = lat->cTheta[nb->bonds[b1]->dir][cn2];
				double E2 = lat->cTheta[nb->bonds[b1]->dir][nv1];
				
				*dE += E2 - E1;
			}
		}
	}
}

void MCTadUpdater::TrialMoveLinear(const MCTad* tad, double* dE)
{
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];

	MCBond* bond1 = tad->bonds[0];
	MCBond* bond2 = tad->bonds[1];

	int cm2 = bond1->dir;
	int cn2 = bond2->dir;

	int en2 = nb1->pos;
			
	if ( lat->nbNN[0][cm2][cn2] > 0 )
	{
		int iv = lat->rngEngine() % lat->nbNN[0][cm2][cn2];
		
		if ( lat->nbNN[2*iv+1][cm2][cn2] >= cm2 ) ++iv;
		
		nv1 = lat->nbNN[2*iv+1][cm2][cn2];
		nv2 = lat->nbNN[2*(iv+1)][cm2][cn2];
		
		vn = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
		
		int b = lat->bitTable[0][vn];

		legal = ( (b == 0) || ( (b == 1) && ( (vn == en2) || (vn == nb2->pos) ) ) );
		
		if ( legal )
		{
			double E1 = 0.;
			double E2 = 0.;
			
			if ( nb1->isLeftEnd() )
			{
				for ( int b2 = 0; b2 < nb2->links; ++b2 )
				{
					if ( nb2->bonds[b2] != bond2 )
					{
						E1 += lat->cTheta[cm2][cn2] + lat->cTheta[cn2][nb2->bonds[b2]->dir];
						E2 += lat->cTheta[nv1][nv2] + lat->cTheta[nv2][nb2->bonds[b2]->dir];
					}
				}
			}
			
			else if ( nb2->isRightEnd() )
			{
				for ( int b1 = 0; b1 < nb1->links; ++b1 )
				{
					if ( nb1->bonds[b1] != bond1 )
					{
						E1 += lat->cTheta[nb1->bonds[b1]->dir][cm2] + lat->cTheta[cm2][cn2];
						E2 += lat->cTheta[nb1->bonds[b1]->dir][nv1] + lat->cTheta[nv1][nv2];
					}
				}
			}
			
			else
			{
				for ( int b1 = 0; b1 < nb1->links; ++b1 )
				{
					if ( nb1->bonds[b1] != bond1 )
					{
						for ( int b2 = 0; b2 < nb2->links; ++b2 )
						{
							if ( nb2->bonds[b2] != bond2 )
							{
								E1 += lat->cTheta[nb1->bonds[b1]->dir][cm2] + lat->cTheta[cm2][cn2] + lat->cTheta[cn2][nb2->bonds[b2]->dir];
								E2 += lat->cTheta[nb1->bonds[b1]->dir][nv1] + lat->cTheta[nv1][nv2] + lat->cTheta[nv2][nb2->bonds[b2]->dir];
							}
						}
					}
				}
			}
			
			*dE = E2 - E1;
		}
	}
}

void MCTadUpdater::AcceptMove(MCTad* tad) const
{
	tad->pos = vn;

	if ( tad->isLeftEnd() )
		tad->bonds[0]->dir = lat->opp[nv2];
	
	else if ( tad->isRightEnd() )
		tad->bonds[0]->dir = nv1;
	
	else
	{
		if ( !tad->isFork() )
		{
			tad->bonds[0]->dir = nv1;
			tad->bonds[1]->dir = nv2;
		}
	}
}
