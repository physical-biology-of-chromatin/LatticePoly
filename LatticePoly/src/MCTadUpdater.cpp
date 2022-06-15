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
	
	else if ( tad->isFork() )
		TrialMoveFork(tad, dE);
	
	else
		TrialMoveLinear(tad, dE);
	
//	if (Rconfinement > 0)
//	{
//		double c = (L-0.5)/2;
//		double d2 = SQR(lat->xyzTable[0][vn]-c)+SQR(lat->xyzTable[1][vn]-c)+SQR(lat->xyzTable[2][vn]-c);
//
//		if ( d2 > SQR(Rconfinement) ) legal = false;
//	}
}

void MCTadUpdater::TrialMoveLeftEnd(const MCTad* tad, double* dE)
{
	MCTad* tad2 = tad->neighbors[1];
	MCBond* bond2 = tad->bonds[1];
	
	int do2 = lat->opp[bond2->dir];
	dn2 = lat->rngEngine() % 11;

	int do1 = std::max(do2, tad2->bonds[1]->dir);
	do2     = std::min(do2, tad2->bonds[1]->dir);
		
	if ( dn2 >= do2 ) ++dn2;
	if ( dn2 >= do1 ) ++dn2;
	
	vn = (dn2 == 0) ? tad2->pos : lat->bitTable[dn2][tad2->pos];
	int b = lat->bitTable[0][vn];

	legal = (b == 0) || ( (b == 1) && (vn == tad2->pos) );
	
	if ( legal )
	{
		do2 = bond2->dir;

		if ( !tad2->isFork() )
		{
			double Eo = lat->cTheta[do2][tad2->bonds[1]->dir];
			double En = lat->cTheta[lat->opp[dn2]][tad2->bonds[1]->dir];
				
			*dE = En - Eo;
		}
	}
}

void MCTadUpdater::TrialMoveRightEnd(const MCTad* tad, double* dE)
{
	MCTad* tad1 = tad->neighbors[0];
	MCBond* bond1 = tad->bonds[0];
	
	int do1 = bond1->dir;
	dn1 = lat->rngEngine() % 11;

	int do2 = std::max(do1, lat->opp[tad1->bonds[0]->dir]);
	do1     = std::min(do1, lat->opp[tad1->bonds[0]->dir]);
		
	if ( dn1 >= do1 ) ++dn1;
	if ( dn1 >= do2 ) ++dn1;
	
	vn = (dn1 == 0) ? tad1->pos : lat->bitTable[dn1][tad1->pos];
	int b = lat->bitTable[0][vn];
	
	legal = (b == 0) || ( (b == 1) && (vn == tad1->pos) );
	
	if ( legal )
	{
		do1 = bond1->dir;

		if ( !tad1->isFork() )
		{
			double Eo = lat->cTheta[tad1->bonds[0]->dir][do1];
			double En = lat->cTheta[tad1->bonds[0]->dir][dn1];
				
			*dE = En - Eo;
		}
	}
}

void MCTadUpdater::TrialMoveLinear(const MCTad* tad, double* dE)
{
	MCTad* tad1 = tad->neighbors[0];
	MCTad* tad2 = tad->neighbors[1];

	int do1 = tad->bonds[0]->dir;
	int do2 = tad->bonds[1]->dir;
			
	if ( lat->nbNN[0][do1][do2] > 0 )
	{
		int iv = lat->rngEngine() % lat->nbNN[0][do1][do2];
		
		if ( lat->nbNN[2*iv+1][do1][do2] >= do1 ) ++iv;
		
		dn1 = lat->nbNN[2*iv+1][do1][do2];
		dn2 = lat->nbNN[2*(iv+1)][do1][do2];
		
		vn = (dn1 == 0) ? tad1->pos : lat->bitTable[dn1][tad1->pos];
		int b = lat->bitTable[0][vn];

		legal = (b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) );
		
		if ( legal )
		{
			double Eo = lat->cTheta[do1][do2];
			double En = lat->cTheta[dn1][dn2];
			
			if ( !tad1->isLeftEnd() && !tad1->isFork() )
			{
				Eo += lat->cTheta[tad1->bonds[0]->dir][do1];
				En += lat->cTheta[tad1->bonds[0]->dir][dn1];
			}

			if ( !tad2->isRightEnd() && !tad2->isFork() )
			{
				Eo += lat->cTheta[do2][tad2->bonds[1]->dir];
				En += lat->cTheta[dn2][tad2->bonds[1]->dir];
			}
			
			*dE = En - Eo;
		}
	}
}

void MCTadUpdater::TrialMoveFork(const MCTad* tad, double* dE)
{
	MCTad* tad1 = tad->neighbors[0];
	MCTad* tad2 = tad->neighbors[1];
	MCTad* tad3 = tad->neighbors[2];
	
	int do1 = tad->bonds[0]->dir;
	int do2 = tad->bonds[1]->dir;
	int do3 = tad->bonds[2]->dir;
		
	if ( lat->nbNN[0][do1][do2] > 0 )
	{
		// Pick new position compatible with bonds 1 & 2
		int iv = lat->rngEngine() % lat->nbNN[0][do1][do2];
		
		if ( lat->nbNN[2*iv+1][do1][do2] >= do1 ) ++iv;
		
		dn1 = lat->nbNN[2*iv+1][do1][do2];
		dn2 = lat->nbNN[2*(iv+1)][do1][do2];
		
		vn = (dn1 == 0) ? tad1->pos : lat->bitTable[dn1][tad1->pos];
		int b = lat->bitTable[0][vn];

		// Check if new position vn complies with occupancy criteria
		bool legal1 = (b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) || (vn == tad3->pos) ) );
		
		// Check if new position is compatible with bond 3 (i.e., vn should be a nearest neighbor of tad 3)
		bool legal2 = false;
		
		if ( legal1 )
		{
			if ( vn == tad3->pos )
			{
				dn3 = 0;
				legal2 = true;
			}
			
			else
			{
				for ( int v = 0; (v < 12) && (!legal2); ++v )
				{
					if ( lat->bitTable[v+1][vn] == tad3->pos )
					{
						// Reverse bond orientation between right fork and rightmost replicated tad for consistency
						dn3 = tad->isRightFork() ? lat->opp[v+1] : v+1;
						legal2 = true;
					}
				}
			}
		}
			
		legal = legal1 && legal2;

		// Compute bending energies, assuming a 0 bending modulus for the forks
		if ( legal )
		{
			double Eo = 0.;
			double En = 0.;

			if ( !tad1->isLeftEnd() && !tad1->isFork() )
			{
				Eo += lat->cTheta[tad1->bonds[0]->dir][do1];
				En += lat->cTheta[tad1->bonds[0]->dir][dn1];
			}

			if ( !tad2->isRightEnd() && !tad2->isFork() )
			{
				Eo += lat->cTheta[do2][tad2->bonds[1]->dir];
				En += lat->cTheta[dn2][tad2->bonds[1]->dir];
			}
				
			if ( tad == tad3->neighbors[1] )
			{
				Eo += lat->cTheta[tad3->bonds[0]->dir][do3];
				En += lat->cTheta[tad3->bonds[0]->dir][dn3];
			}
			
			else
			{
				Eo += lat->cTheta[do3][tad3->bonds[1]->dir];
				En += lat->cTheta[dn3][tad3->bonds[1]->dir];
			}
			
			*dE = En - Eo;
		}
	}
}

void MCTadUpdater::AcceptMove(MCTad* tad) const
{
	tad->pos = vn;

	if ( tad->isLeftEnd() )
		tad->bonds[1]->dir = lat->opp[dn2];
	
	else if ( tad->isRightEnd() )
		tad->bonds[0]->dir = dn1;
	
	else
	{
		tad->bonds[0]->dir = dn1;
		tad->bonds[1]->dir = dn2;
		
		if ( tad->isFork() )
			tad->bonds[2]->dir = dn3;
	}
}
