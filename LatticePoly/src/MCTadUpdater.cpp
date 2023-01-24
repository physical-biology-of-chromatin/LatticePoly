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
}

void MCTadUpdater::TrialMoveLeftEnd(const MCTad* tad, double* dE)
{
	MCTad* tad2 = tad->neighbors[1];
	MCBond* bond2 = tad->bonds[1];
	
	int do2 = lat->opp[bond2->dir];
	dn2 = lat->rngEngine() % 11;//1; //

	//if (dn2 !=1 )
	//	return;

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
	dn1 = lat->rngEngine() % 11; //lat->opp[1];// lat->rngEngine() % 11;

	//if (dn1 != lat->opp[1])
	//	return;

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

void MCTadUpdater::TrialMoveTopo(const MCTad* tadi, double* dT)
{
	legalTopo1 = false;
	legalTopo2 = false;
	
	if (lat->bitTable[0][tadi->pos] ==1)
	{
		if (tadi->isLeftEnd())
		{	
			//std::cout << "*****Left*****" << std::endl;
			MCTad* tadi1 = tadi->neighbors[1];
			MCBond* bondi = tadi->bonds[1];

			int dio1 = lat->opp[bondi->dir]; 
			din2 = lat->rngEngine() % 11;

			int dio2 = std::max(dio1, tadi1->bonds[1]->dir);
			dio1     = std::min(dio1, tadi1->bonds[1]->dir);
						
			if ( din2 >= dio1 ) ++din2;
			if ( din2 >= dio2 ) ++din2;
					
			vin = (din2 == 0) ? tadi1->pos : lat->bitTable[din2][tadi1->pos];
			int b = lat->bitTable[0][vin];

			legalTopo1 = ((b == 1) && (vin != tadi1->pos) && ( vin != tadi->pos) );

			if ( !tadi->isFork() && legalTopo1 )
			{
				dio1 = bondi-> dir;
				double Eo = lat->cTheta[dio2][tadi1->bonds[1]->dir];
				double En = lat->cTheta[lat->opp[din2]][tadi1->bonds[1]->dir];
					
				*dT = En - Eo;
			}
		}
		
		else if (tadi->isRightEnd())
		{
			MCTad* tadi1 = tadi->neighbors[0];
			MCBond* bondi = tadi->bonds[0];

			int dio1 = bondi->dir;
			din1 = lat->rngEngine() % 11;

			int dio2 = std::max(dio1, lat->opp[tadi1->bonds[0]->dir]);
			dio1     = std::min(dio1, lat->opp[tadi1->bonds[0]->dir]);
						
			if ( din1 >= dio1 ) ++din1;
			if ( din1 >= dio2 ) ++din1;
					
			vin = (din1 == 0) ? tadi1->pos :lat->bitTable[din1][tadi1->pos];
			int b = lat->bitTable[0][vin];

			legalTopo1 = ((b == 1) && (vin != tadi->neighbors[0]->pos) && ( vin != tadi->pos) );

			if ( legalTopo1 && !tadi->isFork() )
			{
				dio1 = bondi->dir;

				double Eo = lat->cTheta[tadi1->bonds[0]->dir][dio1];
				double En = lat->cTheta[tadi1->bonds[0]->dir][din1];
				
				*dT = En - Eo;
			}
		}		

		else
		{	
			MCTad* tadi1 = tadi->neighbors[0];
			MCTad* tadi2 = tadi->neighbors[1];

			int dio1 = tadi->bonds[0]->dir;
			int dio2 = tadi->bonds[1]->dir;

			if (lat->nbNN[0][dio1][dio2] > 0)
			{
				int iv = lat->rngEngine() % lat->nbNN[0][dio1][dio2];

				if (lat->nbNN[2 * iv + 1][dio1][dio2] >= dio1) ++iv;

				din1 = lat->nbNN[2 * iv + 1][dio1][dio2];
				din2 = lat->nbNN[2 * (iv + 1)][dio1][dio2];

				vin = (din1 == 0) ? tadi1->pos : lat->bitTable[din1][tadi1->pos];
				int b = lat->bitTable[0][vin];

				legalTopo1 = ( (b == 1) && (vin != tadi->neighbors[0]->pos) && (vin != tadi->neighbors[1]->pos) && (vin != tadi->pos) );

				if ( legalTopo1 )
				{
					double Eo = lat->cTheta[dio1][dio2];
					double En = lat->cTheta[din1][din2];
					
					if ( !tadi1->isLeftEnd() && !tadi1->isFork() )
					{
						Eo += lat->cTheta[tadi1->bonds[0]->dir][dio1];
						En += lat->cTheta[tadi1->bonds[0]->dir][din1];
					}

					if ( !tadi2->isRightEnd() && !tadi2->isFork() )
					{
						Eo += lat->cTheta[dio2][tadi2->bonds[1]->dir];
						En += lat->cTheta[din2][tadi2->bonds[1]->dir];
					}
					
					*dT = En - Eo;
				}		
			}
		}	
	}
}		

void MCTadUpdater::TrialSwapTopo(const MCTad* tadi, const MCTad* tadx, double* dT)
{
	if ( tadx->isLeftEnd() )
	{
		MCTad* tadj2 = tadx->neighbors[1];
		MCBond* bond1 = tadx->bonds[1];
		
		int djo2 = lat->opp[bond1->dir];
		
		for (int jv = 0; jv < 11; ++jv)
		{

			int djo1 = std::max(djo2, tadj2->bonds[1]->dir);
			djo2 = std::min(djo2, tadj2->bonds[1]->dir);

			if (jv >= djo2) ++jv;
			if (jv >= djo1) ++jv;

			vjn = (jv == 0) ? tadj2->pos : lat->bitTable[jv][tadj2->pos];
			legalTopo2 = (vjn == tadi->pos);
			
			if (legalTopo2)
			{
				djn2 = jv;
				
				djo2 = bond1->dir;

				double Eo = lat->cTheta[tadj2->bonds[0]->dir][djo2];
				double En = lat->cTheta[tadj2->bonds[0]->dir][djn2];
				
				*dT += En - Eo;
				break;
			}	
										
		}
							
	}
	
	else if ( tadx->isRightEnd() )
	{
		MCTad* tadj1 = tadx->neighbors[0];
		MCBond* bond1 = tadx->bonds[0];
		
		int djo1 = bond1->dir;
		for (int jv = 0; jv < 11; ++jv)
		{
			int djo2 = std::max(djo1, lat->opp[tadj1->bonds[0]->dir]);
			djo1 = std::min(djo1, lat->opp[tadj1->bonds[0]->dir]);
			if (jv >= djo1) ++jv;
			if (jv >= djo2) ++jv;
									
			vjn = (jv == 0) ? tadj1->pos : lat->bitTable[jv][tadj1->pos];
									
			legalTopo2 = (vjn == tadi->pos);
									
			if (legalTopo2)
			{
				djn1 = jv;
				
				djo1 = bond1->dir;

				double Eo = lat->cTheta[tadj1->bonds[0]->dir][djo1];
				double En = lat->cTheta[tadj1->bonds[0]->dir][djn1];
				
				*dT += En - Eo;				
				break;
			}
		}
	}

	else
	{
		MCTad* tadj1 = tadx->neighbors[0];
		MCTad* tadj2 = tadx->neighbors[1];

		int djo1 = tadx->bonds[0]->dir;
		int djo2 = tadx->bonds[1]->dir;

		if (lat->nbNN[0][djo1][djo2] > 0)
		{
			for (int jv = 0; jv < lat->nbNN[0][djo1][djo2]; ++jv)
			{

				if (lat->nbNN[2 * jv + 1][djo1][djo2] >= djo1) ++jv;

				djn1 = lat->nbNN[2 * jv + 1][djo1][djo2];
				djn2 = lat->nbNN[2 * (jv + 1)][djo1][djo2];

				vjn = (djn1 == 0) ? tadj1->pos : lat->bitTable[djn1][tadj1->pos];

				legalTopo2 = (vjn == tadi->pos);

				if (legalTopo2)
				{
					double Eo = lat->cTheta[djo1][djo2];
					double En = lat->cTheta[djn1][djn2];
				
					if ( !tadj1->isLeftEnd() && !tadj1->isFork() )
					{
						Eo += lat->cTheta[tadj1->bonds[0]->dir][djo1];
						En += lat->cTheta[tadj1->bonds[0]->dir][djn1];
					}

					if ( !tadj2->isRightEnd() && !tadj2->isFork() )
					{
						Eo += lat->cTheta[djo2][tadj2->bonds[1]->dir];
						En += lat->cTheta[djn2][tadj2->bonds[1]->dir];
					}
					
					*dT += En - Eo;		
					break;
				}	
			}
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

void MCTadUpdater::AcceptMoveTopo(MCTad* tadi, MCTad* tadx) const
{
	tadi->pos = vin;
	tadx->pos = vjn;

	if ( tadi->isLeftEnd() )
		tadi->bonds[1]->dir = lat->opp[din2];
	
	else if ( tadi->isRightEnd() )
		tadi->bonds[0]->dir = din1;		
	else
	{
		tadi->bonds[0]->dir = din1;
		tadi->bonds[1]->dir = din2;
	}	
	
	if ( tadx->isLeftEnd() )
		tadx->bonds[1]->dir = lat->opp[djn2];

	else if ( tadx->isRightEnd() )	
		tadx->bonds[0]->dir = djn1;
	else
	{	
		tadx->bonds[0]->dir = djn1;
		tadx->bonds[1]->dir = djn2;
	}

}