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
	reptation_values.clear();

	std::vector<int> first_monomer;
	first_monomer.push_back(tad->pos);
	reptation_values.push_back(first_monomer);
	reptation_step=0;

	if ( tad->isLeftEnd() )
		TrialMoveLeftEnd(tad, dE);
		
	
	else if ( tad->isRightEnd() )
		TrialMoveRightEnd(tad, dE);
	
	else if ( tad->isFork() )
		TrialMoveFork(tad, dE);
	
	else
		TrialMoveLinear(tad, dE);
	
	if(legal)
		SaveSpecialMonomers(tad);
}

void MCTadUpdater::TrialMoveLeftEnd(const MCTad* tad, double* dE)
{

	MCTad* tad2 = tad->neighbors[1];
	MCBond* bond2 = tad->bonds[1];
	int vo =reptation_values[reptation_step][0];
	int dn2=lat->rngEngine() % 13;
	
	int vn = dn2==0? vo :lat->bitTable[dn2][vo];
	
	reptation_values[reptation_step].push_back(vn);
	
	
	
	bool legal1 =false;
	
	if ( vn == tad2->pos )
	{
		
		reptation_values[reptation_step].push_back(0);
		
		legal1 = true;
	}
	else
	{
		for ( int v = 0; (v < 12) ; ++v )
		{
			if ( lat->bitTable[v+1][vn] == tad2->pos )
			{
				
				reptation_values[reptation_step].push_back(v+1);
				legal1 = true;
			}
		}
	}
	if(legal1) //standard move
	{
		
		
		int b = lat->bitTable[0][vn];
		legal = legal = (b == 0) || ( (b == 1) && (vn == tad2->pos) );
		if ( legal )
		{
			int do2 = bond2->dir;
			
			if ( !tad2->isFork() )
			{
				double Eo = lat->cTheta[tad2->bonds[0]->dir][do2];
				double En = lat->cTheta[tad2->bonds[0]->dir][reptation_values[reptation_step][2]];
				
				*dE += En - Eo;
			}
		}else
			reptation_values.pop_back();
	}
	else
	{
		int b = lat->bitTable[0][vn];
		
		if((b == 0) || ( (b == 1) && (vn == tad2->pos)))
		{
			reptation_values[reptation_step].push_back(dn2);
			rept_dir=1;
			TrialReptationMove(tad, rept_dir);
			*dE=0;
		}
		else
			reptation_values.pop_back();
	}
	
}

void MCTadUpdater::TrialMoveRightEnd(const MCTad* tad, double* dE)
{

	MCTad* tad1 = tad->neighbors[0];
	MCBond* bond1 = tad->bonds[0];
	int do1 = bond1->dir;
	int vo =reptation_values[reptation_step][0];
	int dn1=lat->rngEngine() % 13;
	
	int vn = dn1==0? vo :lat->bitTable[dn1][vo];


	
	reptation_values[reptation_step].push_back(vn);



	bool legal1 =false;

	if ( vn == tad1->pos )
	{

		reptation_values[reptation_step].push_back(0);

		legal1 = true;
	}
	else
	{
		for ( int v = 0; (v < 12) ; ++v )
		{
			if ( lat->bitTable[v+1][vn] == tad1->pos )
			{

				reptation_values[reptation_step].push_back(lat->opp[v+1]);
				legal1 = true;
			}
		}
	}
	if(legal1) //standard move
	{


		int b = lat->bitTable[0][vn];
		legal = legal = (b == 0) || ( (b == 1) && (vn == tad1->pos) );
		if ( legal )
		{
			do1 = bond1->dir;
			
			if ( !tad1->isFork() )
			{
				double Eo = lat->cTheta[tad1->bonds[0]->dir][do1];
				double En = lat->cTheta[tad1->bonds[0]->dir][reptation_values[reptation_step][2]];
				
				*dE += En - Eo;
			}
		}else
			reptation_values.pop_back();
	}
	else
	{
		int b = lat->bitTable[0][vn];
		//std::cout << "occupancy =  " <<b<< std::endl;

		if((b == 0) || ( (b == 1) && (vn == tad1->pos)))
		{
			reptation_values[reptation_step].push_back(dn1);
			rept_dir=0;
			TrialReptationMove(tad, rept_dir);
			*dE=0;
		}
		else
			reptation_values.pop_back();
	}
}

void MCTadUpdater::TrialMoveLinear(const MCTad* tad, double* dE)
{
	
	MCTad* tad1 = tad->neighbors[0];
	MCTad* tad2 = tad->neighbors[1];

	int do1 = tad->bonds[0]->dir;
	int do2 = tad->bonds[1]->dir;
	
	int dn1=lat->rngEngine() % 13;
	
	int vo =reptation_values[reptation_step][0];
	int vn = dn1==0? vo :lat->bitTable[dn1][vo];

	reptation_values[reptation_step].push_back(vn);


		
	bool legal1 =false;
	bool legal2 =false;
	
	if ( vn == tad1->pos )
	{
		reptation_values[reptation_step].push_back(0);
		
		legal1 = true;
	}
	else
	{
		for ( int v = 0; (v < 12) ; ++v )
		{
			if ( lat->bitTable[v+1][vn] == tad1->pos )
			{
				
				reptation_values[reptation_step].push_back(lat->opp[v+1]);
				legal1 = true;
			}
		}
	}
	if ( vn == tad2->pos )
	{
		if(!legal1)
			reptation_values[reptation_step].push_back(dn1);
		reptation_values[reptation_step].push_back(0);
		legal2 = true;
	}
	else
	{
		for ( int v = 0; (v < 12) ; ++v )
		{
			if ( lat->bitTable[v+1][vn] == tad2->pos )
			{
				if(!legal1)
					reptation_values[reptation_step].push_back(dn1);
				reptation_values[reptation_step].push_back(v+1);
				legal2 = true;
			}
		}
	}
	//in presence of reptation (b=1) moving collectively is alway forbidden but standard move is allowed
	if(legal1!=legal2 and (tad->pos==tad1->pos or tad->pos==tad2->pos))
	{
		legal1=false;
		legal2=false;
		reptation_values.pop_back();
		
	}
		
	if(legal1 and legal2) //standard move
	{
		
		
		int b = lat->bitTable[0][vn];
		legal  = (b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) );
		if ( legal )
		{
			double Eo = lat->cTheta[do1][do2];
			double En = lat->cTheta[reptation_values[reptation_step][2]][reptation_values[reptation_step][3]];
			
			if ( !tad1->isLeftEnd() && !tad1->isFork() )
			{
				Eo += lat->cTheta[tad1->bonds[0]->dir][do1];
				En += lat->cTheta[tad1->bonds[0]->dir][reptation_values[reptation_step][2]];
			}
			
			if ( !tad2->isRightEnd() && !tad2->isFork() )
			{
				Eo += lat->cTheta[do2][tad2->bonds[1]->dir];
				En += lat->cTheta[reptation_values[reptation_step][3]][tad2->bonds[1]->dir];
			}
			
			*dE = En - Eo;
		}else
			reptation_values.pop_back();
		
	
	}
	else if(legal1 and !legal2)
	{
		int b = lat->bitTable[0][vn];
		//std::cout << "occupancy =  " <<b<< std::endl;
		
		if((b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) ))
		{
			reptation_values[reptation_step].push_back(dn1);
			rept_dir=1;

			TrialReptationMove(tad, rept_dir);

			*dE=0;
		}
		else
			reptation_values.pop_back();
	}
	else if(legal2 and !legal1)
	{
		int b = lat->bitTable[0][vn];
		if((b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) ))
		{
			//reptation_values[reptation_step].push_back(dn1); NO has already loaded dn1 as first bond
			rept_dir=0;

			TrialReptationMove(tad, rept_dir);

			*dE=0;
		}
		else
			reptation_values.pop_back();
	}
}

void MCTadUpdater::TrialMoveFork(const MCTad* tad, double* dE)
{/*
	MCTad* tad1 = tad->neighbors[0];
	MCTad* tad2 = tad->neighbors[1];
	MCTad* tad3 = tad->neighbors[2];
	
	int do1 = tad->bonds[0]->dir;
	int do2 = tad->bonds[1]->dir;
	int do3 = tad->bonds[2]->dir;
		
	if ( lat->nbNN[0][do1][do2] > 0 )
	{
		if(legal==false)
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
	}*/
}


void MCTadUpdater::TrialReptationMove(const MCTad* tad, int dir) 
{
	

	if(tad->neighbors[dir]->isLeftEnd() or tad->neighbors[dir]->isRightEnd())
	{
		legal=true;
		int vo= tad->neighbors[dir]->pos;
		int vn = reptation_values[reptation_step][0];
		
		int dn1 = dir==0 ? reptation_values[reptation_step][2] : reptation_values[reptation_step][3] ;
		auto end_monomer={vo, vn ,dn1};
		reptation_values.push_back(end_monomer);
		

	}
	else
	{
		while(legal==false)
		{

			reptation_step+=1;
			MCTad* reptad = tad->neighbors[dir];
			std::vector<int> next_monomer;
			int vo = reptad->pos;
			
			next_monomer.push_back(vo);
			
			reptation_values.push_back(next_monomer);
			
			MCTad* tad1 = reptad->neighbors[dir];
			//MCBond* bond1 = reptad->bonds[dir];
			//int do1 = bond1->dir;

			//int dn1= tad1->bonds[dir==0?1:0]->dir;


			int vn = tad->pos;
			int dn1;
			int dn2;


			reptation_values[reptation_step].push_back(vn);
			

			if ( vn == tad1->pos )
			{
				reptation_values[reptation_step].push_back(0);
				dn1 = dir ==0 ? 0 : tad->bonds[dir]->dir;
				dn2 = dir ==0 ? tad->bonds[dir]->dir : 0;
				legal = true;
			}
			else
			{
				for ( int v = 0; (v < 12) ; ++v )
				{
					if ( lat->bitTable[v+1][vn] == tad1->pos )
					{
						dn1 = dir ==0 ? lat->opp[v+1] : tad->bonds[dir]->dir;
						dn2 = dir ==0 ? tad->bonds[dir]->dir : lat->opp[v+1];

						legal = true;
					}
				}
			}
			if(legal and (!tad1->isRightEnd() or !tad1->isLeftEnd()))
			{
				reptation_values[reptation_step].push_back(dn1);
				reptation_values[reptation_step].push_back(dn2);
				
				vo= tad1->pos;
				vn = tad1->pos;
				
				dn1 = dir==0 ? tad1->bonds[1]->dir : reptation_values[reptation_step][3];
				dn2 = dir==0 ? reptation_values[reptation_step][2]: tad1->bonds[0]->dir;
				auto last_monomer={vo, vn ,dn1,dn2};
				reptation_values.push_back(last_monomer);
			}
			else
			{
				dn1 = dir ==0 ? reptad->bonds[dir]->dir : tad->bonds[dir]->dir;
				dn2 = dir ==0 ? tad->bonds[dir]->dir : reptad->bonds[dir]->dir;

				reptation_values[reptation_step].push_back(dn1);
				reptation_values[reptation_step].push_back(dn2);
			}
			
			
			tad=tad->neighbors[dir];
			
			
			
			if(!legal)
			{
				if(tad->neighbors[dir]->isLeftEnd() or tad->neighbors[dir]->isRightEnd())
				{

					legal=true;

					vo= tad->neighbors[dir]->pos;
					vn = reptation_values[reptation_step][0];

					dn1 = dir==0 ? reptation_values[reptation_step][2] : reptation_values[reptation_step][3] ;
					auto end_monomer={vo, vn ,dn1};
					reptation_values.push_back(end_monomer);
					
					
				}
			}
		}
	}

}
void MCTadUpdater::SaveSpecialMonomers(const MCTad* tad)
{
	std::cout << "SaveSpecialMonomers 1" << std::endl;

	for ( int i = 0; i < (int) reptation_values.size(); ++i )
	{
		auto shifting_tad = tad;
		if(i!=0)
			shifting_tad=tad->neighbors[rept_dir];
		if(shifting_tad->isChoesin)
		{
			std::vector<int> choesin_vector;
			//vo of choesin
			choesin_vector.push_back(reptation_values[i][0]);
			//vn of choesin
			choesin_vector.push_back(reptation_values[i][1]);
			//search for sister choesin among the reptating tads
			for ( int j = 0; j < (int) reptation_values.size(); ++j )
			{
				auto shifting_tad2 = tad;
				if(j!=0)
					shifting_tad2=tad->neighbors[rept_dir];
				if(shifting_tad2->pos==shifting_tad->choesin_binding_site->pos)
				{
					//vo of choesin bindsite
					choesin_vector.push_back(reptation_values[j][0]);
					//vn of choesin bindingsite
					choesin_vector.push_back(reptation_values[j][1]);
				}
				
			}
			if(choesin_vector.size()==2) //I do not find binding site the new and old position of binding sites are the same
			{
				//vo of choesin bindsite
				choesin_vector.push_back(shifting_tad->choesin_binding_site->pos);
				//vn of choesin bindingsite
				choesin_vector.push_back(shifting_tad->choesin_binding_site->pos);
			}
			
		}
	}
	std::cout << "SaveSpecialMonomers 2" << std::endl;

}
void MCTadUpdater::AcceptMove(MCTad* tad) const
{

	for ( int i = 0; i < (int) reptation_values.size(); ++i )
	{
		if(i!=0)
			tad=tad->neighbors[rept_dir];

		tad->pos = reptation_values[i][1];

		if ( tad->isLeftEnd())
			tad->bonds[1]->dir = reptation_values[i][2];
		
		else if ( tad->isRightEnd())
			tad->bonds[0]->dir = reptation_values[i][2];
		else
		{

			tad->bonds[0]->dir = reptation_values[i][2];
			tad->bonds[1]->dir = reptation_values[i][3];

			
			if ( tad->isFork() )
				tad->bonds[2]->dir = reptation_values[i][4];
		}
	}
}
