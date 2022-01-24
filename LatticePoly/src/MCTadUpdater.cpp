//
//  MCTadUpdater.cpp
//  LatticePoly
//
//  Created by mtortora on 01/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCTadUpdater.hpp"

//MAIN INSTRUCTION: every Trial Move (if succesfull) will update a vector "reptation_values".
//reptation_values contains for each monomer displaced a vector containing : < vo, vn, dn1, dn2, dn3, dir > with dir is the direction needed in order to recover the dislaced monomer from the previus one. The vector can contain only 1 (Left/RigtEnd) or 2 (linear monomr) dn.
//reptating_tads saves the tads involved in the collective movement. NB when I start reptation (including the case I try to move a fork during reptation) I must include one tad's neighbour/neighbours

// MAIN RULE: I can reptate in at most ONE direction

MCTadUpdater::MCTadUpdater(MCLattice* _lat): lat(_lat) {}

void MCTadUpdater::TrialMove(const MCTad* tad, double* dE)
{

	*dE = 0;
	legal = false;
	reptation_values.clear();
	reptating_tads.clear();


	
	
	if ( tad->isLeftEnd() )
		TrialMoveLeftEnd(tad);
		
	
	else if ( tad->isRightEnd() )
		TrialMoveRightEnd(tad);
	
	else if ( tad->isFork() )
		TrialMoveFork(tad );
	
	else
	{
		if(tad->isChoesin!=true)//standard linear move
			TrialMoveLinear(tad);
		else // Choesin binded sites acts as replication forks where the third bonds is with the other anchor
			TrialMoveFork(tad);
	}
}

void MCTadUpdater::TrialMoveLeftEnd(const MCTad* tad)
{

	MCTad* tad2 = tad->neighbors[1];
	int vo = tad->pos;
	int rndir=lat->rngEngine() % 13;
	
	int vn = rndir==0? vo :lat->bitTable[rndir][vo];
	int b = lat->bitTable[0][vn];

	int dn1 = -1;
	
	if ( vn == tad2->pos )
		dn1=0;
	for ( int v = 0; (v < 12) ; ++v )
	{
		if ( lat->bitTable[v+1][vn] == tad2->pos )
			dn1=v+1;
	}
	
	if(dn1 != -1) //standard move no need to save the direction since I won't reptate
	{
		legal = (b == 0) || ( (b == 1) && (vn == tad2->pos) );
		if ( legal )
		{
			reptating_tads.push_back(tad);
			reptation_values.push_back({vo,vn,dn1});
			
		}
		return;
	}
	else if ( (dn1 = -1) and (( b == 0) || ( (b == 1) && (vn == tad2->pos))) )
	{
		//if I break connectivity I can reptate in the only available direction
		reptation_values.push_back({vo,vn,lat->opp[rndir]});
		TrialReptationMove(tad, 1);
	}
	else
		reptation_values.clear();

}

void MCTadUpdater::TrialMoveRightEnd(const MCTad* tad)
{

	MCTad* tad1 = tad->neighbors[0];
	int vo = tad->pos;
	int rndir=lat->rngEngine() % 13;

	int vn = rndir==0? vo :lat->bitTable[rndir][vo];
	int b = lat->bitTable[0][vn];
	
	int dn1 = -1;
	
	if ( vn == tad1->pos )
		dn1=0;
	for ( int v = 0; (v < 12) ; ++v )
	{
		if ( lat->bitTable[v+1][vn] == tad1->pos )
			dn1=lat->opp[v+1];
	}
	
	if(dn1 != -1) //standard move
	{
		legal = (b == 0) || ( (b == 1) && (vn == tad1->pos) );
		if ( legal )
		{
			reptating_tads.push_back(tad);
			reptation_values.push_back({vo,vn,dn1});
		}
		return;
	}
	else if ( (dn1 = -1) and (( b == 0) || ( (b == 1) && (vn == tad1->pos))) )
	{
		reptating_tads.push_back(tad);
		reptation_values.push_back({vo,vn,rndir});
		TrialReptationMove(tad, 0);
	}
	else
		reptation_values.clear();

}

void MCTadUpdater::TrialMoveLinear(const MCTad* tad)
{

	MCTad* tad1 = tad->neighbors[0];
	MCTad* tad2 = tad->neighbors[1];

	
	int vo = tad->pos;
	int rndir = lat->rngEngine() % 13;
	
	int vn = rndir==0? vo :lat->bitTable[rndir][vo];
	int b = lat->bitTable[0][vn];
	
	int dn1 = -1;
	int dn2 = -1;


	if ( vn == tad1->pos )
		dn1 = 0;
	if ( vn == tad2->pos)
		dn2 = 0;
	
	for ( int v = 0; (v < 12) ; ++v )
	{
		if ( lat->bitTable[v+1][vn] == tad1->pos )
			dn1 = lat->opp[v+1];
		if ( lat->bitTable[v+1][vn] == tad2->pos )
			dn2 = v+1;
	}
	
	
	if(dn1 != -1 and dn2 != -1) //standard move
	{
		legal  = (b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) );
		if ( legal )
		{
			reptating_tads.push_back(tad);
			reptation_values.push_back({vo,vn,dn1,dn2});
		}
		return;
	
	}
	
	//in presence of reptation (b=1) moving collectively is alway forbidden: it would lead to double occupancy of non-consecutive monomers
	else if (tad->pos==tad1->pos or tad->pos==tad2->pos)
		return;


	else
	{
		//if only one of the two bons mantains connectivity, reptate to opposite direction
		if((b == 0) || ( (b == 1) && ( (vn == tad1->pos) || (vn == tad2->pos) ) ))
		{
			if(dn1!=-1) //reptating right direction (+1)
			{
				reptating_tads.push_back(tad);
				reptating_tads.push_back(tad->neighbors[0]);
				reptation_values.push_back({vo, vn, dn1, rndir});
				TrialReptationMove(tad, 1);

			}
			else if(dn2!=-1)//reptating right direction (+0)
			{
				reptating_tads.push_back(tad);
				reptating_tads.push_back(tad->neighbors[1]);
				reptation_values.push_back({vo, vn, lat->opp[rndir], dn2});
				TrialReptationMove(tad, 0);
			}
		}
	}
}

void MCTadUpdater::TrialMoveFork(const MCTad* tad)
{
	
	int tad1_pos = tad->neighbors[0]->pos;
	int tad2_pos = tad->neighbors[1]->pos;
	int tad3_pos = tad->neighbors[2]->pos;
	
	int do1 = tad->bonds[0]->dir;
	int do2 = tad->bonds[1]->dir;
	int do3 = tad->bonds[2]->dir;
	std::vector<int> compatible_moves;
	
	if ( lat->nbNN[0][do1][do2] > 0 )
	{
		
		// Pick new position compatible with bonds 1 & 2
		for ( int iv = 0; (iv < lat->nbNN[0][do1][do2]) ; ++iv )
		{
			if ( lat->nbNN[2*iv+1][do1][do2] >= do1 ) ++iv;
		
			int dn1 = lat->nbNN[2*iv+1][do1][do2];
			compatible_moves.push_back(dn1);
		}
	}
	if ( lat->nbNN[0][do2][do3] > 0 )
	{
		// Pick new position compatible with bonds 2 & 3
		for ( int iv = 0; (iv < lat->nbNN[0][do2][do3]) ; ++iv )
		{
			if ( lat->nbNN[2*iv+1][do2][do3] >= do2 ) ++iv;
			
			int dn1 = lat->nbNN[2*iv+1][do2][do3];
			compatible_moves.push_back(dn1);
		}
		
	}
	if ( lat->nbNN[0][do1][do3] > 0 )
	{
		// Pick new position compatible with bonds 1 & 3
		for ( int iv = 0; (iv < lat->nbNN[0][do1][do3]) ; ++iv )
		{
			if ( lat->nbNN[2*iv+1][do1][do3] >= do1 ) ++iv;
			
			int dn1 = lat->nbNN[2*iv+1][do1][do3];
			compatible_moves.push_back(dn1);
		}
		
	}
	int rndir;
	if(compatible_moves.size()>0)
		rndir=compatible_moves[lat->rngEngine() % (int)compatible_moves.size()];
	else
		return;
	
	
	//first move rept_dir is not used. Thus is set to -1
	CheckForkLegal(tad,tad1_pos,tad2_pos,tad3_pos,rndir, -1);
	
	
}

void MCTadUpdater::CheckForkLegal( const MCTad* tad , int tad1_pos, int tad2_pos, int tad3_pos, int dir, int rept_dir)
{
	//function to check the connectivity between 3 point: dir will give me the direction towards MCTad* tad is trying to move to
	
	int dn1=-1;
	int dn2=-1;
	int dn3= -1;

	int vn = dir==0? tad->pos : lat->bitTable[dir][tad->pos];
	int b = lat->bitTable[0][vn];



	if ( vn == tad1_pos )
		dn1=0;
	if ( vn == tad2_pos )
		dn2=0;
	if ( vn == tad3_pos )
		dn3=0;

	for ( int v = 0; (v < 12) ; ++v )
	{
		if ( lat->bitTable[v+1][vn] == tad1_pos )
			dn1=lat->opp[v+1];
		if ( lat->bitTable[v+1][vn] == tad2_pos )
			dn2=v+1;
		if ( lat->bitTable[v+1][vn] == tad3_pos )
			dn3=v+1;
	}
	
	//all three direction have been found compatible: standard move
	if( dn1!=-1 and dn2!=-1 and dn3!=-1 )
	{
		// Check if new position vn complies with occupancy criteria and bonds
		legal = ((b == 0) || ( (b == 1) && ( (vn == tad1_pos) || (vn == tad2_pos) || (vn == tad3_pos) ) )|| reptation_values.size()>0);
		
		if(legal)
		{
			auto fork_standard_move_values = {tad->pos, vn, dn1,dn2,dn3,rept_dir};
			reptation_values.push_back(fork_standard_move_values);
			reptating_tads.push_back(tad);
		}
		return;
	}
	
	//compatible move not found: at least two direction have not been updated
	if( ( dn1== -1 && dn2== -1  ) or (dn1== -1 && dn3== -1 ) or (dn3== -1 && dn2== -1 ) )
	{
		reptation_values.clear();
		return;
	}
	//in presence of reptation (b=1) moving collectively is alway forbidden
	if ((tad->pos==tad1_pos or tad->pos==tad2_pos or tad->pos==tad3_pos) and  reptation_values.size()==0)
	{
		reptation_values.clear();
		return;
	}
	
	//find the reptation direction
	if((b == 0) || ( (b == 1) && ( (vn == tad1_pos) || (vn == tad2_pos) || (vn == tad->pos) ) ) || reptation_values.size()>0)
	{
		if(dn1 == -1 )
		{
			reptating_tads.push_back(tad);
			reptating_tads.push_back(tad->neighbors[1]);
			reptating_tads.push_back(tad->neighbors[2]);


			reptation_values.push_back({tad->pos, vn, lat->opp[dir],dn2,dn3,rept_dir});
			TrialReptationMove(tad,0);
		}
		else if (dn2 == -1)
		{
			reptating_tads.push_back(tad);
			reptating_tads.push_back(tad->neighbors[2]);
			reptating_tads.push_back(tad->neighbors[0]);
			reptation_values.push_back({tad->pos, vn, dn1,dir,dn3,rept_dir});
			TrialReptationMove(tad,1);

		}
		else if (dn3 == -1)
		{
			reptating_tads.push_back(tad);
			reptating_tads.push_back(tad->neighbors[1]);
			reptating_tads.push_back(tad->neighbors[0]);
			reptation_values.push_back({tad->pos, vn, dn1,dn2,dir,rept_dir});
			TrialReptationMove(tad,2);

		}

	}
}

void MCTadUpdater::TrialReptationMove(const MCTad* tad, int dir) 
{
	//iterative process: I will move along the chain until I manage to obtain a legal configuration

	//I cannot try to move again one reptating tad or its neighbourg
	for ( int i = 0; i < (int) reptating_tads.size(); ++i )
		if(tad->neighbors[dir]==reptating_tads[i])
			return;
	
	//if next monomer is a terminal one always conclude succesfully the reptation
	if(tad->neighbors[dir]->isLeftEnd() or tad->neighbors[dir]->isRightEnd())
	{
		legal=true;
		int vo= tad->neighbors[dir]->pos;
		int vn = reptation_values.back()[0];
		
		int dn1 = dir==0 ? reptation_values.back()[2] : reptation_values.back()[3] ;
		auto end_monomer={vo, vn ,dn1, dir};
		reptation_values.push_back(end_monomer);
		reptating_tads.push_back(tad->neighbors[dir]);
		return;
	}
	else if(tad->neighbors[dir]->isChoesin or tad->neighbors[dir]->isFork())
	{
		//Next tad is a Fork or Choesin must go back and check for connectivity consistency with at least 2 direction

		auto next_fork = tad->neighbors[dir];
		
		//find the direction the fork is trying to move to
		int third_dir=-1;
		
		if ( tad->pos == next_fork->pos )
			third_dir=0;
		
		for ( int v = 0; (v < 12) ; ++v )
		{
			if ( lat->bitTable[v+1][next_fork->pos] == tad->pos )
				third_dir=v+1;
		}
		
		// one of the three position is the one left by the previous reptating monomer -> reptation_values.back()[1]
		
		if(next_fork->neighbors[0]==tad)
			CheckForkLegal(next_fork,reptation_values.back()[1], next_fork->neighbors[1]->pos, next_fork->neighbors[2]->pos, third_dir, dir);
		else if(next_fork->neighbors[1]==tad)
			CheckForkLegal(next_fork,next_fork->neighbors[0]->pos,reptation_values.back()[1], next_fork->neighbors[2]->pos, third_dir, dir);
		else if (next_fork->neighbors[2]==tad)
			CheckForkLegal(next_fork, next_fork->neighbors[0]->pos, next_fork->neighbors[1]->pos,reptation_values.back()[1], third_dir, dir);

		return;
	}
	else
	{

		while(legal==false)
		{
			//find orientation for next step of loop when I find a fork.
			//here the next monomer is ALWAYS part of a linear connection so I have to determine if the direction of reptation is either 0 or 1. To determine the direction I find the one which does not bring me back to the fork or choesin
			int nextdir;
			if(tad->isFork() or tad->isChoesin )
				nextdir = (tad->neighbors[dir]->neighbors[0]==tad) ? 1 :0 ;
			else
				nextdir=dir;
			
			//tad that is reptating
			MCTad* reptad = tad->neighbors[dir];
			int vo = reptad->pos;

			//consecutive tad in the reptating direction necessary to check for connectivity
			MCTad* tad1 = reptad->neighbors[nextdir];


			for ( int i = 0; i < (int) reptating_tads.size(); ++i )
				if(reptad==reptating_tads[i])
					return;
			
			int vn = tad->pos;
			int dn1=-1;
			int dn2=-1;

			
			//Connectivity check
			for ( int v = 0; (v < 12) ; ++v )
			{

				if ( lat->bitTable[v+1][vn] == tad1->pos )
				{
					dn1 = nextdir ==0 ? lat->opp[v+1] : tad->bonds[dir]->dir;
					dn2 = nextdir ==0 ? tad->bonds[dir]->dir : lat->opp[v+1];
					legal = true;
				}

			}
			

			if(!legal)
			{
				// not legal I need to load the reptation_values vector befor moving to the next monomer: check for the special cases
				
				if(!tad->isLeftEnd() and !tad->isRightEnd())
				{
					dn1 = nextdir ==0 ? lat->opp[tad->bonds[dir]->dir] : reptation_values.back()[3];
					dn2 = nextdir ==0 ? lat->opp[reptation_values.back()[2]] : tad->bonds[dir]->dir;
				}
				else
				{
					dn1 = nextdir ==0 ? lat->opp[tad->bonds[dir]->dir] : reptation_values.back()[2];
					dn2 = nextdir ==0 ? lat->opp[reptation_values.back()[2]] : tad->bonds[dir]->dir;
				}
		
				if(reptad->neighbors[nextdir]->isLeftEnd() or reptad->neighbors[nextdir]->isRightEnd())
				{
					//Succesfully conclude reptation
					reptation_values.push_back({vo, vn, dn1, dn2, nextdir});
					reptating_tads.push_back(reptad);

					legal=true;
					vo= reptad->neighbors[nextdir]->pos;
					vn = reptation_values.back()[0];
					
					dn1 = dir==0 ? reptation_values.back()[2] : reptation_values.back()[3] ;
					auto end_monomer={vo, vn ,dn1,nextdir};
					reptation_values.push_back(end_monomer);
					reptating_tads.push_back(reptad->neighbors[nextdir]);
					return;
					
				}
				else if(reptad->neighbors[nextdir]->isChoesin or reptad->neighbors[nextdir]->isFork())
				{
					//As above:
					//Next tad is a Fork or Choesin must go back and check for connectivity consistency with at least 2 direction

					for ( int i = 0; i < (int) reptating_tads.size(); ++i )
						if(reptad->neighbors[nextdir]==reptating_tads[i])
							return;
					
					reptation_values.push_back({vo, vn, dn1, dn2, dir});
					reptating_tads.push_back(reptad);
					auto next_fork = reptad->neighbors[nextdir];

					//iteratevely find a new reptation moves compatible with the fork. The standard CheckForkLegal function is used but it is needed to find which 2 out of the 3 monomer involved I am moving
					int third_dir=-1;
					if ( reptad->pos == next_fork->pos )
						third_dir=0;
					for ( int v = 0; (v < 12) ; ++v )
						if ( lat->bitTable[v+1][next_fork->pos] == reptad->pos )
							third_dir=v+1;
					

					if(next_fork->neighbors[0]==reptad)
						CheckForkLegal(next_fork,reptation_values.back()[1],next_fork->neighbors[1]->pos,next_fork->neighbors[2]->pos,third_dir,nextdir);
					else if(next_fork->neighbors[1]==reptad)
						CheckForkLegal(next_fork,next_fork->neighbors[0]->pos,reptation_values.back()[1],next_fork->neighbors[2]->pos,third_dir,nextdir);
					else if (next_fork->neighbors[2]==reptad)
						CheckForkLegal(next_fork, next_fork->neighbors[0]->pos,next_fork->neighbors[1]->pos,reptation_values.back()[1],third_dir,nextdir);

					return;
				}
			}
			else
			{
				//move is legal! update with either 2, or 3 direction (left/right end cannot enter to the stage of the loop, return before)
				//NB vo and vn are the same since I found the legal move
				reptation_values.push_back({vo, vn, dn1, dn2, dir});
				reptating_tads.push_back(reptad);

				//push the fork in the vector
				if(tad1->isFork())
				{
					vo= tad1->pos;
					vn = tad1->pos;
					auto next_fork = reptad->neighbors[nextdir];
					int third_dir = nextdir==0 ? lat->opp[dn2] : dn1 ;

					if(next_fork->neighbors[0]==reptad)
						reptation_values.push_back({vo, vn, third_dir, tad1->bonds[1]->dir,tad1->bonds[2]->dir ,nextdir});

					else if(next_fork->neighbors[1]==reptad)
						reptation_values.push_back({vo, vn, tad1->bonds[0]->dir, third_dir, tad1->bonds[2]->dir, nextdir});
					else
						reptation_values.push_back({vo, vn, tad1->bonds[0]->dir,tad1->bonds[1]->dir, third_dir, nextdir});
					
					reptating_tads.push_back(reptad->neighbors[nextdir]);
					return;
				}
				else
				{
					//push linear monomer in the vector
					vo= tad1->pos;
					vn = tad1->pos;
					dn1 = nextdir==0 ? tad1->bonds[1]->dir : reptation_values.back()[3];
					dn2 = nextdir==0 ? reptation_values.back()[2]: tad1->bonds[0]->dir;
					auto last_monomer={vo, vn ,dn1,dn2,nextdir};
					reptation_values.push_back(last_monomer);
					reptating_tads.push_back(reptad->neighbors[nextdir]);
					return;
				}
			}
			//If function didn't return I'm still in a linear part of the chain: adding to the reptating_values the rept_tad information
			reptation_values.push_back({vo, vn, dn1, dn2, dir});
			reptating_tads.push_back(tad->neighbors[dir]);

			//update for next step of loop
			tad=tad->neighbors[dir];
			dir=nextdir;



		}
	}

}

void MCTadUpdater::AcceptMove(MCTad* tad) const
{
	//accept move by updatind the direction. I need to recontruct the reptating portion of the chain by following the direction given by the last element of the reptation_values vector

	for ( int i = 0; i < (int) reptation_values.size(); ++i )
	{

		if(i != 0)
			tad=tad->neighbors[reptation_values[i].back()];
		
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
