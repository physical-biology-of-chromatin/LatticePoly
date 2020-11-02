//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>

#include "MCReplicPoly.hpp"


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	Replicate(Nchain/2, Nchain/2 +2);
	
	Update();

	CreateOrigins();

	std::cout << "finish init with "<<activeforks.size()<<  std::endl;

}

void MCReplicPoly::CreateOrigins()
{
	if(RandomOrigin==true)
	{
		std::vector<int> arr;
		for(int i = 2; i < Nchain-3; ++i)
			arr.push_back(i);
		
		std::random_shuffle ( arr.begin(), arr.end() );
		for(int i = 0; i < Norigins; ++i)
			origins.push_back(arr[i]);
	}
}

void MCReplicPoly::Replicate(int tOrig, int tEnd)
{
	auto origin = tadConf.begin() + tOrig;
	auto end = tadConf.begin() + tEnd;
	
	if ( !origin->isLeftEnd())
	{
		activeforks.push_back(tOrig);
		origin->replstatus=-1;
	}
	
	if ( !end->isRightEnd() )
	{
		activeforks.push_back(tEnd);
		end->replstatus=+1;

	}
	if ( !origin->isFork())
	{

		auto nextFork = std::find_if(origin, end+1, [](const MCTad& t){return t.isFork();});
		end = nextFork-1;
		
		int dist = (int) std::distance(origin, end);
		
		if ( dist > 1 )
		{
			ReplicateTADs(origin, end);
			ReplicateBonds(origin, end);
		}
	}
	
}

void MCReplicPoly::ReplicateTADs(std::vector<MCTad>::iterator origin, std::vector<MCTad>::iterator end)
{
	MCTad tadReplic;
	

	
	if ( origin->isLeftEnd() )
	{
		tadReplic = *origin;
		tadReplic.replstatus=2;
		tadReplic.sisterID= 0;
		++ReplTable[0][origin->pos];
		++ReplTable[1][tadReplic.pos];

		
		tadConf.push_back(tadReplic);
	}
	
	for ( auto tad = origin+1; tad != end; ++tad )
	{
		tad->replstatus=-2;
		tadReplic = *tad;
		tadReplic.replstatus=2;
		tadReplic.sisterID= tad->bonds[0]->id2;
		
		++ReplTable[0][tad->pos];
		++ReplTable[1][tadReplic.pos];
		
		tadConf.push_back(tadReplic);

	}
	
	if ( end->isRightEnd())
	{
		end->replstatus=2;
		tadReplic.replstatus=2;
		tadReplic.sisterID= Ntad-1;
		
		++ReplTable[0][end->pos];
		++ReplTable[1][tadReplic.pos];

		
		tadConf.push_back(tadReplic);
	}
}

void MCReplicPoly::ReplicateBonds(std::vector<MCTad>::iterator origin, std::vector<MCTad>::iterator end)
{
	MCBond bondReplic;

	int Nreplic = (int) tadConf.size() - Ntad;
	
	int tOrig = (int) std::distance(tadConf.begin(), origin);
	int tEnd = (int) std::distance(tadConf.begin(), end);
	
	MCBond* bond = origin->isLeftEnd() ? origin->bonds[0] : origin->bonds[1];

	if ( !origin->isLeftEnd() )
	{
		bondReplic.id1 = tOrig;
		bondReplic.id2 = Ntad;
				
		bondReplic.dir = bond->dir;
	
		tadTopo.push_back(bondReplic);
		
		++bond;
	}
				
	for ( int t = 0; t < Nreplic-1; ++t )
	{
		bondReplic.id1 = t+Ntad;
		bondReplic.id2 = t+Ntad+1;
		
		bondReplic.dir = bond->dir;
		
		tadTopo.push_back(bondReplic);

		++bond;
	}
	
	if ( !end->isRightEnd() )
	{
		bondReplic.id1 = Ntad+Nreplic-1;
		bondReplic.id2 = tEnd;
		
		bondReplic.dir = bond->dir;

		tadTopo.push_back(bondReplic);
		
		++bond;
	}
	

}

void MCReplicPoly::Update()
{

	if ( (int) tadTopo.size() > 0 )
	{
		for ( auto bond = tadTopo.begin()+Nbond; bond != tadTopo.end(); ++bond )
		{
			CreateBond(*bond);
			
		}
		Nbond = (int) tadTopo.size();
	}
	
	if ( (int) tadConf.size() > 0 )
	{
		for ( auto tad = tadConf.begin()+Ntad; tad != tadConf.end(); ++tad )
		{

			if ( tad->type == 1 )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
					++hetTable[vi];
				}
			}
			
			++lat->bitTable[0][tad->pos];
		}
		
		Ntad = (int) tadConf.size();
		
	}
	
	
	/*for ( int k = 0; k < Ntad+1; ++k ){
		auto monomer = tadConf.begin() + k;

		for ( int j = 0; j < monomer->links;++j ){
			for ( int  s= 0; s < monomer->neighbors[j]->links; ++s ){
				std::cout << "monomer with links "<< monomer->links<<" at pos " << k << " has a neig that connect "<< monomer->neighbors[j]->bonds[s]->id1<<" and "<< monomer->neighbors[j]->bonds[s]->id2<<" with direction "<< monomer->neighbors[j]->bonds[s]->dir<< std::endl;
			}
			std::cout << "END MONOMER"<< std::endl;

		}
	}*/
}


void MCReplicPoly::MoveFork(int forkID, int i)
{
	MCTad tadReplic;
	MCTad tadReplic2;
	MCBond bondReplic;
	
	auto fork = tadConf.begin() + forkID;

	if(fork->replstatus < 0){
		if (forkID!=1)
		{
			int previousmonomerID =fork->bonds[0]->id1 ;
			auto previousmonomer = tadConf.begin() + previousmonomerID;
			int nextmonomerID = fork->bonds[2]->id2;
			auto nextmonomer = tadConf.begin() + nextmonomerID;
			if(previousmonomer->replstatus==0)
			{
				tadReplic = *fork;

				tadConf.push_back(tadReplic);
				fork->links=2;

				nextmonomer->neighbors[0]=&tadConf.at(Ntad);
				nextmonomer->bonds[0]->dir = fork->bonds[2]->dir;
				nextmonomer->bonds[0]->id1 = Ntad;
				
				auto newmonomer = tadConf.begin() + Ntad;
				newmonomer->links=1;
				newmonomer->neighbors[0]=fork->neighbors[2];
				newmonomer->bonds[0]=fork->bonds[2];
				newmonomer->bonds[0]->id1=Ntad;
				newmonomer->bonds[0]->dir=fork->bonds[2]->dir;
				newmonomer->replstatus=2;
				newmonomer->sisterID= fork->bonds[0]->id2;
				

				previousmonomer->neighbors[2] = &tadReplic;
				bondReplic.id1=previousmonomerID;
				bondReplic.id2= Ntad;
				bondReplic.dir=fork->bonds[0]->dir;
				tadTopo.push_back(bondReplic);
				
				Update();
				

				activeforks[i]=forkID-1;
				fork->replstatus=-2;
				previousmonomer->replstatus= -1;
				std::swap(tadConf.back().bonds[0],tadConf.back().bonds[1]);
				std::swap(tadConf.back().neighbors[0],tadConf.back().neighbors[1]);

				++ReplTable[0][fork->pos];
				++ReplTable[1][newmonomer->pos];

				/*for ( int k = 0; k < Ntad+1; ++k ){
					auto monomer = tadConf.begin() + k;
					
					for ( int j = 0; j < monomer->links;++j ){
						for ( int  s= 0; s < monomer->neighbors[j]->links; ++s ){
							std::cout << "monomer with links "<< monomer->links<<" at pos " << k << " has a neig that connect "<< monomer->neighbors[j]->bonds[s]->id1<<" and "<< monomer->neighbors[j]->bonds[s]->id2<<" with direction "<< monomer->neighbors[j]->bonds[s]->dir<< std::endl;
						}
						std::cout << "END MONOMER"<< std::endl;
					}
				}*/
			}else{
				auto previousfork = find(activeforks.begin(),activeforks.end(), previousmonomer->bonds[1]->id1);
				int previousfork_index = (int) std::distance(activeforks.begin(), previousfork);
				MoveFork(previousmonomer->bonds[1]->id1, previousfork_index);
			}
		}
	}
	if(fork->replstatus > 0){
		if (forkID != Nchain-2)
		{
			int previousmonomerID = fork->bonds[2]->id1;
			auto previousmonomer = tadConf.begin() + previousmonomerID;
			int nextmonomerID =fork->bonds[1]->id2 ;
			auto nextmonomer = tadConf.begin() + nextmonomerID;
			if(nextmonomer->replstatus==0)
			{
				tadReplic = *fork;


				tadConf.push_back(tadReplic);

				
				previousmonomer->bonds[1]->dir = fork->bonds[2]->dir;
				previousmonomer->bonds[1]->id2 = Ntad;
				previousmonomer->neighbors[1]=&tadConf.at(Ntad);
				
				auto newmonomer = tadConf.begin() + Ntad;
				newmonomer->links=1;
				newmonomer->neighbors[0]=fork->neighbors[2];
				newmonomer->bonds[0]=fork->bonds[2];
				newmonomer->bonds[0]->id2=Ntad;
				newmonomer->bonds[0]->dir=fork->bonds[2]->dir;
				newmonomer->replstatus=2;
				newmonomer->sisterID= fork->bonds[0]->id2;

				
				
				
				nextmonomer->neighbors[2] = &tadReplic;
				bondReplic.id2=nextmonomerID;
				bondReplic.id1= Ntad;
				bondReplic.dir=fork->bonds[1]->dir;
				tadTopo.push_back(bondReplic);
				
				
				Update();
				
				fork->links=2;
				activeforks[i]=forkID+1;
				fork->replstatus=-2;
				nextmonomer->replstatus= +1;
				
				++ReplTable[0][fork->pos];
				++ReplTable[1][newmonomer->pos];
				
				/*for ( int k = 0; k < Ntad+1; ++k ){
					auto monomer = tadConf.begin() + k;
					
					for ( int j = 0; j < monomer->links;++j ){
						for ( int  s= 0; s < monomer->neighbors[j]->links; ++s ){
							std::cout << "monomer with links "<< monomer->links<<" at pos " << k << " has a neig that connect "<< monomer->neighbors[j]->bonds[s]->id1<<" and "<< monomer->neighbors[j]->bonds[s]->id2<<" with direction "<< monomer->neighbors[j]->bonds[s]->dir<< std::endl;
						}
						std::cout << "END MONOMER"<< std::endl;
					}
				}*/

			}else{
				std::cout << "FORK MERGING +1"<< std::endl;

				double rnd = lat->rngDistrib(lat->rngEngine);
				if(rnd<0.5){
					


					tadReplic = *fork;
					tadConf.push_back(tadReplic);
		
					previousmonomer->neighbors[1]= &tadConf.at(Ntad);
					previousmonomer->bonds[1]->id2=Ntad;
					
					auto newmonomer1=tadConf.begin() + Ntad;
					newmonomer1->links=1;
					newmonomer1->neighbors[0]=fork->neighbors[2];
					newmonomer1->bonds[0]=fork->bonds[2];
					newmonomer1->bonds[0]->dir=fork->bonds[2]->dir;
					newmonomer1->replstatus=2;
					newmonomer1->sisterID= fork->bonds[0]->id2;


					tadReplic2 = *nextmonomer;
					tadConf.push_back(tadReplic2);

					bondReplic.id1= Ntad;
					bondReplic.id2=Ntad+1;
					bondReplic.dir=fork->bonds[1]->dir;
					tadTopo.push_back(bondReplic);

					nextmonomer->neighbors[2]->bonds[0]->id1=Ntad+1;
					nextmonomer->neighbors[2]->neighbors[0]=&tadConf.at(Ntad+1);
					auto newmonomer2=tadConf.begin() + Ntad+1;
					
					newmonomer2->links=1;
					newmonomer2->neighbors[0]=nextmonomer->neighbors[2];
					newmonomer2->bonds[0]=nextmonomer->bonds[2];
					newmonomer2->bonds[0]->dir=nextmonomer->bonds[2]->dir;
					newmonomer2->replstatus=2;
					newmonomer2->sisterID= nextmonomer->bonds[0]->id2;


					Update();
					
					std::swap(tadConf.back().bonds[0],tadConf.back().bonds[1]);
					std::swap(tadConf.back().neighbors[0],tadConf.back().neighbors[1]);

					//togliere


					activeforks.erase(activeforks.begin()+i);
					auto nextmonomer_iterator = find(activeforks.begin(),activeforks.end(), nextmonomerID);
					int nextmonomer_index = (int) std::distance(activeforks.begin(), nextmonomer_iterator);
					activeforks.erase(activeforks.begin()+nextmonomer_index);

					fork->links=2;
					nextmonomer->links=2;
					fork->replstatus=-2;
					nextmonomer->replstatus=-2;

					++ReplTable[0][fork->pos];
					++ReplTable[0][nextmonomer->pos];
					++ReplTable[1][newmonomer1->pos];
					++ReplTable[1][newmonomer2->pos];
					
					/*for ( int k = 0; k < Ntad+1; ++k ){
						auto monomer = tadConf.begin() + k;
						
						for ( int j = 0; j < monomer->links;++j ){
							for ( int  s= 0; s < monomer->neighbors[j]->links; ++s ){
								std::cout << "monomer with links "<< monomer->links<<" at pos " << k << " has a neig that connect "<< monomer->neighbors[j]->bonds[s]->id1<<" and "<< monomer->neighbors[j]->bonds[s]->id2<<" with direction "<< monomer->neighbors[j]->bonds[s]->dir<< std::endl;
							}
							std::cout << "END MONOMER"<< std::endl;
						}
					}*/
				}
			}
		}
	}
}


void MCReplicPoly::CreateFork()
{
	if(origins.size()>0)
	{
		std::cout << "ERROR "<<  std::endl;

		int t1 = lat->rngEngine() % (int) origins.size();
		int neworigin = origins.at(t1);
		std::cout << neworigin <<  std::endl;

		if(tadConf.at(neworigin-1).replstatus==0 and tadConf.at(neworigin+1).replstatus==0)
		{
			
			Replicate(neworigin-1, neworigin+1);
			Update();
			origins.erase(origins.begin()+t1);


		}else{
			
			origins.erase(origins.begin()+t1);

		}

	}

}

double MCReplicPoly::GetEffectiveEnergy() const
{
	double dE=MCHeteroPoly::GetEffectiveEnergy();

	if ( Jpp > 0. )
	{
		if ( tadTrial->replstatus == -2){
			return 	dE +Jpp * (ReplTable[0][tadUpdater->vo]-ReplTable[0][tadUpdater->vn]);
		}
		
		else if ( tadTrial->replstatus == 2)
			return 	dE +Jpp * (ReplTable[1][tadUpdater->vo]-ReplTable[1][tadUpdater->vn]);
		
}	
	return 0.;
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();

	
	if ( tadTrial->replstatus==-2)
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			--ReplTable[0][vi1];
			++ReplTable[0][vi2];
		}
	}else{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
		
			--ReplTable[1][vi1];
			++ReplTable[1][vi2];
		}
	}
}
