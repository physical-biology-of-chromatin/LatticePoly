//
//  MCReplicPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#include "MCReplicPoly.hpp"

#include <iterator>
#include <algorithm>


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	
	MCsteps=0;
	MCrepl=0;
	activeForks.reserve(Nchain);

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	/*
	for ( int i=Nchain/2 ; i < Nchain/2 +100 ; i++)
	{
		while(tadConf.at(i).status==0){
			Replicate(&tadConf.at(i));
		}
	}
	std::vector<int> originsvector;
	for (int i=0; i<Nchain; ++i) originsvector.push_back(i);
	std::random_shuffle ( originsvector.begin(), originsvector.end() );
	for ( int i=0 ; i < 93; i++)
	{
		origins.push_back(originsvector[i]);
	
	}
*/
	origins={0,7,12,17,36,40,68,76,80,99,110,126,170,185,188,203,205,253,263,281,326,348,355,370,381,387,
		404,444,454,503,511,515,562,576,598,602,644,676,687,703,719,723,731,737,754,780,813,817,
		826,846,888,904,927,932,963,992,1021,1042,1046,1082,1103,1108,1123,1158,1163,1169,1189,1199,1202,1204,1219};
	
 mrt={1000,1000,1000,0.9208379238843918,0.7043074965476991,0.7438885271549225,0.7962751537561417,0.8440044820308685,0.8253782838582993,0.4947620034217834,0.642608255147934,0.8812578544020653,0.5261937975883484,0.6554138362407684,0.6670551002025604,0.7322472929954529,0.6728757321834564,0.3876601457595825,0.20954644680023196,0.6647264659404755,0.34575146436691284,0.2060534954071045,0.3015140891075134,0.16298043727874756,0.32828861474990845,0.4039586782455444,0.4051220417022705,0.07683342695236206,0.2630963921546936,0.7741559892892838,0.7334106266498566,0.7415598928928375,0.6938305497169495,0.8661236464977264,0.6821883320808411,0.6880099475383759,0.7229337096214294,0.9208379238843918,0.850989431142807,1000,0.22118771076202395,0.1548311710357666,0.0,0.09313195943832396,0.637950986623764,0.9289871901273729,0.3725259304046631,0.3597203493118286,0.6123398542404175,0.551804929971695,0.7881258875131607,0.6682194173336029,0.0861470103263855,0.18509960174560547,0.8672879636287689,0.6949939131736755,0.7660067230463028,0.5506406128406525,0.7834695726633072,0.4761348366737366,0.8416768163442612,0.7415598928928375,0.7334106266498566,0.7462162077426909,0.6821883320808411,0.5250294804573059,0.8125727027654648,0.850989431142807,0.8195576518774033,0.8428401648998259,0.8800935298204422};
	

	std::cout << "Origins :";
	for (std::vector<int>::iterator it=origins.begin(); it!=origins.end(); ++it)
	std::cout << ' ' << *it;


}

void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
}

void MCReplicPoly::OriginMove()
{
	if ( origins.size() > 0 and MCsteps> (Nrelax)*Ninter )
	{

		auto originsCopy =origins;
		auto mrtCopy =mrt;

		std::vector<int> indexes; //create a indexes vector
		indexes.reserve(originsCopy.size());
		for (int i = 0; i < (int)originsCopy.size(); ++i)
			indexes.push_back(i); //populate
		std::random_shuffle(indexes.begin(), indexes.end()); // randomize
		
		for ( int i=0 ; i < (int)indexes.size(); i++) //for every element in indexes
		{
			MCTad* origin = &tadConf[originsCopy[indexes[i]]]; //select origin taf
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < (11- int(Nfork/2) + 0.5)*originRate*exp(-0*mrtCopy[indexes[i]]) and origin->status==0)
			{
				Replicate(origin);
				std::vector<int>::iterator itr = std::find(origins.begin(), origins.end(), originsCopy[indexes[i]]);

				origins.erase(origins.begin()+std::distance(origins.begin(), itr));
				mrt.erase(mrt.begin()+ std::distance(origins.begin(), itr));
				
			}
		}
	}
	MCsteps+=1;
}
void MCReplicPoly::ForkMove()
{
	if ( Nfork > 0)
	{
		auto activeForksCopy =activeForks;
		for ( int i=0 ; i < (int)activeForksCopy.size(); i++)
		{
			MCTad* fork = activeForks[i];
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < replicRate and fork->isFork())
			{
				Replicate(fork);
			}
		}
	}
}
void MCReplicPoly::Replicate(MCTad* tad)
{
	if (tad->isRightEnd() || tad->isLeftEnd())
		return;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
		
	double rnd = lat->rngDistrib(lat->rngEngine);
	
	if ( !tad->isFork() )
	{
		// Can't replicate tad if it's already adjacent to a fork
		if ( nb1->isFork() || nb2->isFork() )
			return;
		
		else
		{
			// Replicate extremities at half the normal rate
			if ( nb1->isLeftEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}

			else
			{
				activeForks.push_back(nb1);
				UpdateReplTable(nb1);
			}
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
			{
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);
			}
		}
	}
	
	else
	{
		// Replicating left fork means displacing it to its left neighbor, so we need to check if it's already a fork or chain end
		if ( tad->isLeftFork() )
		{
			if ( nb1->isLeftFork() || nb1->isRightEnd() )
				// Probably should never happen, do nothing
				return;
			
			if ( nb1->isRightFork() || nb1->isLeftEnd() )
			{
				// Merge forks/replicate extremities at half the normal rate
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb1);
			}
		}
		
		// Same for right forks
		else if ( tad->isRightFork() )
		{
			if ( nb2->isRightFork() || nb2->isLeftEnd() )
				return;
			
			if  ( nb2->isLeftFork() || nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
		
			else
			{
				activeForks.push_back(nb2);
			}
		}
		
		// Delete old forks
		auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
		activeForks.erase(fork);
		
		if ( nb1->isFork() || nb2->isFork() )
		{
			auto fork2 = std::find(activeForks.begin(), activeForks.end(), nb1->isFork() ? nb1 : nb2);
			activeForks.erase(fork2);
		}
	}
	
	ReplicateTADs(tad);
	ReplicateBonds(tad);

	Update();
	UpdateReplTable(tad);

	
}

void MCReplicPoly::ReplicateTADs(MCTad* tad)
{
	
	MCTad tadReplic;
	
	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];

	// Replicate left end/fork, if applicable
	if ( nb1->isLeftEnd() || nb1->isRightFork() )
	{
		tadReplic = *nb1;
		tadConf.push_back(tadReplic);
		nb1->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb1);


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);

	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		tadConf.push_back(tadReplic);
		nb2->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb2);


	}
}

void MCReplicPoly::ReplicateBonds(MCTad* tad)
{
	MCBond bondReplic1;
	MCBond bondReplic2;

	MCTad* nb1 = tad->neighbors[0];
	MCTad* nb2 = tad->neighbors[1];
	
	// Create/modify bonds between relevant neighbors and replicated tad
	MCBond* bond1 = tad->isRightFork() ? tad->bonds[2] : tad->bonds[0];
	MCBond* bond2 = tad->isLeftFork() ? tad->bonds[2] : tad->bonds[1];
	
	bondReplic1.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad : bond1->id1;
	bondReplic1.id2 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
				
	bondReplic1.dir = bond1->dir;

	bondReplic2.id1 = (nb1->isLeftEnd() || nb1->isRightFork()) ? Ntad+1 : Ntad;
	bondReplic2.id2 = (nb2->isRightEnd() || nb2->isLeftFork()) ? Ntad+1 : bond2->id2;
		
	bondReplic2.dir = bond2->dir;
	
	// For right forks, update bond1 to link left (replicated) neighbor to new tad
	if ( tad->isRightFork() )
	{
		bond1->id2 = bondReplic1.id2;
		
		CreateBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			CreateBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		CreateBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			CreateBond(*bond3);
			UnsetFork(nb1);
		}
	}
	
	else
		tadTopo.push_back(bondReplic2);
}

void MCReplicPoly::UnsetFork(MCTad* tad)
{
	if ( tad->isFork() )
	{
		tad->bonds[2] = nullptr;
		tad->neighbors[2] = nullptr;
		
		--tad->links;
	}
}

void MCReplicPoly::Update()
{
	// Update bonds
	if ( (int) tadTopo.size() > Nbond )
	{
		for ( auto bond = tadTopo.begin()+Nbond; bond != tadTopo.end(); ++bond )
			CreateBond(*bond);
		
		Nbond = (int) tadTopo.size();
	}
	
	// Update tads
	if ( (int) tadConf.size() > Ntad )
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
	
	Nfork = (int) activeForks.size();
}

double MCReplicPoly::GetEffectiveEnergy() const
{
	if ( Jf > 0.  )
	{
		
		if ( tadTrial->isFork()){
			return 	MCHeteroPoly::GetEffectiveEnergy() +Jf * (ReplTable[0][tadUpdater->vo]-ReplTable[0][tadUpdater->vn]);
		}
		else{
			return 	MCHeteroPoly::GetEffectiveEnergy();
		}
	}
	if ( Jpair > 0.  )
	{
		if (tadTrial->status == 1)
		{
			double Jbott1=0.0;
			double Jbott2=0.0;

			int SistPos = tadConf.at(tadTrial->SisterID).pos;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				if(SistPos==vi2 and tadTrial->SisterID%1==0){
					Jbott2=Jpair;
				}
				if(SistPos==vi1 and tadTrial->SisterID%1==0){
					Jbott1=Jpair;
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		if (tadTrial->status == -1)
		{
			double Jbott1=0.0;
			double Jbott2=0.0;

			int SistPos = tadConf.at(tadTrial->SisterID).pos;
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				if(SistPos==vi2 and tadConf.at(tadTrial->SisterID).SisterID%1==0){
					Jbott2=Jpair;
				}
				if(SistPos==vi1 and tadConf.at(tadTrial->SisterID).SisterID%1==0){
					Jbott1=Jpair;
				}
			}


			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
	}
	return MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->isFork()) //increase energy at fork site
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			--ReplTable[0][vi1];
			++ReplTable[0][vi2];
		}
	}
	
	if( tadTrial->isLeftEnd()==false and tadTrial->isRightEnd()==false) //increase energy at fork's neighbouring sites,first check if terminal monomers to avoid segmentation errors
	{
		if ( tadTrial->neighbors[0]->isFork() or tadTrial->neighbors[1]->isFork())
		{
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				
				--ReplTable[0][vi1];
				++ReplTable[0][vi2];
			}
		}
	}else if(tadTrial->isLeftEnd())
		{
		if ( tadTrial->neighbors[1]->isFork())
		{
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				
				--ReplTable[0][vi1];
				++ReplTable[0][vi2];
			}
		}
		}else if(tadTrial->isRightEnd())
		{
			if ( tadTrial->neighbors[0]->isFork())
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
					int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
					
					--ReplTable[0][vi1];
					++ReplTable[0][vi2];
				}
			}
		}

	
	if ( tadTrial->status == -1)
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];

			--ReplTable[1][vi1];
			++ReplTable[1][vi2];
			
		}
	}
	if ( tadTrial->status == 1)
	{
	
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			
			--ReplTable[2][vi1];
			++ReplTable[2][vi2];
		}
	}
}

void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
	
	if(tad->status == -1)
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[1][vi];
		}
	}
	if(tad->status == 1)
	{
	for ( int v = 0; v < 13; ++v )
		{

			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[2][vi];
		}
	}
	
	if(tad->isRightEnd() || tad->isLeftEnd())
	{
	for ( int v = 0; v < 13; ++v )
	{
		int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
	
		--ReplTable[0][vi];
		}
	}
	else{
		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
			
			++ReplTable[0][vi];
		}
	}
	 
}
std::vector<double3> MCReplicPoly::GetPBCConf()
{
	std::vector<MCTad*> leftEnds;
	std::vector<MCTad*> builtTads;
	
	std::vector<double3> conf(Ntad);
	
	builtTads.reserve(Ntad);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			conf[t][i] = lat->xyzTable[i][tadConf[t].pos];
		
		if ( tadConf[t].isLeftEnd() )
			leftEnds.push_back(&tadConf[t]);
	}
	
	// Grow chains recursively, starting from their respective left extremities
	auto leftEnd = leftEnds.begin();
	
	while ( (int) builtTads.size() < Ntad )
	{
		MCTad *tad1, *tad2;
		tad1 = *leftEnd;
		
		bool builtTad1 = (std::find(builtTads.begin(), builtTads.end(), tad1) != builtTads.end());

		if ( !builtTad1 )
		{
			builtTads.push_back(tad1);

			// Traverse main branch
			while ( (tad2 = tad1->neighbors[1]) )
			{
				bool builtTad2 = (std::find(builtTads.begin(), builtTads.end(), tad2) != builtTads.end());

				if ( !builtTad2 )
				{
					BuildPBCPair(builtTads, conf, tad1, tad2);
					
					// Traverse side branches
					if ( tad2->isFork() )
					{
						MCTad *tad3, *tad4;
						tad3 = tad2->neighbors[2];
						
						BuildPBCPair(builtTads, conf, tad2, tad3);
					
						while ( (tad4 = (tad2->isLeftFork() ? tad3->neighbors[1] : tad3->neighbors[0])) )
						{
							BuildPBCPair(builtTads, conf, tad3, tad4);
						
							if ( tad4->isFork() )
								break;
							
							tad3 = tad4;
						}
					}
				}
				
				tad1 = tad2;
			}
		}
		
		++leftEnd;
	}
	
	int chainNum = Ntad / Nchain;
	int chainLength = (chainNum == 1) ? Ntad : Nchain;
	
	std::vector<double3> centers(chainNum);
	
	for ( int c = 0; c < chainNum; ++c )
	{
		auto end1 = conf.begin() + c*chainLength;
		auto end2 = conf.begin() + (c+1)*chainLength;

		centers[c] = GetPBCCenterMass(end1, end2);
	}
	
	centerMass = {0., 0., 0.};
	
	for ( int c = 0; c < chainNum; ++c )
	{
		for ( int i = 0; i < 3; ++i )
			centerMass[i] += centers[c][i] / ((double) chainNum);
	}
	
	return conf;
}

void MCReplicPoly::BuildPBCPair(std::vector<MCTad*>& builtTads, std::vector<double3>& conf, MCTad* tad1, MCTad* tad2)
{
	int id1 = (int) std::distance(tadConf.data(), tad1);
	int id2 = (int) std::distance(tadConf.data(), tad2);
	
	FixPBCPair(conf, id1, id2);
	
	builtTads.push_back(tad2);
}
