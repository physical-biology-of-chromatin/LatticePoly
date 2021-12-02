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
#include <iostream>
#include <fstream>
#include <sstream>


MCReplicPoly::MCReplicPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCReplicPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	
	MCsteps=0;
	MCrepl=0;
	activeForks.reserve(Nchain);
	neigh=true;

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
	//chr4
	/*
	origins={0,7,12,17,36,40,68,76,80,99,110,126,170,185,188,203,205,253,263,281,326,348,355,370,381,387,
		404,444,454,503,511,515,562,576,598,602,644,676,687,703,719,723,731,737,754,780,813,817,
		826,846,888,904,927,932,963,992,1021,1042,1046,1082,1103,1108,1123,1158,1163,1169,1189,1199,1202,1204,1219};
	
 mrt={1000,1000,1000,0.9208379238843918,0.7043074965476991,0.7438885271549225,0.7962751537561417,0.8440044820308685,0.8253782838582993,0.4947620034217834,0.642608255147934,0.8812578544020653,0.5261937975883484,0.6554138362407684,0.6670551002025604,0.7322472929954529,0.6728757321834564,0.3876601457595825,0.20954644680023196,0.6647264659404755,0.34575146436691284,0.2060534954071045,0.3015140891075134,0.16298043727874756,0.32828861474990845,0.4039586782455444,0.4051220417022705,0.07683342695236206,0.2630963921546936,0.7741559892892838,0.7334106266498566,0.7415598928928375,0.6938305497169495,0.8661236464977264,0.6821883320808411,0.6880099475383759,0.7229337096214294,0.9208379238843918,0.850989431142807,1000,0.22118771076202395,0.1548311710357666,0.0,0.09313195943832396,0.637950986623764,0.9289871901273729,0.3725259304046631,0.3597203493118286,0.6123398542404175,0.551804929971695,0.7881258875131607,0.6682194173336029,0.0861470103263855,0.18509960174560547,0.8672879636287689,0.6949939131736755,0.7660067230463028,0.5506406128406525,0.7834695726633072,0.4761348366737366,0.8416768163442612,0.7415598928928375,0.7334106266498566,0.7462162077426909,0.6821883320808411,0.5250294804573059,0.8125727027654648,0.850989431142807,0.8195576518774033,0.8428401648998259,0.8800935298204422};
	*/
	
	//chr 7
	origins={1,6,14,20,26,51,55,89,91,94,130,133,137,149,163,185,192,213,215,220,228,230,255,270,282,311,336,365,385,388,407,431,442,454,459,460,465,473,485,497,501,523,527,572,594,612,622,636,64,667,677,690,706,710,733,776,782,799,801,804,829,850,859,866,871};
	
	weights={0.0008441 , 0.00113861, 0.00615295, 0.0011208 , 0.00126901,0.01631902, 0.00450673, 0.00759752, 0.00806823, 0.00269195,0.0199772 , 0.01153813, 0.002372  , 0.00223142, 0.03786356,0.00247823, 0.00170537, 0.00148846, 0.00186694, 0.00142549,0.04605837, 0.04407694, 0.00312641, 0.00132944, 0.00448002,0.0357918 , 0.07323744, 0.00183768, 0.05248039, 0.08644023,0.06210324, 0.00209721, 0.00130972, 0.02321365, 0.01614855,0.01131486, 0.00320274, 0.00343364, 0.00220153, 0.00144076,0.00138033, 0.00346672, 0.04149439, 0.04672691, 0.00159914,0.00221997, 0.06408213, 0.00273266, 0.00133516, 0.08104678,0.00311496, 0.00147574, 0.01425362, 0.07897057, 0.00378222,0.00140005, 0.00499844, 0.01212906, 0.01337517, 0.00358758,0.00097005, 0.00169265, 0.00099994, 0.00379431, 0.00139114};
	
		
	
	//origins={20,40,60,80,100,120,140,160,180};
	//origins={10,30,50,70,90,110,130,150,170,190,20,40,60,80,100,120,140,160,180};
	

	//mrt={0.8,0.6,0.4,0.19999999999999996,0.0,0.19999999999999996,0.3999999999999999,0.6000000000000001,0.8};
	//mrt={0.9,0.8,0.7,0.6,0.5,0.4,0.30000000000000004,0.19999999999999996,0.09999999999999998,0.0,0.10000000000000009,0.19999999999999996,0.30000000000000004,0.3999999999999999,0.5,0.600000000000001,0.7,0.8,0.8999999999999999};
	//CAR={5,15,25,35,45,65,75,85,95,105,125,135,155,175,185,195,10,30,50,70,90,110,130,150,170,190,20,40,60,80,100,120,140,160,180};
	
	for (int i = 0; i < (int)origins.size(); ++i)
		tadConf[origins[i]].isCAR=true;
	//origins={100};
		
}

void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	

}

void MCReplicPoly::OriginMove()
{
	
	if(Ntad>=int(.95*Nchain+Nchain))
	{
		
		std::ostringstream streamObj;
		streamObj << originRate;
		std::string strObj = streamObj.str();

		std::ofstream outfile("repltime"+std::to_string(Ndf)+ "_" + strObj, std::ios_base::app | std::ios_base::out);

		outfile << MCsteps << std::endl;


		
		
		exit(0);
		
	}

	if ( origins.size() > 0 and MCsteps> (Nrelax)*Ninter )
	{
		auto originsCopy =origins;
		auto weightsCopy =weights;

		std::vector<int> indexes; //create a indexes vector
		indexes.reserve(originsCopy.size());
		for (int i = 0; i < (int)originsCopy.size(); ++i)
			indexes.push_back(i); //populate
		std::random_shuffle(indexes.begin(), indexes.end()); // randomize
		
		for ( int i=0 ; i < (int)indexes.size(); i++) //for every element in indexes
		{
			MCTad* origin = &tadConf[originsCopy[indexes[i]]]; //select origin taf
			double rndReplic = lat->rngDistrib(lat->rngEngine);
			if ( rndReplic < (Ndf- int(double(Nfork)/2 + 0.5))*originRate*weightsCopy[indexes[i]] and origin->status==0)
			{

				
				
				Replicate(origin);

				
				
				std::vector<int>::iterator itr = std::find(origins.begin(), origins.end(), originsCopy[indexes[i]]);

				origins.erase(origins.begin()+std::distance(origins.begin(), itr));
				weights.erase(weights.begin()+ std::distance(origins.begin(), itr));
				
			}
		}
	}
	MCsteps+=1;
}
void MCReplicPoly::ForkMove()
{
	if ( Nfork > 0 )
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
	
	//origin replication
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
				UpdateReplTable(nb1);//increase energy around new fork
			}
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
			{
				activeForks.push_back(nb2);
				UpdateReplTable(nb2);//increase energy around new fork
			}
		}
	}
	
	else //fork displacement
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
				UpdateReplTable(nb1);//increase energy around new fork
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
				UpdateReplTable(nb2);//increase energy around new fork

			}
		}
		
		// Delete old forks
		auto fork = std::find(activeForks.begin(), activeForks.end(), tad);
		activeForks.erase(fork);
		UpdateReplTable(tad);

		
		if ( nb1->isFork() || nb2->isFork() )
		{
			auto fork2 = std::find(activeForks.begin(), activeForks.end(), nb1->isFork() ? nb1 : nb2);
			activeForks.erase(fork2);
			if(nb1->isFork()){
				UpdateReplTable(nb1);//since it's a fork about to be replicated, energy decreases
			}
			else{
				UpdateReplTable(nb2);//since it's a fork about to be replicated, energy decreases
			}
			
		}
	}
	
	ReplicateTADs(tad);
	ReplicateBonds(tad);

	Update();

	
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
		
		if(nb1->isCAR)
		{

			tadConf.back().isCAR=true;
			tadConf.back().isChoesin=true;
			nb1->isChoesin=true;
			nb1->choesin_binding_site = &tadConf.back();
			tadConf.back().choesin_binding_site=nb1;
		}


	}
	
	// Replicate TAD
	tadReplic = *tad;
	tadConf.push_back(tadReplic);
	tad->SisterID= (int) tadConf.size()-1;
	tadConf.back().SisterID = (int) std::distance(tadConf.data(), tad);
	
	if(tad->isCAR)
	{

		tadConf.back().isCAR=true;
		tadConf.back().isChoesin=true;
		tad->isChoesin=true;
		tad->choesin_binding_site = &tadConf.back();
		tadConf.back().choesin_binding_site=tad;
	}

	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		tadConf.push_back(tadReplic);
		nb2->SisterID= (int) tadConf.size()-1;
		tadConf.back().SisterID = (int) std::distance(tadConf.data(), nb2);

		
		if(nb2->isCAR)
		{

			tadConf.back().isCAR=true;
			tadConf.back().isChoesin=true;
			nb2->isChoesin=true;
			nb2->choesin_binding_site = &tadConf.back();
			tadConf.back().choesin_binding_site=nb2;
		}

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
			//Update the choesin status: if in the original chain is CAR so it is in the new chain.
			if(tadConf[tad->SisterID].isCAR)
			{
				tad->isCAR=true;
				tad->isChoesin=true;
				tadConf[tad->SisterID].isChoesin=true;
				tad->choesin_binding_site=&tadConf[tad->SisterID];
				tadConf[tad->SisterID].choesin_binding_site = &tadConf[tadConf[tad->SisterID].SisterID];
				
			}
		}
		
		Ntad = (int) tadConf.size();
	}
	
	Nfork = (int) activeForks.size();
}
void MCReplicPoly::MoveChoesin(MCTad* tad)
{

	
}
double MCReplicPoly::GetEffectiveEnergy() const
{
	if ( Jf > 0.  )
	{
		if (tadTrial->isFork() and neigh==true)
		{
			
			std::vector<int> indexes; //create a indexes vector
			indexes.reserve(activeForks.size());
			for (int i = 0; i < (int)activeForks.size(); ++i)
				indexes.push_back(i); //populate
			std::random_shuffle(indexes.begin(), indexes.end()); // randomize
			
			double Jf1=0.0;
			double Jf2=0.0;
			
			for ( int i = 0; i < (int) indexes.size(); ++i )
			{
				int forkpos = activeForks[indexes[i]]->pos;
				if(forkpos!=tadUpdater->vo)
				{
					for ( int dir = 0; dir < 3; ++dir )
						if(abs(lat->xyzTable[dir][forkpos]-lat->xyzTable[dir][tadUpdater->vo])>6)
							break;
					
					bool neigh1=false;
					bool neigh2=false;
					for ( int v = 0; v < 13; ++v )
					{
						int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
						int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
						
						for ( int v1 = 0; v1 < 13; ++v1)
						{
							int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
							int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
							
							if(forkpos==vi2 )
								neigh2=true;

							if(forkpos==vi1)
								neigh1=true;
						}
						if(neigh1==neigh2 and neigh2==true)
							break;
					}
					if(neigh2==true)
						Jf2+=Jf;
					if(neigh1==true)
						Jf1+=Jf;
					if(Jf2==Jf1 and Jf1==Jf)
						break;
					
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jf2+Jf1;
		}
		if (tadTrial->isFork() and neigh==false)
		{
			
			std::vector<int> indexes; //create a indexes vector
			indexes.reserve(activeForks.size());
			for (int i = 0; i < (int)activeForks.size(); ++i)
				indexes.push_back(i); //populate
			std::random_shuffle(indexes.begin(), indexes.end()); // randomize
			
			double Jf1=0.0;
			double Jf2=0.0;
			
			for ( int i = 0; i < (int) indexes.size(); ++i )
			{
				int forkpos = activeForks[indexes[i]]->pos;
				if(forkpos!=tadUpdater->vo)
				{
					for ( int dir = 0; dir < 3; ++dir )
						if(abs(lat->xyzTable[dir][forkpos]-lat->xyzTable[dir][tadUpdater->vo])>6)
							break;
					
					bool neigh1=false;
					bool neigh2=false;
					for ( int v = 0; v < 13; ++v )
					{
						int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
						int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];

							
						if(forkpos==vn )
							neigh2=true;
							
						if(forkpos==vo)
							neigh1=true;
					}
					if(neigh1==neigh2 and neigh2==true)
						break;
				if(neigh2==true)
					Jf2+=Jf;
				if(neigh1==true)
					Jf1+=Jf;
				}
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jf2+Jf1;
		}
	}
	if ( Jpair > 0.  )
	{
		if (tadTrial->isChoesin == true and neigh==true)
		{


			double Jbott1=0.0;
			double Jbott2=0.0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];


				
				for ( int v1 = 0; v1 < 13; ++v1)
				{
					int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
					int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
				
					if(tadTrial->choesin_binding_site->pos==vi2 ){
						Jbott2=Jpair;
					}
					if(tadTrial->choesin_binding_site->pos==vi1){
						Jbott1=Jpair;
					}
				}
				if(Jbott1==Jbott2 and Jbott2==Jpair)
					break;
			}

			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		if (tadTrial->isChoesin == true and neigh==false)
		{
			
			
			double Jbott1=0.0;
			double Jbott2=0.0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
				int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];

				if(tadTrial->choesin_binding_site->pos==vn )
					Jbott2=Jpair;
				
				if(tadTrial->choesin_binding_site->pos==vo)
					Jbott1=Jpair;

			
				if(Jbott1==Jbott2 and Jbott2==Jpair)
					break;
			}
			return 	MCHeteroPoly::GetEffectiveEnergy() -Jbott2+Jbott1;
		}
		
	}
	
	return 	MCHeteroPoly::GetEffectiveEnergy();
}

void MCReplicPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();

	if ( tadTrial->isFork()) //increase energy at fork site
	{
		std::vector<int> updatedposition1;
		std::vector<int> updatedposition2;

		
		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ? tadUpdater->vo: lat->bitTable[v][tadUpdater->vo];
			int vn =(v == 0) ? tadUpdater->vn: lat->bitTable[v][tadUpdater->vn];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				int vi2 = (v1 == 0) ? vn : lat->bitTable[v1][vn];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					--ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
				if ( std::find(updatedposition2.begin(), updatedposition2.end(), vi2) != updatedposition2.end() )
				{
					++ReplTable[0][vi2];
					updatedposition2.push_back(vi2);
					
				}
			}
		}
	}
	
}

void MCReplicPoly::UpdateReplTable(MCTad* tad)
{
/*
	if(tad->isFork())
	{
		std::vector<int> updatedposition1;
		
		
		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					--ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
			}
		}
	}
	else
	{
		std::vector<int> updatedposition1;

		for ( int v = 0; v < 13; ++v )
		{
			int vo =(v == 0) ?  tad->pos : lat->bitTable[v][tad->pos];
			
			for ( int v1 = 0; v1 < 13; ++v1)
			{
				int vi1 = (v1 == 0) ? vo: lat->bitTable[v1][vo];
				if ( std::find(updatedposition1.begin(), updatedposition1.end(), vi1) != updatedposition1.end() )
				{
					++ReplTable[0][vi1];
					updatedposition1.push_back(vi1);
				}
			}
		}
	}*/
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
