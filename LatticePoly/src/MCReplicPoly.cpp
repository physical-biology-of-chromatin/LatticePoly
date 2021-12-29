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
	
	activeForks.reserve(Nchain);

	// Locate existing forks
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->isFork() )
			activeForks.push_back(&(*tad));
	}
		
	// Set origin locations	
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( (tad->status == 0) && !tad->isLeftEnd() && !tad->isRightEnd() )
			inactiveOrigins.push_back(&(*tad));
	}
	
	Nfork = (int) activeForks.size();
	Norigin = (int) inactiveOrigins.size();
}

void MCReplicPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);

	if ( Norigin > 0 )
	{
		// Nucleate replication bubble at random inactive origin with rate originRate
		double rndOrigin = lat->rngDistrib(lat->rngEngine);
	
		if ( rndOrigin < originRate / (double) Ntad )
		{
			int o = lat->rngEngine() % Norigin;
			MCTad* tad = inactiveOrigins[o];
		
			Replicate(tad);
		}
	}
	
	if ( Nfork > 0 )
	{
		// Move forks (i.e. replicate them) with rate replicRate
		double rndReplic = lat->rngDistrib(lat->rngEngine);

		if ( rndReplic < replicRate / (double) Ntad )
		{
			std::vector<MCTad*> _activeForks = activeForks;
			
			for ( auto tadptr = _activeForks.begin(); tadptr != _activeForks.end(); ++tadptr )
			{
				// Check whether fork has already coalesced
				if ( (*tadptr)->isFork() )
					Replicate(*tadptr);
			}
		}
	}
}

void MCReplicPoly::Replicate(MCTad* tad)
{
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
				activeForks.push_back(nb1);
			
			if ( nb2->isRightEnd() )
			{
				if ( rnd < 0.5 )
					return;
			}
			
			else
				activeForks.push_back(nb2);
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
				activeForks.push_back(nb1);
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
				activeForks.push_back(nb2);
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
	}
	
	// Replicate TAD
	tadReplic = *tad;
		
	tadConf.push_back(tadReplic);
	
	// Same for right end/fork
	if ( nb2->isRightEnd() || nb2->isLeftFork() )
	{
		tadReplic = *nb2;
		
		tadConf.push_back(tadReplic);
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
		
		SetBond(*bond1);
		UnsetFork(tad);
		
		// Merge forks if necessary
		if ( nb2->isLeftFork() )
		{
			MCBond* bond3 = nb2->bonds[2];
			bond3->id1 = bondReplic2.id2;
			
			SetBond(*bond3);
			UnsetFork(nb2);
		}
	}
	
	else
		tadTopo.push_back(bondReplic1);
			
	// Same for left forks
	if ( tad->isLeftFork() )
	{
		bond2->id1 = bondReplic2.id1;
		
		SetBond(*bond2);
		UnsetFork(tad);
		
		if ( nb1->isRightFork() )
		{
			MCBond* bond3 = nb1->bonds[2];
			bond3->id2 = bondReplic1.id1;
			
			SetBond(*bond3);
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
			SetBond(*bond);
		
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
	
	// Update origins
	inactiveOrigins.erase(std::remove_if(inactiveOrigins.begin(), inactiveOrigins.end(), [](const MCTad* tad){return tad->status != 0;}),
						  inactiveOrigins.end());
	
	// Update fork/origin counters
	Nfork = (int) activeForks.size();
	Norigin = (int) inactiveOrigins.size();
}

vtkSmartPointer<vtkPolyData> MCReplicPoly::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCHeteroPoly::GetVTKData();
	
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto status = vtkSmartPointer<vtkIntArray>::New();
	auto sisterID = vtkSmartPointer<vtkIntArray>::New();
	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);
	
	sisterID->SetName("Sister ID");
	sisterID->SetNumberOfComponents(1);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		int state = tadConf[t].status;
		int id = tadConf[t].sisterID;
		
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
				
		forks->InsertNextValue(fork);
		status->InsertNextValue(state);
		sisterID->InsertNextValue(id);
	}
		
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(status);
	polyData->GetPointData()->AddArray(sisterID);
	
	return polyData;
}

void MCReplicPoly::SetVTKData(vtkSmartPointer<vtkPolyData> polyData)
{
	MCHeteroPoly::SetVTKData(polyData);

	vtkDataArray* statusData = polyData->GetPointData()->GetArray("Replication status");
	vtkDataArray* sisterData = polyData->GetPointData()->GetArray("Sister ID");
	
	for ( int t = 0; t < Ntad; ++t )
	{
		tadConf[t].status = (int) statusData->GetComponent(t, 0);
		tadConf[t].sisterID = (int) sisterData->GetComponent(t, 0);
	}
}
