//
//  MCPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>

#include "MCPoly.hpp"


MCPoly::MCPoly(MCLattice* _lat): lat(_lat)
{
	tadUpdater = new MCTadUpdater(lat);
}

MCPoly::~MCPoly()
{
	delete tadUpdater;
}

void MCPoly::Init(int Ninit)
{
	tadConf.reserve(2*Nchain);
	tadTopo.reserve(2*Nchain);
	
	std::fill(centerMass.begin(), centerMass.end(), (double3) {0., 0., 0.});

	if ( RestartFromFile )
		FromVTK(Ninit);
	else
		GenerateHedgehog(L/2);

	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		SetBond(*bond);
	
	std::cout << "Running with initial polymer density " << Ntad / ((double) Ntot) << std::endl;
	std::cout << "Using " << Ntad << " TADs, including main chain of length " << Nchain << std::endl;
}

void MCPoly::SetBond(MCBond& bond)
{
	MCTad* tad1 = &tadConf[bond.id1];
	MCTad* tad2 = &tadConf[bond.id2];
	
	int id1 = ( !bond.isSet && (tad1->links == 2) ) ? 2 : 1;
	int id2 = ( !bond.isSet && (tad2->links == 2) ) ? 2 : 0;

	tad1->neighbors[id1] = tad2;
	tad2->neighbors[id2] = tad1;

	tad1->bonds[id1] = &bond;
	tad2->bonds[id2] = &bond;
	
	if ( !bond.isSet || (tad1->links == 0) ) ++tad1->links;
	if ( !bond.isSet || (tad2->links == 0) ) ++tad2->links;
	
	bond.isSet = true;
}

void MCPoly::GenerateHedgehog(int lim)
{
	Ntad = Nchain;
	Nbond = Nchain-1;
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
	
	for ( int t = 0; t < Ntad; ++t )
		tadConf[t].sisterID = t;

	for ( int b = 0; b < Nbond; ++b )
	{
		tadTopo[b].id1 = b;
		tadTopo[b].id2 = b+1;
	}
	
	int turn1[7];
	
	turn1[0] = 12;
	turn1[1] = 12;
	turn1[2] = 1;
	turn1[3] = 1;
	turn1[4] = 11;
	turn1[5] = 11;
	turn1[6] = 2;

	int turn2[7];
	
	turn2[0] = 12;
	turn2[1] = 1;
	turn2[2] = 1;
	turn2[3] = 11;
	turn2[4] = 11;
	turn2[5] = 2;
	turn2[6] = 2;
	
	int vi = 2*CUB(L) + SQR(L) + L/2; // Set to lat->rngEngine() % Ntot for random chromosome placement
	
	tadConf[0].pos = vi;
	lat->bitTable[0][vi] = 1;
	
	int ni = 1;
	
	for ( int i = 0; i < lim; ++i )
	{
		for ( int j = 0; j < 7; ++j )
		{
			int turn = ((i % 2) == 0) ? turn1[j] : turn2[j];
			
			tadTopo[ni-1].dir = turn;
			tadConf[ni].pos = lat->bitTable[turn][tadConf[ni-1].pos];
			
			lat->bitTable[0][tadConf[ni].pos] = 1;
			
			++ni;
		}
		
		tadTopo[ni-1].dir = 10;
		tadConf[ni].pos = lat->bitTable[10][tadConf[ni-1].pos];
		
		lat->bitTable[0][tadConf[ni].pos] = 1;
		
		++ni;
	}
	
	--ni;
	
	while ( ni < Nbond )
	{
		int t = lat->rngEngine() % ni;
		int iv = lat->rngEngine() % lat->nbNN[0][0][tadTopo[t].dir];
		
		int nd1 = lat->nbNN[2*iv+1][0][tadTopo[t].dir];
		int nd2 = lat->nbNN[2*(iv+1)][0][tadTopo[t].dir];
		
		int en2 = tadConf[t].pos;
		int v1 = (nd1 == 0) ? en2 : lat->bitTable[nd1][en2];
		
		int b = lat->bitTable[0][v1];
					
		if ( b == 0 )
		{
			for ( int i = ni+1; i > t+1; --i )
			{
				tadConf[i].pos = tadConf[i-1].pos;
				tadTopo[i].dir = tadTopo[i-1].dir;
			}
			
			tadConf[t+1].pos = v1;
			
			tadTopo[t].dir = nd1;
			tadTopo[t+1].dir = nd2;

			lat->bitTable[0][v1] = 1;
			
			++ni;
		}
	}
	
	id_cut1 = 1015;     
    tadTopo.erase(tadTopo.begin() + id_cut1);
    --Nbond;

	// Random Walk Configuration
    // Ntad = Nchain;
	// Nbond = Nchain-1;
	
	// tadConf.resize(Ntad);
	// tadTopo.resize(Nbond);
	
	// for ( int t = 0; t < Ntad; ++t )
	// 	tadConf[t].sisterID = t;

	// for ( int b = 0; b < Nbond; ++b )
	// {
	// 	tadTopo[b].id1 = b;
	// 	tadTopo[b].id2 = b+1;
	// }

    // int vi = 2*CUB(L) + SQR(L) + L/2; // Set to lat->rngEngine() % Ntot for random chromosome placement
	// tadConf[0].pos = vi;
	// lat->bitTable[0][vi] = 1;
	
    // int stuck = 0;
    // int j = 1;
	// for ( int i = 1; i < Nchain; ++i )
	// {
	// 	while (j == i)
    //     {
    //         int dir = lat->rngEngine() % 13;
    //         while (dir == 0)  // To not have double occupancy in the initial config
    //             dir = lat->rngEngine() % 13;
                
	// 	    int previouspos = tadConf[i-1].pos;
	// 	    int next_pos= dir==0? previouspos : lat->bitTable[dir][previouspos];
    //         //int next_pos = lat->bitTable[dir][previouspos];
    //         if (lat->bitTable[0][next_pos] == 0) // Excluded volume
    //             {
	// 	    		lat->bitTable[0][next_pos] = 1;
	// 	    		tadTopo[i-1].dir = dir;
	// 	    		tadConf[i].pos=next_pos;
	// 				stuck = 0;
    //         		j++;
    //             }
    //         else 
    //             stuck++;
    //         if (stuck == 12) // trashing the configuration and resetting 
    //             {
	// 				for ( int k = 1; k < i; ++k )
	// 					{
	// 						tadTopo[k].dir = 0;
	// 						tadConf[k].pos = 0;      
	// 					}
					
	// 				for ( int k = 0; k < Ntot; ++k )
	// 					{
	// 						lat->bitTable[0][k] = 0;   
	// 					}      
	// 				std::cout << "Trashing and resetting at " << i << std::endl;
	// 				i = 0;
	// 				j = 1; 
	// 				break; 
    //             }       	
	//     }
	// }	
}

void MCPoly::TrialMove(double* dE)
{
	int t = lat->rngEngine() % Ntad;
	tadTrial = &tadConf[t];
	
	tadUpdater->TrialMove(tadTrial, dE);
	*dE = tadUpdater->legal ? *dE : 0.;
}

void MCPoly::AcceptMove()
{
	tadUpdater->AcceptMove(tadTrial);
		
	--lat->bitTable[0][tadUpdater->vo];
	++lat->bitTable[0][tadUpdater->vn];

}

void MCPoly::TrialMoveTopo(double* dT)
{
	int ti = lat->rngEngine() % Ntad;
	tadi = &tadConf[ti];
	
	tadUpdater->TrialMoveTopo(tadi, dT); 
	
	if ( tadUpdater -> legalTopo1 )
	{
		for (int t = 0; t < Ntad; ++t)
		{
			if( tadConf.at(t).pos == tadUpdater -> vin  )
			{
				tadx = &tadConf[t];
				tadUpdater->TrialSwapTopo(tadi,tadx, dT);
				break;
			}
		}
	}
	
}

void MCPoly::AcceptMoveTopo()
{
	tadUpdater->AcceptMoveTopo(tadi,tadx);
	
	// //Connectivity Check
	// for (int i = 0; i < Ntad-1; ++i)
	// {
	// 	bool check=false;
	// 	for ( int dir = 1; dir < 13 ; ++dir )
	// 	{
	// 		if(tadConf.at(i).pos==lat->bitTable[dir][tadConf.at(i+1).pos])
	// 			check=true;
	// 		if(tadConf.at(i).pos==tadConf.at(i+1).pos)
	// 		{
	// 			check=true;
	// 		}	
	// 	}		
	// 	if(check==false)
	// 		throw std::runtime_error("OH NO TOPO STUPIDO");
	// }
}

void MCPoly::ToVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	vtkSmartPointer<vtkPolyData> polyData = GetVTKData();
	
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName(path.c_str());
	writer->SetInputData(polyData);
	
 	writer->Write();
}

void MCPoly::FromVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	std::cout << "Starting from polymer configuration file " << path << std::endl;

	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(path.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	
	SetVTKData(polyData);
	
	auto lastBond = std::find_if(tadTopo.begin(), tadTopo.end(), [](const MCBond& b){return b.id2 != b.id1+1;});
	int length = (int) std::distance(tadTopo.begin(), lastBond) + 1;
	
	if ( length != Nchain )
		throw std::runtime_error("MCPoly: Found incompatible main chain dimension " + std::to_string(length));
}

vtkSmartPointer<vtkPolyData> MCPoly::GetVTKData()
{
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto lines = vtkSmartPointer<vtkCellArray>::New();
	
	std::vector<double3> conf = BuildUnfoldedConf();
	
	for ( int t = 0; t < Ntad; ++t )
		points->InsertNextPoint(conf[t][0], conf[t][1], conf[t][2]);
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
	{
		auto line = vtkSmartPointer<vtkLine>::New();
		
		line->GetPointIds()->SetId(0, bond->id1);
		line->GetPointIds()->SetId(1, bond->id2);
		
		lines->InsertNextCell(line);
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	
	return polyData;
}

void MCPoly::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	Ntad = (int) polyData->GetNumberOfPoints();
	Nbond = (int) polyData->GetNumberOfLines();
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
	
	vtkCellArray* lineData = polyData->GetLines();

	for ( int t = 0; t < Ntad; ++t )
	{
		double point[3];
		
		polyData->GetPoint(t, point);
		
		int chainNum = Ntad / Nchain;
		
		if ( id_cut1 > 0 )
			chainNum = 2;	

		int chainId = (chainNum == 1) ? 0 : t / Nchain;
		
		if ( id_cut1 > 0 && chainNum == 2 )
			chainId = ( t < Nchain/2) ? 0 : 1;

		int chainLength = (chainNum == 1) ? Ntad : Nchain/chainNum;

		for ( int i = 0; i < 3; ++i )
		{
			centerMass[chainId][i] += point[i] / ((double) chainLength);

			while ( point[i] >= L ) point[i] -= L;
			while ( point[i] < 0 )  point[i] += L;
		}

		int ixp = (int) 1*point[0];
		int iyp = (int) 2*point[1];
		int izp = (int) 4*point[2];
			
		tadConf[t].pos = ixp + iyp*L + izp*L2;
			
		++lat->bitTable[0][tadConf[t].pos];
	}
	
	for ( int b = 0; b < Nbond; ++b )
	{
		auto cell = vtkSmartPointer<vtkIdList>::New();
		
		lineData->GetCellAtId(b, cell);

		int t1 = (int) cell->GetId(0);
		int t2 = (int) cell->GetId(1);

		tadTopo[b].id1 = t1;
		tadTopo[b].id2 = t2;
		
		if ( tadConf[t1].pos == tadConf[t2].pos )
			tadTopo[b].dir = 0;
		
		else
		{
			for ( int v = 0; v < 12; ++v )
			{
				if ( lat->bitTable[v+1][tadConf[t1].pos] == tadConf[t2].pos )
				{
					tadTopo[b].dir = v+1;
					break;
				}
			}
		}
	}
}

std::vector<double3> MCPoly::BuildUnfoldedConf()
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
					builtTads.push_back(tad2);
					FixPBCPair(conf, tad1, tad2);

					// Traverse side branches
					if ( tad2->isFork() )
					{
						MCTad *tad3, *tad4;
						tad3 = tad2->neighbors[2];
						
						builtTads.push_back(tad3);
						FixPBCPair(conf, tad2, tad3);

						while ( (tad4 = (tad2->isLeftFork() ? tad3->neighbors[1] : tad3->neighbors[0])) )
						{
							builtTads.push_back(tad4);
							FixPBCPair(conf, tad3, tad4);

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
	
	// Translate chain center(s) of mass back into the appropriate box
	FixPBCCenterMass(conf);
	
	return conf;
}

void MCPoly::FixPBCPair(std::vector<double3>& conf, MCTad* tad1, MCTad* tad2)
{
	int id1 = (int) std::distance(tadConf.data(), tad1);
	int id2 = (int) std::distance(tadConf.data(), tad2);
	
	double3* pos1 = &conf[id1];
	double3* pos2 = &conf[id2];

	for ( int i = 0; i < 3; ++i )
	{
		double deltaTad = (*pos2)[i] - (*pos1)[i];
		
		while ( std::abs(deltaTad) > L/2. )
		{
			double pbcShift = std::copysign(L, deltaTad);

			(*pos2)[i] -= pbcShift;
			deltaTad -= pbcShift;
		}
	}
}

void MCPoly::FixPBCCenterMass(std::vector<double3>& conf)
{
	int chainNum = Ntad / Nchain;
	if ( id_cut1 > 0 )
			chainNum = 2;		
	
	int chainLength = (chainNum == 1) ? Ntad : Nchain/chainNum;
	
	for ( int chainId = 0; chainId < chainNum; ++chainId )
	{
		double3 newCenter = {0., 0., 0.};
		double3* oldCenter = centerMass.begin() + chainId;
		
		auto end1 = conf.begin() + chainId*chainLength;
		auto end2 = conf.begin() + (chainId+1)*chainLength;
		
		bool isSetOldCenter = std::any_of(oldCenter->begin(), oldCenter->end(), [](double x){return x != 0.;});

		for ( int i = 0; i < 3; ++i )
		{
			for ( auto tadPos = end1; tadPos != end2; ++tadPos )
				newCenter[i] += (*tadPos)[i] / ((double) chainLength);

			// Translate chain center of mass into the same box as the previous conformation, if it exists
			if ( isSetOldCenter )
			{
				double deltacenterMass = newCenter[i] - (*oldCenter)[i];
			
				while ( std::abs(deltacenterMass) > L/2. )
				{
					double pbcShift = std::copysign(L, deltacenterMass);
				
					for ( auto tadPos = end1; tadPos != end2; ++tadPos )
						(*tadPos)[i] -= pbcShift;
				
					deltacenterMass -= pbcShift;
					newCenter[i] -= pbcShift;
				}
			}
		}
		
		*oldCenter = newCenter;
	}
}
