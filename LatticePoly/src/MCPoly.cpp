//
//  MCPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//
#include <set>
#include <iterator>
#include <algorithm>
#include "MCPoly.hpp"
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>





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
	{
		if(Centromere!=0)
			GenerateHedgehog(L/2);
		else
			GenerateRing(L/2);
	}
	

	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		SetBond(*bond);
	

	
	
	std::cout << "Running with initial polymer density " << Ntad / ((double) Ntot) << std::endl;
	std::cout << "Using " << Ntad << " TADs, including main chain of length " << Nchain << std::endl;
	
	
	//mbds=vtkSmartPointer<vtkMultiBlockDataSet>::New();


	// Set the number of timesteps

	

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
		tadConf[t].SisterID = t;
	
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
	
}

void MCPoly::GenerateRing(int lim)
{
	Ntad = Nchain;
	Nbond = Nchain-1;
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
	
	for ( int t = 0; t < Ntad; ++t )
		tadConf[t].SisterID = t;
	
	for ( int b = 0; b < Nbond; ++b )
	{
		tadTopo[b].id1 = b;
		tadTopo[b].id2 = b+1;
	}
	

	
	int vi = 2*CUB(L) + SQR(L) + L/2; // Set to lat->rngEngine() % Ntot for random chromosome placement
	

	
	int ni = 0;
	std::vector<int> turns ={5,2,lat->opp[5],lat->opp[2]};
	for ( int i = 0; i < 4; ++i )
	{
		int turn=turns[i];
		for ( int j = 0; j < lim-1; ++j )
		{
			if(i==0 and j==0)
			{
				tadConf[0].pos = vi;
				lat->bitTable[0][vi] = 1;
				++ni;

			}
			else{
				tadTopo[ni-1].dir = turn;
				tadConf[ni].pos = lat->bitTable[turn][tadConf[ni-1].pos];
				
				lat->bitTable[0][tadConf[ni].pos] = 1;
				
				++ni;
			}
		}
	}
		
	
	--ni;

	while ( ni < Nbond )
	{
		int t = lat->rngEngine() % ni;
		while(t==0 and t==ni)
			t = lat->rngEngine() % ni;

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

	tadConf.back().binding_site=&tadConf.at(0);
	tadConf.at(0).binding_site=&tadConf.back();
	//tadConf.back().neighbors[1]=&tadConf.at(0);
	//tadConf.at(0).neighbors[0]=&tadConf.back();
	//tadConf.at(0).bonds[0]->dir=lat->opp[2];
	//tadConf.back().bonds[1]->dir=lat->opp[2];


	
	
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



void MCPoly::ToVTK(int frame)
{


	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;

	vtkSmartPointer<vtkPolyData> polyData = GetVTKData();

	auto writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer2->SetFileName(path.c_str());
	writer2->SetInputData(polyData);
	
 	writer2->Write();
	
	/*
	vtkSmartPointer<vtkPolyData> polyData = GetVTKData();
	
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->ShallowCopy(polyData);
	
	// Set the time step using the vtkDataSetAttributes object
	vtkSmartPointer<vtkDoubleArray> timeArray = vtkSmartPointer<vtkDoubleArray>::New();
	timeArray->SetName("Time");
	timeArray->InsertNextValue(frame); // Set the time step to the current index
	pd->GetFieldData()->AddArray(timeArray);
	
	// Add the polydata to the multi-block dataset
	mbds->SetBlock(frame, pd);*/
	
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

		int chainId = (chainNum == 1) ? 0 : t / Nchain;
		int chainLength = (chainNum == 1) ? Ntad : Nchain;

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

void MCPoly::Update_rcms_before_separation()
{
	std::cout << "UPDATING BEFORE"<<std::endl;

	auto conf=BuildUnfoldedConf();
	
	
	// create the two separate chromatid in the rcms computetion
	std::vector<double3> SC1;
	std::vector<double3> SC2;
	
	std::cout << conf.size()<<std::endl;

	for ( int i = 0; i < (int)conf.size(); ++i )
	{
		if(tadConf.at(i).status==0 or tadConf.at(i).status==-1)
			SC1.push_back( conf.at(i));
		if(tadConf.at(i).status==0 or tadConf.at(i).status==1)
			SC2.push_back( conf.at(i));
		
	}
	if((int)SC1.size()!=Nchain or (int)SC2.size()!=Nchain)
		throw std::runtime_error("ERROR " + std::to_string((int)SC1.size()));
	
	
	/*std::cout << "old_cm 1"<<std::endl;
	for ( int i = 0; i <3; ++i )
		std::cout << centerMass.at(0)[i]<<"  ";
	std::cout << ""<<std::endl;

	std::cout << "old_cm 2"<<std::endl;
	for ( int i = 0; i <3; ++i )
		std::cout << centerMass.at(1)[i]<<"  ";
	std::cout << ""<<std::endl;*/




	
	double3 Rcm_SC1 ={0,0,0};
	double3 Rcm_SC2 = {0,0,0};
	for ( int i = 0; i <3; ++i )
		for ( int k = 0; k < (int)SC2.size(); ++k )
		{
			Rcm_SC1[i]+=SC1[k][i]/(int)SC2.size();
			Rcm_SC2[i]+=SC2[k][i]/(int)SC2.size();
			
		}
	
	centerMass.at(0)=Rcm_SC1;
	centerMass.at(1)=Rcm_SC2;

	//std::cout << "new_cm 1"<<std::endl;
	//for ( int i = 0; i <3; ++i )
	//	std::cout << Rcm_SC1[i]<<"  ";
	//std::cout << ""<<std::endl;
	
	//std::cout << "new_cm 2"<<std::endl;
	//for ( int i = 0; i <3; ++i )
	//	std::cout << Rcm_SC2[i]<<"  ";
	//std::cout << ""<<std::endl;
	
	double distance=0;
	for ( int i = 0; i <3; ++i )
		distance=distance+SQR(Rcm_SC1[i]-Rcm_SC2[i]);
	
	//std::cout << "distance = "<< sqrt(distance)<< std::endl;
		

}
void MCPoly::FixPBCCenterMass(std::vector<double3>& conf)
{
	int chainNum =  Ntad / Nchain;
	
	int chainLength = (chainNum == 1) ? Ntad : Nchain;
	
	for ( int chainId = 0; chainId < chainNum; ++chainId )
	{
		double3 newCenter = {0., 0., 0.};
		double3* oldCenter = centerMass.begin() + chainId;
		
		

		auto end1 = conf.begin() + chainId*chainLength;
		auto end2 = conf.begin() + (chainId+1)*chainLength;
		
		bool isSetOldCenter = std::any_of(oldCenter->begin(), oldCenter->end(), [](double x){return x != 0.;});
		if ( !isSetOldCenter and chainId>0)
		{
			//NB works only for two chains
			centerMass.at(chainId)=centerMass.at(chainId-1);
			oldCenter = centerMass.begin() + chainId;
		}
		isSetOldCenter = std::any_of(oldCenter->begin(), oldCenter->end(), [](double x){return x != 0.;});


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
bool MCPoly::PrintCohesins()
{

		
		
	std::ofstream outfile_trans(outputDir+"/cohesion_pattern_trans.res", std::ios_base::app | std::ios_base::out);
	std::ofstream outfile_cis1(outputDir+"/cohesion_pattern_cis1.res", std::ios_base::app | std::ios_base::out);
	std::ofstream outfile_cis2(outputDir+"/cohesion_pattern_cis2.res", std::ios_base::app | std::ios_base::out);
	std::vector<int> check;
	std::cout << "PRINTING COHESINS" << std::endl;
	for ( int i = 0; i < Nchain ; ++i )
		if(tadConf.at(i).isCohesin)
		{
			if(tadConf.at(i).binding_site->status!=tadConf.at(i).status)
			{
				std::cout << "Cohesion: SC1 bound at " << i<< "with SC2 at "<<tadConf.at(i).binding_site->SisterID << std::endl;
				outfile_trans << i << std::endl;
				outfile_trans << tadConf.at(i).binding_site->SisterID << std::endl;

			}
			else
			{
				std::cout << "Looping: anchor at " << i<< " binding with anchor at "<<(int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;
				outfile_cis1 << i << std::endl;
				outfile_cis1 << (int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;

			}
			check.push_back((int) std::distance(tadConf.data(), tadConf.at(i).binding_site));
			
		}
	for ( int i = Nchain; i < Ntad ; ++i )
		if(tadConf.at(i).isCohesin)
		{
			
			if(tadConf.at(i).binding_site->status!=tadConf.at(i).status)
			{
				std::cout << "Cohesion: SC2 bound at " << (int) tadConf.at(i).SisterID << "with SC1 at "<< (int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;

			}
			else
			{
				std::cout << "Looping: anchor at " << i<< " binding with anchor at "<<(int) std::distance(tadConf.data(), tadConf.at(i).binding_site) << std::endl;
				outfile_cis2 << tadConf.at(i).SisterID << std::endl;
				outfile_cis2 << tadConf.at((int) std::distance(tadConf.data(), tadConf.at(i).binding_site)).SisterID << std::endl;

			}
			check.push_back((int) std::distance(tadConf.data(), tadConf.at(i).binding_site));
			
			
		}
	for (std::vector<int>::iterator it=check.begin(); it!=check.end(); ++it)
		std::cout << ' ' << *it<< ',';
	std::cout << '\n';

	std::set<int> setOfNumbers(check.begin(), check.end());
	if (setOfNumbers.size() == check.size())
		std::cout<<"Vector has only unique values" <<std::endl;
	else
		std::cout<<"Vector is not unique" <<std::endl;

	
	return 0;
	/*bool bool_v=0;
	for ( int i = 0; i < Nchain ; ++i )
		if(tadConf.at(i).isCohesin)
			if(tadConf.at(i).binding_site!=check[i])
			{
				std::cout << "CHANGED COORDINATE OF MONOMER " << i<< std::endl;
				check[i]=tadConf.at(i).binding_site;
				bool_v=1;
			}
	return bool_v;*/

}
