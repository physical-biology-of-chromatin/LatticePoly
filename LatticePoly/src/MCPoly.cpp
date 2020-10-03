//
//  MCPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <iterator>
#include <algorithm>

#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

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
	tadConf.reserve(2*Ntot);
	tadTopo.reserve(2*Ntot);
	
	if ( RestartFromFile )
		FromVTK(Ninit);
	else
		GenerateRandom(L/2);

	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		CreateBond(*bond);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		MCTad* tad = &tadConf[t];
		
		if ( tad->links == 1 )
		{
			if ( t == tad->bonds[0]->id1 )
				tad->setLeftEnd();
			else if ( t == tad->bonds[0]->id2 )
				tad->setRightEnd();
		}
	}
	
	std::fill(centreMass.begin(), centreMass.end(), 0.);
	
	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			centreMass[i] += lat->xyzTable[i][tadConf[t].pos] / Ntad;
	}
	
	std::cout << "Running with initial polymer density " << Ntad / ((double) Ntot) << std::endl;
	std::cout << "Using " << Ntad << " TADs, including main chain of length " << Nchain << std::endl;
}

void MCPoly::CreateBond(MCBond& bond)
{
	int id1 = bond.id1;
	int id2 = bond.id2;

	MCTad* tad1 = &tadConf[id1];
	MCTad* tad2 = &tadConf[id2];
	
	tad1->neighbors[tad1->links] = tad2;
	tad2->neighbors[tad2->links] = tad1;

	tad1->bonds[tad1->links] = &bond;
	tad2->bonds[tad2->links] = &bond;
	
	++tad1->links;
	++tad2->links;
}

void MCPoly::GenerateRandom(int lim)
{
	Ntad = Nchain;
	Nbond = Nchain-1;
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
	
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
	
	int vi = 2*CUB(L) + SQR(L) + L/2; // Set to rngEngine() % Ntot for random chromosome placement
	
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
		
		int nv1 = lat->nbNN[2*iv+1][0][tadTopo[t].dir];
		int nv2 = lat->nbNN[2*(iv+1)][0][tadTopo[t].dir];
		
		int en2 = tadConf[t].pos;
		int v1 = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
		
		int b = lat->bitTable[0][v1];
					
		if ( b == 0 )
		{
			for ( int i = ni+1; i > t+1; --i )
			{
				tadConf[i].pos = tadConf[i-1].pos;
				tadTopo[i].dir = tadTopo[i-1].dir;
			}
			
			tadConf[t+1].pos = v1;
			
			tadTopo[t].dir = nv1;
			tadTopo[t+1].dir = nv2;

			lat->bitTable[0][v1] = 1;
			
			++ni;
		}
	}
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
	tadUpdater->AcceptMovePos(tadTrial);
	
	--lat->bitTable[0][tadUpdater->vo];
	++lat->bitTable[0][tadUpdater->vn];
}

void MCPoly::ToVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto lines = vtkSmartPointer<vtkCellArray>::New();
	
	auto types = vtkSmartPointer<vtkIntArray>::New();
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto contour = vtkSmartPointer<vtkFloatArray>::New();

	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	forks->SetName("Fork");
	forks->SetNumberOfComponents(1);
	
	contour->SetName("Contour");
	contour->SetNumberOfComponents(1);
	
	std::vector<double3> confPBC(Ntad);
	double3 centreMassPBC = {0.,0.,0.};

	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			confPBC[t][i] = lat->xyzTable[i][tadConf[t].pos];
		
		if ( t > 0 )
		{
			for ( int i = 0; i < 3; ++i )
			{
				double deltaTad = confPBC[t][i] - confPBC[t-1][i];

				while ( std::abs(deltaTad) > L/2. )
				{
					double pbcShift = std::copysign(L, deltaTad);

					confPBC[t][i] -= pbcShift;
					deltaTad -= pbcShift;
				}
			}
		}
		
		for ( int i = 0; i < 3; ++i )
			centreMassPBC[i] += confPBC[t][i] / ((double) Ntad);
	}
	
	for ( int i = 0; i < 3; ++i )
	{
		double deltaCentreMass = centreMassPBC[i] - centreMass[i];
		
		while ( std::abs(deltaCentreMass) > L/2. )
		{
			double pbcShift = std::copysign(L, deltaCentreMass);
			
			for ( int t = 0; t < Ntad; ++t )
				confPBC[t][i] -= pbcShift;
			
			deltaCentreMass -= pbcShift;
			centreMassPBC[i] -= pbcShift;
		}
		
		centreMass[i] = centreMassPBC[i];
	}

	for ( int t = 0; t < Ntad; ++t )
	{
		int type = tadConf[t].type;
		int fork = tadConf[t].isFork();
		
		double curvAbs = t / ((double) Ntad-1);
		
		points->InsertNextPoint(confPBC[t][0], confPBC[t][1], confPBC[t][2]);
		
		types->InsertNextValue(type);
		forks->InsertNextValue(fork);
		contour->InsertNextValue(curvAbs);
	}
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
	{
		if ( (bond->id1 >= 0) && (bond->id2 >= 0) )
		{
			auto line = vtkSmartPointer<vtkLine>::New();
			
			line->GetPointIds()->SetId(0, bond->id1);
			line->GetPointIds()->SetId(1, bond->id2);
		
			lines->InsertNextCell(line);
		}
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	
	polyData->GetPointData()->AddArray(types);
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(contour);

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
	
	vtkCellArray* lineData = polyData->GetLines();
	vtkDataArray* typeData = polyData->GetPointData()->GetArray("TAD type");

	Ntad = (int) polyData->GetNumberOfPoints();
	Nbond = (int) polyData->GetNumberOfLines();
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
			
	for ( int t = 0; t < Ntad; ++t )
	{
		double point[3];
		
		polyData->GetPoint(t, point);
		
		tadConf[t].type = (int) typeData->GetComponent(t, 0);
		
		for ( int i = 0; i < 3; ++i )
		{
			centreMass[i] += point[i] / ((double) Ntad);

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
		auto bond = vtkSmartPointer<vtkIdList>::New();
		
		lineData->GetCellAtId(b, bond);

		int t1 = (int) bond->GetId(0);
		int t2 = (int) bond->GetId(1);

		tadTopo[b].id1 = t1;
		tadTopo[b].id2 = t2;
		
		if ( tadConf[t1].pos == tadConf[t2].pos )
			tadTopo[b].dir = 0;
		
		for ( int v = 0; v < 12; ++v )
		{
			if ( lat->bitTable[v+1][tadConf[t1].pos] == tadConf[t2].pos )
			{
				tadTopo[b].dir = v+1;
				break;
			}
		}
	}
	
	auto lastBond = std::find_if(tadTopo.begin(), tadTopo.end(), [](const MCBond& b){return b.id2 != b.id1+1;});
	int length = (int) std::distance(tadTopo.begin(), lastBond) + 1;
	
	if ( length != Nchain )
		throw std::runtime_error("MCPoly: Found incompatible main chain dimension " + std::to_string(length));
}
