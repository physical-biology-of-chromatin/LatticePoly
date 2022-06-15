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
	
	std::fill(centerMass.begin(), centerMass.end(), 0.);

	if ( RestartFromFile )
		FromVTK(Ninit);
	else if (Rconfinement > 0)
		GenerateRandom(Rconfinement/2);
	else
		GenerateRandom(L/2);

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
	
	if ( !bond.isSet || (tad1->links == 0) )
		++tad1->links;
	if ( !bond.isSet || (tad2->links == 0) )
		++tad2->links;
	
	bond.isSet = true;
}

void MCPoly::GenerateRandom(int lim)
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
	if (Rconfinement > 0)
	{
		double c = (L-0.5)/2;

		for ( int t = 0; t < Ntad; ++t )
		{
			double d2 = SQR(lat->xyzTable[0][tadConf[t].pos]-c)+SQR(lat->xyzTable[1][tadConf[t].pos]-c)+SQR(lat->xyzTable[2][tadConf[t].pos]-c);

			if ( d2 > SQR(Rconfinement) ) throw std::runtime_error("Confinement is screwed up");
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
	tadUpdater->AcceptMove(tadTrial);
	
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
	
	auto types = vtkSmartPointer<vtkFloatArray>::New();
	auto forks = vtkSmartPointer<vtkIntArray>::New();
	auto status = vtkSmartPointer<vtkIntArray>::New();
	auto painters = vtkSmartPointer<vtkFloatArray>::New();
	auto sisterIDs = vtkSmartPointer<vtkIntArray>::New();

	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	forks->SetName("Fork type");
	forks->SetNumberOfComponents(1);
	
	status->SetName("Replication status");
	status->SetNumberOfComponents(1);

	painters->SetName("Painter status");
	painters->SetNumberOfComponents(1);
	
	sisterIDs->SetName("Sister ID");
	sisterIDs->SetNumberOfComponents(1);
	
	std::vector<double3> conf = GetPBCConf();

	for ( int t = 0; t < Ntad; ++t )
	{
		double type = tadConf[t].type;
		int state = tadConf[t].status;
		int id = tadConf[t].sisterID;
		
		double painter = tadConf[t].painter;
		
		int fork = tadConf[t].isFork() ? (tadConf[t].isLeftFork() ? -1 : 1) : 0;
		
		points->InsertNextPoint(conf[t][0], conf[t][1], conf[t][2]);
		
		types->InsertNextValue(type);
		forks->InsertNextValue(fork);
		status->InsertNextValue(state);
		painters->InsertNextValue(painter);
		sisterIDs->InsertNextValue(id);
	}
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
	{
		auto line = vtkSmartPointer<vtkLine>::New();
		
		line->GetPointIds()->SetId(0, bond->id1);
		line->GetPointIds()->SetId(1, bond->id2);
	
		lines->InsertNextCell(line);
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	
	polyData->GetPointData()->AddArray(types);
	polyData->GetPointData()->AddArray(forks);
	polyData->GetPointData()->AddArray(status);
	polyData->GetPointData()->AddArray(painters);
	polyData->GetPointData()->AddArray(sisterIDs);

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
	vtkDataArray* statusData = polyData->GetPointData()->GetArray("Replication status");
	vtkDataArray* sisterData = polyData->GetPointData()->GetArray("Sister ID");
	vtkDataArray* painterData = polyData->GetPointData()->GetArray("Painter status");
	
    Ntad = (int) polyData->GetNumberOfPoints();
	Nbond = (int) polyData->GetNumberOfLines();
	
	tadConf.resize(Ntad);
	tadTopo.resize(Nbond);
			
	for ( int t = 0; t < Ntad; ++t )
	{
		double point[3];
		
		polyData->GetPoint(t, point);
		
		tadConf[t].type = (double) typeData->GetComponent(t, 0);
		tadConf[t].status = (int) statusData->GetComponent(t, 0);
		tadConf[t].sisterID = (int) sisterData->GetComponent(t, 0);
		tadConf[t].painter = (double) painterData->GetComponent(t, 0);

		for ( int i = 0; i < 3; ++i )
		{
			centerMass[i] += point[i] / ((double) Ntad);

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
	
	auto lastBond = std::find_if(tadTopo.begin(), tadTopo.end(), [](const MCBond& b){return b.id2 != b.id1+1;});
	int length = (int) std::distance(tadTopo.begin(), lastBond) + 1;
	
	if ( length != Nchain )
		throw std::runtime_error("MCPoly: Found incompatible main chain dimension " + std::to_string(length));
}

std::vector<double3> MCPoly::GetPBCConf()
{
	std::vector<double3> conf(Ntad);

	for ( int t = 0; t < Ntad; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			conf[t][i] = lat->xyzTable[i][tadConf[t].pos];
	}
	
	for ( auto bond = tadTopo.begin(); bond != tadTopo.end(); ++bond )
		FixPBCPair(conf, bond->id1, bond->id2);
	
	SetPBCCenterMass(conf.begin(), conf.end(), &centerMass);
	
	return conf;
}

void MCPoly::FixPBCPair(std::vector<double3>& conf, int id1, int id2)
{
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

void MCPoly::SetPBCCenterMass(std::vector<double3>::iterator end1, std::vector<double3>::iterator end2, double3* oldCenter)
{
	double3 newCenter = {0., 0., 0.};
	
	bool isSetOldCenter = std::any_of(oldCenter->begin(), oldCenter->end(), [](double x){return x != 0.;});
	bool isSetCenterMass = std::any_of(centerMass.begin(), centerMass.end(), [](double x){return x != 0.;});

	for ( int i = 0; i < 3; ++i )
	{
		for ( auto tadPos = end1; tadPos != end2; ++tadPos )
			newCenter[i] += (*tadPos)[i] / ((double) std::distance(end1, end2));

		if ( isSetOldCenter || isSetCenterMass )
		{
			double deltacenterMass = isSetOldCenter ? newCenter[i] - (*oldCenter)[i] : newCenter[i] - centerMass[i];
		
			// Translate chain center of mass into the same box as the previous conformation
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
