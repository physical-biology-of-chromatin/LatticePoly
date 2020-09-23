//
//  MCPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "MCPoly.hpp"


MCPoly::MCPoly(MCLattice* _lat): lat(_lat)
{
	tad = new MCTad(lat);
}

MCPoly::~MCPoly()
{
	delete tad;
}

void MCPoly::Init(int Ninit)
{
	for ( int t = 0; t < Nchain; ++t )
	{
		tadType[t] = 0;
		tadConf[t] = -1;
		
		if ( t > 0 )
			tadBond[t-1] = -1;
	}

	for ( int i = 0; i < 3; ++i )
		centreMass[i] = 0.;
	
	if ( RestartFromFile )
		FromVTK(Ninit);
	else
		GenerateRandom(L/2);
	
	std::cout << "Running with polymer density " << Nchain / ((double) Ntot) << std::endl;
}

void MCPoly::GenerateRandom(int lim)
{
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
	
	tadConf[0] = vi;
	lat->bitTable[0][vi] = 1;
	
	int ni = 1;
	
	for ( int i = 0; i < lim; ++i )
	{
		for ( int j = 0; j < 7; ++j )
		{
			int turn = ((i % 2) == 0) ? turn1[j] : turn2[j];
			
			tadBond[ni-1] = turn;
			tadConf[ni] = lat->bitTable[turn][tadConf[ni-1]];
			
			lat->bitTable[0][tadConf[ni]] = 1;
			
			++ni;
		}
		
		tadBond[ni-1] = 10;
		tadConf[ni] = lat->bitTable[10][tadConf[ni-1]];
		
		lat->bitTable[0][tadConf[ni]] = 1;
		
		++ni;
	}
	
	--ni;
	
	while ( ni < Nchain-1 )
	{
		int t = lat->rngEngine() % ni;
		int iv = lat->rngEngine() % lat->nbNN[0][0][tadBond[t]];
		
		int nv1 = lat->nbNN[2*iv+1][0][tadBond[t]];
		int nv2 = lat->nbNN[2*(iv+1)][0][tadBond[t]];
		
		int en2 = tadConf[t];
		int v1 = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
		
		int b = lat->bitTable[0][v1];
					
		if ( b == 0 )
		{
			for ( int i = ni+1; i > t+1; --i )
			{
				tadConf[i] = tadConf[i-1];
				tadBond[i] = tadBond[i-1];
			}
			
			tadConf[t+1] = v1;
			
			tadBond[t] = nv1;
			tadBond[t+1] = nv2;

			lat->bitTable[0][v1] = 1;
			
			++ni;
		}
	}

	for ( int t = 0; t < Nchain; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			centreMass[i] += lat->xyzTable[i][tadConf[t]] / Nchain;
	}
}

void MCPoly::TrialMove(double* dE)
{
	tad->TrialMovePos(tadConf, tadBond);
	
	*dE = tad->legal ? tad->dE : 0.;
}

void MCPoly::AcceptMove()
{
	tad->AcceptMovePos(tadConf, tadBond);
	
	--lat->bitTable[0][tad->vo];
	++lat->bitTable[0][tad->vn];
}

void MCPoly::ToVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "poly%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto lines = vtkSmartPointer<vtkCellArray>::New();
	
	auto types = vtkSmartPointer<vtkIntArray>::New();
	auto contour = vtkSmartPointer<vtkFloatArray>::New();

	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	contour->SetName("Contour");
	contour->SetNumberOfComponents(1);
	
	double confPBC[3][Nchain];
	
	double centreMassPBC[3] = {0.,0.,0.};

	for ( int t = 0; t < Nchain; ++t )
	{
		for ( int i = 0; i < 3; ++i )
			confPBC[i][t] = lat->xyzTable[i][tadConf[t]];
		
		if ( t > 0 )
		{
			for ( int i = 0; i < 3; ++i )
			{
				double deltaTad = confPBC[i][t] - confPBC[i][t-1];

				while ( std::abs(deltaTad) > L/2. )
				{
					double pbcShift = std::copysign(L, deltaTad);

					confPBC[i][t] -= pbcShift;
					deltaTad -= pbcShift;
				}
			}
			
			auto line = vtkSmartPointer<vtkLine>::New();
			
			line->GetPointIds()->SetId(0, t-1);
			line->GetPointIds()->SetId(1, t);
			
			lines->InsertNextCell(line);
		}
		
		for ( int i = 0; i < 3; ++i )
			centreMassPBC[i] += confPBC[i][t] / ((double) Nchain);
	}
	
	for ( int i = 0; i < 3; ++i )
	{
		double deltaCentreMass = centreMassPBC[i] - centreMass[i];
		
		while ( std::abs(deltaCentreMass) > L/2. )
		{
			double pbcShift = std::copysign(L, deltaCentreMass);
			
			for ( int t = 0; t < Nchain; ++t )
				confPBC[i][t] -= pbcShift;
			
			deltaCentreMass -= pbcShift;
			centreMassPBC[i] -= pbcShift;
		}
		
		centreMass[i] = centreMassPBC[i];
	}

	for ( int t = 0; t < Nchain; ++t )
	{
		int type = tadType[t];
		
		double curvAbs = t / ((double) Nchain-1);
		
		points->InsertNextPoint(confPBC[0][t], confPBC[1][t], confPBC[2][t]);
		
		types->InsertNextValue(type);
		contour->InsertNextValue(curvAbs);
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	
	polyData->GetPointData()->AddArray(types);
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
	
	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(path.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	vtkDataArray* typeData = polyData->GetPointData()->GetArray("TAD type");

	vtkIdType Npoints = polyData->GetNumberOfPoints();
		
	if ( Npoints != Nchain )
		throw std::runtime_error("MCPoly: Found polymer configuration file with incompatible dimension " + std::to_string(Npoints));
	else
		std::cout << "Starting from polymer configuration file " << path << std::endl;
		
	for ( int t = 0; t < Nchain; ++t )
	{
		double point[3];
		
		polyData->GetPoint(t, point);
		
		tadType[t] = (int) typeData->GetComponent(t, 0);
		
		for ( int i = 0; i < 3; ++i )
		{
			centreMass[i] += point[i] / ((double) Nchain);

			while ( point[i] >= L ) point[i] -= L;
			while ( point[i] < 0 )  point[i] += L;
		}

		int ixp = (int) 1*point[0];
		int iyp = (int) 2*point[1];
		int izp = (int) 4*point[2];
		
		tadConf[t] = ixp + iyp*L + izp*L2;
		
		++lat->bitTable[0][tadConf[t]];
		
		if ( t > 0 )
		{
			if ( tadConf[t] == tadConf[t-1] )
				tadBond[t-1] = 0;

			else
			{
				for ( int v = 0; v < 12; ++v )
				{
					if ( lat->bitTable[v+1][tadConf[t-1]] == tadConf[t] )
					{
						tadBond[t-1] = v+1;
						break;
					}
				}
			}
		}
	}
}
