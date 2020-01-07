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

void MCPoly::Init(std::mt19937_64& rngEngine)
{
	int lim = L/2;
	int turn1[7], turn2[7];
	
	turn1[0] = 12;
	turn1[1] = 12;
	turn1[2] = 1;
	turn1[3] = 1;
	turn1[4] = 11;
	turn1[5] = 11;
	turn1[6] = 2;

	turn2[0] = 12;
	turn2[1] = 1;
	turn2[2] = 1;
	turn2[3] = 11;
	turn2[4] = 11;
	turn2[5] = 2;
	turn2[6] = 2;
	
	for ( int i = 0; i < Nchain; i++ )
	{
		tadType[i] = 0;
		tadConf[i] = -1;
		tadNbId[i] = -1;
	}
	
	// int idx = rngEngine() % Ntot;
	int idx = 2*CUB(L) + SQR(L) + L/2;
	
	tadConf[0] = idx;
	lat->bitTable[0][idx] = 1;
	
	int ni = 1;
	
	for ( int i = 0; i < lim; i++ )
	{
		for ( int j = 0; j < 7; j++ )
		{
			int turn = ((i % 2) == 0) ? turn1[j] : turn2[j];
			
			tadNbId[ni-1] = turn;
			tadConf[ni] = lat->bitTable[turn][tadConf[ni-1]];
			
			lat->bitTable[0][tadConf[ni]] = 1;
			
			ni++;
		}
		
		tadNbId[ni-1] = 10;
		tadConf[ni] = lat->bitTable[10][tadConf[ni-1]];
		
		lat->bitTable[0][tadConf[ni]] = 1;
		
		ni++;
	}
	
	ni--;
	
	while ( ni < Nchain-1 )
	{
		int t  = rngEngine() % ni;
		int iv = rngEngine() % lat->nbNN[0][0][tadNbId[t]];
		
		int nv1 = lat->nbNN[2*iv+1][0][tadNbId[t]];
		int nv2 = lat->nbNN[2*(iv+1)][0][tadNbId[t]];
		
		int en2 = tadConf[t];
		
		int v = (nv1 == 0) ? en2 : lat->bitTable[nv1][en2];
		int b = lat->bitTable[0][v];
					
		if ( b == 0 )
		{
			for ( int i = ni+1; i > t+1; i-- )
			{
				tadConf[i] = tadConf[i-1];
				tadNbId[i] = tadNbId[i-1];
			}
			
			tadConf[t+1] = v;
			tadNbId[t+1] = nv2;
			tadNbId[t] = nv1;

			lat->bitTable[0][v] = 1;
			
			ni++;
		}
	}
	
	centreMass[0] = 0.;
	centreMass[1] = 0.;
	centreMass[2] = 0.;

	for ( int i = 0; i < Nchain; i++ )
	{
		for ( int j = 0; j < 3; j++ )
			centreMass[j] += lat->xyzTable[j][tadConf[i]] / Nchain;
	}
	
	std::cout << "Running with polymer density " << Nchain / ((double) Ntot) << std::endl;
}

void MCPoly::TrialMoveTAD(std::mt19937_64& rngEngine, double* dE)
{
	*dE = 0.;
	
	tad->Init();
	tad->RandomMove(rngEngine, tadConf, tadNbId);
	
	if ( tad->legal )
		*dE = tad->dE;
}

void MCPoly::AcceptMoveTAD()
{
	lat->bitTable[0][tad->en] -= 1;
	lat->bitTable[0][tad->v2] += 1;
	
	if ( tad->n == 0 )
	{
		tadConf[0] = tad->v2;
		tadNbId[0] = lat->opp[tad->iv];
	}
	
	else if ( tad->n == Nchain-1 )
	{
		tadConf[Nchain-1] = tad->v2;
		tadNbId[Nchain-1] = tad->iv;
	}
	
	else
	{
		tadConf[tad->n] = tad->v2;
		tadNbId[tad->n] = tad->nv2;
		tadNbId[tad->n-1] = tad->nv1;
	}
}

void MCPoly::ToVTK(int idx)
{
	char buf[128];
	sprintf(buf, "%05d", idx);
	
	std::string filename = outputDir + "/poly" + buf + ".vtp";
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	vtkSmartPointer<vtkIntArray> types = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkFloatArray> contour = vtkSmartPointer<vtkFloatArray>::New();

	types->SetName("TAD type");
	types->SetNumberOfComponents(1);
	
	contour->SetName("Contour");
	contour->SetNumberOfComponents(1);
	
	double confPBC[3][Nchain];
	double centreMassPBC[3] = {0.,0.,0.};

	for ( int i = 0; i < Nchain; i++ )
	{
		for ( int j = 0; j < 3; j++ )
			confPBC[j][i] = lat->xyzTable[j][tadConf[i]];
		
		if ( i > 0 )
		{
			for ( int j = 0; j < 3; j++ )
			{
				double deltaTad = confPBC[j][i] - confPBC[j][i-1];

				while ( std::abs(deltaTad) > L/2. )
				{
					double pbcShift = std::copysign(L, deltaTad);

					confPBC[j][i] -= pbcShift;
					deltaTad -= pbcShift;
				}
			}
			
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			
			line->GetPointIds()->SetId(0, i-1);
			line->GetPointIds()->SetId(1, i);
			
			lines->InsertNextCell(line);
		}
		
		centreMassPBC[0] += confPBC[0][i] / ((double) Nchain);
		centreMassPBC[1] += confPBC[1][i] / ((double) Nchain);
		centreMassPBC[2] += confPBC[2][i] / ((double) Nchain);
	}
	
	for ( int i = 0; i < 3; i++ )
	{
		double deltaCentreMass = centreMassPBC[i] - centreMass[i];
		
		while ( std::abs(deltaCentreMass) > L/2. )
		{
			double pbcShift = std::copysign(L, deltaCentreMass);
			
			for ( int j = 0; j < Nchain; j++ )
				confPBC[i][j] -= pbcShift;
			
			deltaCentreMass -= pbcShift;
			centreMassPBC[i] -= pbcShift;
		}
		
		centreMass[i] = centreMassPBC[i];
	}

	for ( int i = 0; i < Nchain; i++ )
	{
		int type = tadType[i];
		double curvAbs = i / ((double) Nchain-1);
		
		points->InsertNextPoint(confPBC[0][i], confPBC[1][i], confPBC[2][i]);
		
		types->InsertNextValue(type);
		contour->InsertNextValue(curvAbs);
	}
	
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polydata->SetPoints(points);
	polydata->SetLines(lines);
	
	polydata->GetPointData()->AddArray(types);
	polydata->GetPointData()->AddArray(contour);

	writer->SetFileName(filename.c_str());
	writer->SetInputData(polydata);
	
 	writer->Write();
}
