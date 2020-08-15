//
//  MCLiqLattice.cpp
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

#include "MCLiqLattice.hpp"


void MCLiqLattice::Init(int Ninit)
{
	nLiq = 0;
	
	MCLattice::Init(Ninit);

	for ( int i = 0; i < Ntot; i++ )
	{
		spinTable[i] = 0;
		spinIdTable[i] = -1;
	}
	
	if ( RestartFromFile )
		FromVTK(Ninit);

	else
	{
		if ( InitDrop )
			GenerateDroplets();
		else
			GenerateRandom();
		
		for ( int i = 0; i < nLiq; i++ )
		{
			disp initDisp;
			
			initDisp.dx = 0.;
			initDisp.dy = 0.;
			initDisp.dz = 0.;

			spinDisp.push_back(initDisp);
		}
		
		for ( int i = 0; i < Ntot; i++ )
		{
			if ( spinTable[i] > 0 )
			{
				spinIdTable[i] = (int) spinConf.size();
				spinConf.push_back(i);
			}
		}
	}
	
	std::cout << "Set up lattice with fixed liquid density " << nLiq / ((double) Ntot) << std::endl;
}

void MCLiqLattice::GenerateDroplets()
{
	int centers[3][Ndrop];
	int r = std::floor(R)+1; // Set to 1 to allow initial droplets to cross PBCs
	
	for ( int i = 0; i < Ndrop; i++ )
	{
		centers[0][i] = (rngEngine() % (L-2*r+1)) + r;
		centers[1][i] = (rngEngine() % (L-2*r+1)) + r;
		centers[2][i] = (rngEngine() % (L-2*r+1)) + r;
	}
	
	for ( int i = 0; i < Ntot; i++ )
	{
		for ( int j = 0; j < Ndrop; j++ )
		{
			double dx = xyzTable[0][i] - centers[0][j];
			double dy = xyzTable[1][i] - centers[1][j];
			double dz = xyzTable[2][i] - centers[2][j];

			if      ( SQR(dx-L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

			else if ( SQR(dx)   + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;

			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[i] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[i] = 1;
			
			if ( spinTable[i] > 0 )
			{
				nLiq++;
				break;
			}
		}
	}
}

void MCLiqLattice::GenerateRandom()
{
	int nLiqTot = std::floor(Ntot*Ldens);
	
	while ( nLiq < nLiqTot )
	{
		int idx = rngEngine() % Ntot;
		
		if ( spinTable[idx] == 0 )
		{
			spinTable[idx]++;
			nLiq++;
		}
	}
}

void MCLiqLattice::TrialMove(double* dE)
{
	int v1 = rngEngine() % nLiq;
	int v2 = rngEngine() % 12;

	idx1 = spinConf[v1];
	idx2 = bitTable[v2+1][idx1];
	
	*dE = GetSpinEnergy();	
}

void MCLiqLattice::AcceptMove()
{
	int spin1 = spinTable[idx1];
	int spin2 = spinTable[idx2];
	
	if ( (spin1 == 1) || (spin2 == 1) )
	{
		DisplaceSpins();

		spinTable[idx1] = spin2;
		spinTable[idx2] = spin1;
						
		if ( (spin1 == 1) && (spin2 == 0) )
		{
			int id1 = spinIdTable[idx1];
		
			spinIdTable[idx1] = -1;
			spinIdTable[idx2] = id1;
			
			spinConf[id1] = idx2;
		}
		
		else if ( (spin1 == 0) && (spin2 == 1) )
		{
			int id2 = spinIdTable[idx2];
			
			spinIdTable[idx1] = id2;
			spinIdTable[idx2] = -1;
				
			spinConf[id2] = idx1;
		}
		
		else
		{
			int id1 = spinIdTable[idx1];
			int id2 = spinIdTable[idx2];
			
			spinIdTable[idx1] = id2;
			spinIdTable[idx2] = id1;
			
			spinConf[id1] = idx2;
			spinConf[id2] = idx1;
		}
	}
}

void MCLiqLattice::DisplaceSpins()
{	
	double dx = xyzTable[0][idx2] - xyzTable[0][idx1];
	double dy = xyzTable[1][idx2] - xyzTable[1][idx1];
	double dz = xyzTable[2][idx2] - xyzTable[2][idx1];
	
	if ( std::abs(dx) > L/2. ) dx -= std::copysign(L, dx);
	if ( std::abs(dy) > L/2. ) dy -= std::copysign(L, dy);
	if ( std::abs(dz) > L/2. ) dz -= std::copysign(L, dz);
	
	if ( spinTable[idx1] == 1 )
	{
		int id1 = spinIdTable[idx1];

		spinDisp[id1].dx += dx;
		spinDisp[id1].dy += dy;
		spinDisp[id1].dz += dz;
	}
	
	if ( spinTable[idx2] == 1 )
	{
		int id2 = spinIdTable[idx2];

		spinDisp[id2].dx -= dx;
		spinDisp[id2].dy -= dy;
		spinDisp[id2].dz -= dz;
	}
}

double MCLiqLattice::GetSpinEnergy() const
{
	double dE = 0.;
	
	if ( spinTable[idx1] != spinTable[idx2] )
	{
		for ( int i = 0; i < 12; i++ )
		{
			if ( bitTable[i+1][idx1] != idx2 )
				dE += spinTable[bitTable[i+1][idx1]];
			if ( bitTable[i+1][idx2] != idx1 )
				dE -= spinTable[bitTable[i+1][idx2]];
		}
	
		dE *= Jll * (spinTable[idx1]-spinTable[idx2]);
	}
	
	return dE;
}

double MCLiqLattice::GetCouplingEnergy(const int tadHetTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( spinTable[idx1] != spinTable[idx2] )
	{
		for ( int i = 0; i < 13; i++ )
		{
			int v1 = (i == 0) ? idx1 : bitTable[i][idx1];
			int v2 = (i == 0) ? idx2 : bitTable[i][idx2];
			
			E1 -= tadHetTable[v1];
			E2 -= tadHetTable[v2];
		}
	}
	
	return Jlp * (E2-E1);
}

void MCLiqLattice::ToVTK(int frame)
{
	char buf[128];
	sprintf(buf, "%05d", frame);
	
	std::string filename = outputDir + "/liq" + buf + ".vtp";
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	
	vtkSmartPointer<vtkFloatArray> liqDensity = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> liqDisplacement = vtkSmartPointer<vtkFloatArray>::New();
	
	liqDensity->SetName("Density");
	liqDensity->SetNumberOfComponents(1);
	
	liqDisplacement->SetName("Displacement");
	liqDisplacement->SetNumberOfComponents(3);
		
	for ( int i = 0; i < nLiq; i++ )
	{
		int idx = spinConf[i];
		double aveDensity = 0.;

		for ( int j = 0; j < 12; j++ )
			aveDensity += spinTable[bitTable[j+1][idx]] / 12.;
		
		double x = xyzTable[0][idx];
		double y = xyzTable[1][idx];
		double z = xyzTable[2][idx];
		
		double dx = spinDisp[i].dx;
		double dy = spinDisp[i].dy;;
		double dz = spinDisp[i].dz;;
				
		points->InsertNextPoint(x, y, z);
	
		liqDisplacement->InsertNextTuple3(dx, dy, dz);
		liqDensity->InsertNextValue(aveDensity);
	}
	
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	
	polyData->GetPointData()->AddArray(liqDensity);
	polyData->GetPointData()->AddArray(liqDisplacement);

	writer->SetFileName(filename.c_str());
	writer->SetInputData(polyData);
	
	writer->Write();
}

void MCLiqLattice::FromVTK(int frame)
{
	char buf[128];
	sprintf(buf, "%05d", frame);
	
	std::string filename = outputDir + "/liq" + buf + ".vtp";
	
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(filename.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	vtkDataArray* dispData = polyData->GetPointData()->GetArray("Displacement");

	nLiq = std::floor(Ntot*Ldens);
	vtkIdType Npoints = polyData->GetNumberOfPoints();
		
	if ( Npoints != nLiq )
		throw std::runtime_error("MCLiqLattice: Found liquid configuration file with incompatible dimension " + std::to_string(Npoints));
	else
		std::cout << "Starting from liquid configuration file " << filename << std::endl;
	
	for ( int i = 0; i < nLiq; i++ )
	{
		disp initDisp;
		double point[3];
		
		polyData->GetPoint(i, point);
		
		initDisp.dx = dispData->GetComponent(i, 0);
		initDisp.dy = dispData->GetComponent(i, 1);
		initDisp.dz = dispData->GetComponent(i, 2);

		int ixp = (int) 1*point[0];
		int iyp = (int) 2*point[1];
		int izp = (int) 4*point[2];
		
		int idx = ixp + iyp*L + izp*L2;
		
		spinTable[idx]++;
		spinIdTable[idx] = i;
		
		spinConf.push_back(idx);
		spinDisp.push_back(initDisp);
	}
}
