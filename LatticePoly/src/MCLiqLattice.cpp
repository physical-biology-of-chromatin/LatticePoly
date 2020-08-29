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
	MCLattice::Init(Ninit);

	for ( int vi = 0; vi < Ntot; ++vi )
	{
		spinTable[vi] = 0;
		lookupTable[vi] = -1;
	}
	
	if ( RestartFromFile )
		FromVTK(Ninit);

	else
	{
		nLiq = 0;
		
		if ( InitDrop )
			GenerateDroplets();
		else
			GenerateRandom();
		
		for ( int i = 0; i < nLiq; ++i )
		{
			disp initDisp;
			
			initDisp.dx = 0.;
			initDisp.dy = 0.;
			initDisp.dz = 0.;

			spinDisp.push_back(initDisp);
		}
		
		for ( int vi = 0; vi < Ntot; ++vi )
		{
			if ( spinTable[vi] > 0 )
			{
				lookupTable[vi] = (int) spinConf.size();
				
				spinConf.push_back(vi);
			}
		}
	}
	
	std::cout << "Set up lattice with fixed liquid density " << nLiq / ((double) Ntot) << std::endl;
}

void MCLiqLattice::GenerateDroplets()
{
	int centers[3][Ndrop];
	
	int r = std::floor(R) + 1; // Set to 1 to allow initial droplets to cross PBCs
	
	for ( int i = 0; i < Ndrop; ++i )
	{
		centers[0][i] = (rngEngine() % (L-2*r+1)) + r;
		centers[1][i] = (rngEngine() % (L-2*r+1)) + r;
		centers[2][i] = (rngEngine() % (L-2*r+1)) + r;
	}
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		for ( int i = 0; i < Ndrop; ++i )
		{
			double dx = xyzTable[0][vi] - centers[0][i];
			double dy = xyzTable[1][vi] - centers[1][i];
			double dz = xyzTable[2][vi] - centers[2][i];

			if      ( SQR(dx-L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx-L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx)   + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx)   + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy-L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;

			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy)   + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz-L) < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz)   < SQR(R) ) spinTable[vi] = 1;
			else if ( SQR(dx+L) + SQR(dy+L) + SQR(dz+L) < SQR(R) ) spinTable[vi] = 1;
			
			if ( spinTable[vi] > 0 )
			{
				++nLiq;
				break;
			}
		}
	}
}

void MCLiqLattice::GenerateRandom()
{
	while ( nLiq < std::floor(Ntot*Ldens) )
	{
		int vi = rngEngine() % Ntot;
		
		if ( spinTable[vi] == 0 )
		{
			++spinTable[vi];
			++nLiq;
		}
	}
}

void MCLiqLattice::TrialMove(double* dE)
{
	int n = rngEngine() % nLiq;
	int v = rngEngine() % 12;

	v1 = spinConf[n];
	v2 = bitTable[v+1][v1];
	
	*dE = GetSpinEnergy();	
}

void MCLiqLattice::AcceptMove()
{
	DisplaceSpins();
					
	if ( spinTable[v2] == 0 )
	{
		int id1 = lookupTable[v1];
	
		lookupTable[v1] = -1;
		lookupTable[v2] = id1;
		
		spinConf[id1] = v2;
	}
	
	else
	{
		int id1 = lookupTable[v1];
		int id2 = lookupTable[v2];
		
		lookupTable[v1] = id2;
		lookupTable[v2] = id1;
		
		spinConf[id1] = v2;
		spinConf[id2] = v1;
	}
	
	spinTable[v1] = spinTable[v2];
	spinTable[v2] = 1;
}

void MCLiqLattice::DisplaceSpins()
{	
	double dx = xyzTable[0][v2] - xyzTable[0][v1];
	double dy = xyzTable[1][v2] - xyzTable[1][v1];
	double dz = xyzTable[2][v2] - xyzTable[2][v1];
	
	if ( std::abs(dx) > L/2. ) dx -= std::copysign(L, dx);
	if ( std::abs(dy) > L/2. ) dy -= std::copysign(L, dy);
	if ( std::abs(dz) > L/2. ) dz -= std::copysign(L, dz);
	
	int id1 = lookupTable[v1];

	spinDisp[id1].dx += dx;
	spinDisp[id1].dy += dy;
	spinDisp[id1].dz += dz;
	
	if ( spinTable[v2] == 1 )
	{
		int id2 = lookupTable[v2];

		spinDisp[id2].dx -= dx;
		spinDisp[id2].dy -= dy;
		spinDisp[id2].dz -= dz;
	}
}

double MCLiqLattice::GetSpinEnergy() const
{
	if ( Jll > 0. )
	{
		if ( spinTable[v2] == 0 )
		{
			double dE = 0.;

			for ( int v = 0; v < 12; ++v )
			{
				if ( bitTable[v+1][v1] != v2 )
					dE += spinTable[bitTable[v+1][v1]];
				
				if ( bitTable[v+1][v2] != v1 )
					dE -= spinTable[bitTable[v+1][v2]];
			}
			
			return Jll * dE;
		}
	}
	
	return 0.;
}

double MCLiqLattice::GetCouplingEnergy(const int hetTable[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( spinTable[v2] == 0 )
			return Jlp * (hetTable[v1]-hetTable[v2]);
	}
	
	return 0.;
}

void MCLiqLattice::ToVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "liq%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto liqDensity = vtkSmartPointer<vtkFloatArray>::New();
	auto liqDisplacement = vtkSmartPointer<vtkFloatArray>::New();
	
	liqDensity->SetName("Density");
	liqDensity->SetNumberOfComponents(1);
	
	liqDisplacement->SetName("Displacement");
	liqDisplacement->SetNumberOfComponents(3);
		
	for ( int i = 0; i < nLiq; ++i )
	{
		int vi = spinConf[i];
		double aveDensity = 0.;

		for ( int v = 0; v < 12; ++v )
			aveDensity += spinTable[bitTable[v+1][vi]] / 12.;
		
		double x = xyzTable[0][vi];
		double y = xyzTable[1][vi];
		double z = xyzTable[2][vi];
		
		double dx = spinDisp[i].dx;
		double dy = spinDisp[i].dy;;
		double dz = spinDisp[i].dz;;
				
		points->InsertNextPoint(x, y, z);
	
		liqDensity->InsertNextValue(aveDensity);
		liqDisplacement->InsertNextTuple3(dx, dy, dz);
	}
	
	auto polyData = vtkSmartPointer<vtkPolyData>::New();
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	polyData->SetPoints(points);
	
	polyData->GetPointData()->AddArray(liqDensity);
	polyData->GetPointData()->AddArray(liqDisplacement);

	writer->SetFileName(path.c_str());
	writer->SetInputData(polyData);
	
	writer->Write();
}

void MCLiqLattice::FromVTK(int frame)
{
	char fileName[32];
	sprintf(fileName, "liq%05d.vtp", frame);
	
	std::string path = outputDir + "/" + fileName;
	
	auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

	reader->SetFileName(path.c_str());
	reader->Update();
	
	vtkPolyData* polyData = reader->GetOutput();
	vtkDataArray* dispData = polyData->GetPointData()->GetArray("Displacement");

	nLiq = (int) polyData->GetNumberOfPoints();
	
	if ( (InitDrop == 0) && (nLiq != std::floor(Ntot*Ldens)) )
		throw std::runtime_error("MCLiqLattice: Found liquid configuration file with incompatible dimension " + std::to_string(nLiq));
	
	std::cout << "Starting from liquid configuration file " << path << std::endl;
	
	for ( int i = 0; i < nLiq; ++i )
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
		
		int vi = ixp + iyp*L + izp*L2;
		
		lookupTable[vi] = i;
		
		spinConf.push_back(vi);
		spinDisp.push_back(initDisp);
		
		++spinTable[vi];
	}
}
