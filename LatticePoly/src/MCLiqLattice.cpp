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

	nLiq = 0;
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		spinTable[vi] = 0;
		lookupTable[vi] = -1;
	}
	
	if ( RestartFromFile )
		FromVTK(Ninit);

	else
	{
		if ( InitDrop )
			GenerateDroplets();
		else
			GenerateRandom();
		
		spinConf.resize(nLiq);
		spinDisp.resize(nLiq);

		std::fill(spinDisp.begin(), spinDisp.end(), (double3) {0., 0., 0.});

		int ctr = 0;
		
		for ( int vi = 0; vi < Ntot; ++vi )
		{
			if ( spinTable[vi] > 0 )
			{
				lookupTable[vi] = ctr;
				spinConf[ctr] = vi;
				
				++ctr;
			}
		}
	}
	
	std::cout << "Set up lattice with fixed liquid density " << nLiq / ((double) Ntot) << std::endl;
}

void MCLiqLattice::GenerateDroplets()
{
	std::vector<double3> centers(Ndrop);
	
	int r = std::floor(R) + 1; // Set to 1 to allow initial droplets to cross PBCs
	
	for ( int i = 0; i < Ndrop; ++i )
	{
		centers[i][0] = (rngEngine() % (L-2*r+1)) + r;
		centers[i][1] = (rngEngine() % (L-2*r+1)) + r;
		centers[i][2] = (rngEngine() % (L-2*r+1)) + r;
	}
	
	for ( int vi = 0; vi < Ntot; ++vi )
	{
		for ( int i = 0; i < Ndrop; ++i )
		{
			double dx = xyzTable[0][vi] - centers[i][0];
			double dy = xyzTable[1][vi] - centers[i][1];
			double dz = xyzTable[2][vi] - centers[i][2];

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
	for ( int i = 0; i < 3; ++i )
	{
		double disp = xyzTable[i][v2] - xyzTable[i][v1];

		if ( std::abs(disp) > L/2. )
			disp -= std::copysign(L, disp);

		int id1 = lookupTable[v1];
		spinDisp[id1][i] += disp;
		
		if ( spinTable[v2] == 1 )
		{
			int id2 = lookupTable[v2];
			spinDisp[id2][i] -= disp;
		}
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

double MCLiqLattice::GetCouplingEnergy(const double hetTable[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( spinTable[v2] == 0 )
			return Jlp * (hetTable[v1]-hetTable[v2]);
	}
	
	return 0.;
}

double MCLiqLattice::GetCouplingEnergyPainter(const double hetTable[Ntot], const double painterTable[Ntot] ) const
{

	double dE = MCLiqLattice::GetCouplingEnergy( hetTable );
	

	if ( ( Jlpp > 0. ) )
	{
		if ( spinTable[v2] == 0 )
			dE += Jlpp * (painterTable[v1]-painterTable[v2]);
	}	
    
	if ( ( EV > 0. ) )
	{
	        if ( spinTable[v2] == 0)
		       dE += EV * (bitTable[0][v2]-bitTable[0][v1]);
	}	

	return dE;
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
		
		double dx = spinDisp[i][0];
		double dy = spinDisp[i][1];
		double dz = spinDisp[i][2];
				
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
	
	spinConf.reserve(nLiq);
	spinDisp.reserve(nLiq);

	if ( (InitDrop == 0) && (nLiq != std::floor(Ntot*Ldens)) )
		throw std::runtime_error("MCLiqLattice: Found liquid configuration file with incompatible dimension " + std::to_string(nLiq));
	
	std::cout << "Starting from liquid configuration file " << path << std::endl;
	
	for ( int i = 0; i < nLiq; ++i )
	{
		double3 initDisp;
		double point[3];
		
		polyData->GetPoint(i, point);
		
		for ( int j = 0; j < 3; ++j )
			initDisp[j] = dispData->GetComponent(i, j);

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
