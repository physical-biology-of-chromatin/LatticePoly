//
//  MCHeteroPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <sstream>

#include "MCHeteroPoly.hpp"


MCHeteroPoly::MCHeteroPoly(MCLattice* _lat): MCPoly(_lat) {}

void MCHeteroPoly::Init(int Ninit)
{
	MCPoly::Init(Ninit);
	
	for ( int vi = 0; vi < Ntot; ++vi )
		hetTable[vi] = 0;
	
	if ( !RestartFromFile )
	{
		std::ifstream domainFile(domainPath);

		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + domainPath);
				
		std::string line;

		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			int d2;
			
			if ( ss >> d1 >> d2 )
			{
				if ( (d1 >= 0) && (d2 >= 0) && (d1 <= Nchain) && (d2 <= Nchain) )
				{
					for ( int t = std::min(d1, d2); t < std::max(d1, d2); ++t )
						tadConf[t].type = 1;
				}
				
				else
					throw std::runtime_error("MCHeteroPoly: Found inconsistent domain boundaries '" + line + "' in file " + domainPath);
			}
			
			else
				throw std::runtime_error("MCHeteroPoly: Bad line '" + line + "' in file " + domainPath);
		}
		
		domainFile.close();
	}
		
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->type == 1 )
		{
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
				
				++hetTable[vi];
			}
		}
	}
}

void MCHeteroPoly::AcceptMove()
{
	for ( int v = 0; v < Ntot; ++v )
		if ( hetTable[v] == -1 )
			std::cout <<"het -1" << std::endl;
	
	
	MCPoly::AcceptMove();
	
	if ( tadTrial->type == 1 )
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			--hetTable[vi1];
			++hetTable[vi2];
		}
	}
		
}

double MCHeteroPoly::GetEffectiveEnergy() const
{
	if ( Jpp > 0. )
	{
		if ( tadTrial->type == 1 )
			return Jpp * (hetTable[tadUpdater->vo]-hetTable[tadUpdater->vn]);
	}
	
	return 0.;
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( tadTrial->type == 1 )
		{
			double dE = 0.;
		
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
				dE += spinTable[vi1];
				dE -= spinTable[vi2];
			}
		
			return Jlp * dE;
		}
	}
	
	return 0.;
}


vtkSmartPointer<vtkPolyData> MCHeteroPoly::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCPoly::GetVTKData();
	
	auto type = vtkSmartPointer<vtkIntArray>::New();

	type->SetName("TAD type");
	type->SetNumberOfComponents(1);
	
	for ( int t = 0; t < Ntad; ++t )
		type->InsertNextValue(tadConf[t].type);
		
	polyData->GetPointData()->AddArray(type);

	return polyData;
}

void MCHeteroPoly::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	MCPoly::SetVTKData(polyData);
	
	vtkDataArray* type = polyData->GetPointData()->GetArray("TAD type");

	for ( int t = 0; t < Ntad; ++t )
		tadConf[t].type = (int) type->GetComponent(t, 0);
}
