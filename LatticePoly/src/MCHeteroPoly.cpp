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
	

	
	if ( !RestartFromFile )
	{
		std::ifstream domainFile(domainPath);

		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + domainPath);
				
		std::string line;
		int line_id=0;
		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			int d2;

			if ( ss >> d1 >>d2 )
			{
				
				tadConf[line_id].type = d2;
				std::cout << "at monomer " << line_id << " type "<< d2<< std::endl;

				tadConf[line_id].domain = d1;
				std::cout << "at monomer " << line_id << " domain "<< d1<< std::endl;

				++line_id;
			}
			
		}
		
		domainFile.close();
	}
	
	MCHeteroPoly::BuildHetTable();

	
}



void MCHeteroPoly::BuildHetTable()
{
	for ( int vi = 0; vi < Ntot; ++vi )
		for ( int k = 0; k < 2; ++k )
			hetTable_domain[k][vi] = 0;
	
	for ( int vi = 0; vi < Ntot; ++vi )
		for ( int k = 0; k < 4; ++k )
			hetTable_tads[k][vi] = 0;
	
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->type != 0 )
		{
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
				
				++hetTable_tads[tad->type-1][vi];
			}
		}

		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
				
			++hetTable_domain[tad->domain][vi];
		}		
	}
	

}
void MCHeteroPoly::AcceptMove()
{

	MCPoly::AcceptMove();
	
	if ( tadTrial->type != 0 )
	{

		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			--hetTable_tads[tadTrial->type-1][vi1];
			++hetTable_tads[tadTrial->type-1][vi2];
		}
		
	}
	

	for ( int v = 0; v < 13; ++v )
	{

		int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
		int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
		
		--hetTable_domain[tadTrial->domain][vi1];
		++hetTable_domain[tadTrial->domain][vi2];
	}
	

		
}

double MCHeteroPoly::GetEffectiveEnergy() const
{
	if ( Jaa > 0. or Jbb > 0. or Jtad_b > 0. or Jtad_a > 0. )
	{
		double inter_domain_energy=0;
		double J_domain = tadTrial->domain==0 ? Jaa : Jbb;
		double J_tad =  tadTrial->domain==0 ? Jtad_a : Jtad_b;
		double domain_energy= J_domain * (hetTable_domain[tadTrial->domain][tadUpdater->vo]-hetTable_domain[tadTrial->domain][tadUpdater->vn]);
		if(tadTrial->domain==0)
			 inter_domain_energy= Jab * (hetTable_domain[1][tadUpdater->vo]-hetTable_domain[1][tadUpdater->vn]);
		else
			 inter_domain_energy= Jab * (hetTable_domain[0][tadUpdater->vo]-hetTable_domain[1][tadUpdater->vn]);


		if ( tadTrial->type != 0 )
			return inter_domain_energy+domain_energy+J_tad * (hetTable_tads[tadTrial->type-1][tadUpdater->vo]-hetTable_tads[tadTrial->type-1][tadUpdater->vn]);
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
	auto domain = vtkSmartPointer<vtkIntArray>::New();

	type->SetName("TAD type");
	type->SetNumberOfComponents(1);
	
	domain->SetName("Domain type");
	domain->SetNumberOfComponents(1);

	
	for ( int t = 0; t < Ntad; ++t )
	{
		type->InsertNextValue(tadConf[t].type);
		domain->InsertNextValue(tadConf[t].domain);
	}

		
	polyData->GetPointData()->AddArray(type);
	polyData->GetPointData()->AddArray(domain);


	return polyData;
}

void MCHeteroPoly::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	MCPoly::SetVTKData(polyData);
	
	vtkDataArray* type = polyData->GetPointData()->GetArray("TAD type");
	vtkDataArray* domain = polyData->GetPointData()->GetArray("Domain type");


	for ( int t = 0; t < Ntad; ++t )
	{
		tadConf[t].type = (int) type->GetComponent(t, 0);
		tadConf[t].domain = (int) domain->GetComponent(t, 0);

	}
}
