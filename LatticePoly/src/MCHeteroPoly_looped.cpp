//
//  MCHeteroPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <sstream>

#include "MCHeteroPoly_looped.hpp"


MCHeteroPoly_looped::MCHeteroPoly_looped(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCHeteroPoly_looped::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);


	if ( !RestartFromFile )
	{
		std::ifstream domainFile(InsulatorPath);
		
		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + InsulatorPath);
		
		std::string line;
		int line_id=0;
		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			int d2;

			
			
			if ( ss >> d1 >>d2)
			{
			tadConf[d1].insulator_type.push_back(d2);
			std::cout << "at monomer " << d1 << " type "<< d2<<  std::endl;
				
				
				
			++line_id;
			}
			
		}
		
		domainFile.close();
	}
	
	MCHeteroPoly_looped::BuildLoopTable();
	
}



void MCHeteroPoly_looped::BuildLoopTable()
{
	
	for ( int vi = 0; vi < Ntot; ++vi )
		for ( int k = 0; k < 9; ++k )
			hetTable_insulator[k][vi] = 0;
	
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		if ( tad->insulator_type.size() != 0 )
			{
			for ( int type_id = 0; type_id< tad->insulator_type.size();  ++type_id )
			{
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					++hetTable_insulator[tad->insulator_type.at(type_id)][vi];
				}
			}
		}
	}
}
void MCHeteroPoly_looped::AcceptMove()
{

	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->insulator_type.size() != 0 )
	{
		for ( int type_id = 0; type_id< tadTrial-> insulator_type.size();  ++type_id )
		{
			for ( int v = 0; v < 13; ++v )
			{
			
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				
				--hetTable_domain[tadTrial->insulator_type.at(type_id)][vi1];
				++hetTable_domain[tadTrial->insulator_type.at(type_id)][vi2];

			}
		}
	}
}
	
	

double MCHeteroPoly_looped::GetEffectiveEnergy() const
{
	double E_copolymer = 0;
	
	double cis_insulators_energy=0;
	double trans_insulators_energy=0;
	
	if ( Jaa > 0. or Jbb > 0. or Jtad_b > 0. or Jtad_a > 0. )
	{
		double domain_energy=0;
		double inter_domain_energy=0;
		double tads_energy=0;
		double J_domain = tadTrial->domain==1 ? Jaa : Jbb;
		double J_tad =  tadTrial->domain==1 ? Jtad_a : Jtad_b;
		
		if(tadTrial->domain>=0) //energy withn domain 0 B or 1 A
			domain_energy= J_domain * (hetTable_domain[tadTrial->domain][tadUpdater->vo]-hetTable_domain[tadTrial->domain][tadUpdater->vn]);
		if(tadTrial->domain==0)//energy between  AB (if I move B)
			inter_domain_energy= Jab * (hetTable_domain[1][tadUpdater->vo]-hetTable_domain[1][tadUpdater->vn]);
		if(tadTrial->domain==1)//energy between  AB (if I move A)
			inter_domain_energy= Jab * (hetTable_domain[0][tadUpdater->vo]-hetTable_domain[1][tadUpdater->vn]);
		
		if ( tadTrial->type != -1 )//energy within tad
			tads_energy = J_tad * (hetTable_tads[tadTrial->type][tadUpdater->vo]-hetTable_tads[tadTrial->type][tadUpdater->vn]);
		
		E_copolymer=domain_energy+inter_domain_energy+tads_energy;
	}

	if (  J_insulator_cis>0 or J_insulator_trans>0)
	{

		
		if ( tadTrial->insulator_type.size() != 0 )
			for ( int type_id = 0; type_id< tadTrial-> insulator_type.size();  ++type_id )
					cis_insulators_energy =cis_insulators_energy + J_insulator_cis * (hetTable_insulator[tadTrial-> insulator_type.size()][tadUpdater->vo]-hetTable_insulator[tadTrial-> insulator_type.size()][tadUpdater->vn]);

		if(J_insulator_trans>0) //energy between insulators of diffent kinds
		{
			for ( int k = 0; k < 9; ++k )//Number of insulators
				if(std::find( tadTrial->insulator_type.begin(), tadTrial->insulator_type.end(),k) ==  tadTrial->insulator_type.end())
					trans_insulators_energy = trans_insulators_energy+J_insulator_trans*(hetTable_insulator[k][tadUpdater->vo]-hetTable_insulator[k][tadUpdater->vn]);
		}
	}
	
	return E_copolymer+cis_insulators_energy+trans_insulators_energy;
}


double MCHeteroPoly_looped::GetCouplingEnergy(const int spinTable[Ntot]) const
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




vtkSmartPointer<vtkPolyData> MCHeteroPoly_looped::GetVTKData()
{
	vtkSmartPointer<vtkPolyData> polyData = MCHeteroPoly::GetVTKData();
	
	auto insulator = vtkSmartPointer<vtkIntArray>::New();
	
	insulator->SetName("insulator");
	insulator->SetNumberOfComponents(1);
	
	
	
	
	for ( int t = 0; t < Ntad; ++t )
	{
		if(tadConf[t].insulator_type.size()!=0)
			insulator->InsertNextValue(1);
		else
			insulator->InsertNextValue(0);
			
	}
	
	
	polyData->GetPointData()->AddArray(insulator);
	
	
	
	return polyData;
}

void MCHeteroPoly_looped::SetVTKData(const vtkSmartPointer<vtkPolyData> polyData)
{
	MCHeteroPoly::SetVTKData(polyData);//Does not work
}
