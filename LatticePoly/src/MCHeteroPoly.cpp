//
//  MCHeteroPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <utility>

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
		
		std::string line;
		std::vector<std::pair<int, double>> domains;

		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + domainPath);
		
		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			double dcharge;
            
			if ( ss >> d1 >> dcharge )
			{
				if ( d1 < Nchain )
				{
                    tadConf[d1].type = dcharge;
					//std::cout <<"het "<<tadConf[d1].type  << std::endl;
				}	
				else
					throw std::runtime_error("MCHeteroPoly: Found inconsistent values greater than chain length '" + line + "' in file " + domainPath);
			}
            
			else
				throw std::runtime_error("MCHeteroPoly: Bad line '" + line + "' in file " + domainPath);
        }
		domainFile.close();
	}
		
	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
			
		if ( tad->type != 0 )
		{	
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
				hetTable[vi] += tad->type;
				
			}
		}	
		
	}

	//for ( int v = 0; v < Nchain; ++v )
	//	std::cout <<tadConf[v].type << std::endl;
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
			
			//std::cout <<"before "<<hetTable[vi1]<< std::endl;
			hetTable[vi1] -= tadTrial->type;
			hetTable[vi2] += tadTrial->type;
			//std::cout <<"after "<<hetTable[vi1]<< std::endl;
		}
	}

	//double tothet=0;
	//for ( int v = 0; v < Ntot; ++v )
	//	tothet+=hetTable[v];
	//std::cout <<tothet << std::endl;
}

double MCHeteroPoly::GetEffectiveEnergy() const
{
	if ( Jpp > 0. )
	{
		if ( tadTrial->type != 0 )
		{
			//std::cout <<Jpp* (hetTable[tadUpdater->vo]-hetTable[tadUpdater->vn]) * tadTrial->type << std::endl;
			return Jpp * (hetTable[tadUpdater->vo]-hetTable[tadUpdater->vn])* tadTrial->type;
		}	
			
	}
	
	return 0.;
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( tadTrial->type != 0 ) 
		{
			double dE = 0.;
		
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
				dE += spinTable[vi1];
				dE -= spinTable[vi2];
			}
		
			return Jlp * dE * tadTrial->type;
		}
	}
	
	return 0.;
}
