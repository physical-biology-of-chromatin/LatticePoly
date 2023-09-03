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
	{
		hetNeighborhood[vi] = 0;
		hetTable[vi] = 0;
	}

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
			hetTable[tad->pos] += tad -> type;

			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
					
				hetNeighborhood[vi] += tad->type;
				
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

		hetTable[tadUpdater->vo] -= tadTrial -> type;
		hetTable[tadUpdater->vn] += tadTrial -> type;
	
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			//std::cout <<"before "<<hetNeighborhood[vi1]<< std::endl;
			hetNeighborhood[vi1] -= tadTrial->type;
			hetNeighborhood[vi2] += tadTrial->type;
			//std::cout <<"after "<<hetNeighborhood[vi1]<< std::endl;
		}
	}

	//double tothet=0;
	//for ( int v = 0; v < Ntot; ++v )
	//	tothet+=hetNeighborhood[v];
	//std::cout <<tothet << std::endl;
}

double MCHeteroPoly::GetEffectiveEnergy() const
{
	if ( Jpp > 0. )
	{
		if ( tadTrial->type != 0 )
		{
			//std::cout <<Jpp* (hetNeighborhood[tadUpdater->vo]-hetNeighborhood[tadUpdater->vn]) * tadTrial->type << std::endl;
			return Jpp * (hetNeighborhood[tadUpdater->vo]-hetNeighborhood[tadUpdater->vn])* tadTrial->type;
		}	
			
	}
	
	return 0.;
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot], const int spinNeighborhood[Ntot]) const
{
	if ( Jlp > 0. )
	{
		if ( tadTrial->type != 0 ) 
		{
			double dN = 0.;
		
			for ( int v = 1; v < 13; ++v )
			{	

				if ( lat->bitTable[v][tadUpdater->vo] != tadUpdater->vn && spinTable[lat->bitTable[v][tadUpdater->vo]] != 0.)
				{
					int c = 0;

					for ( int i = 0; i<13; ++i)
					{
						if ( lat->bitTable[v][tadUpdater->vo] == lat->bitTable[i][tadUpdater->vn])
						{
							c = 1;
						}
					}
					if (c == 0)
					{
						dN += ((Jlp_Valency < hetNeighborhood[lat->bitTable[v][tadUpdater->vo]]) ? 0. : 1.);
					}
				}

				if (lat->bitTable[v][tadUpdater->vn] != tadUpdater->vo && spinTable[lat->bitTable[v][tadUpdater->vn]] != 0.)
				{

					int c = 0;

					for ( int i = 0; i<13; ++i)
					{
						if ( lat->bitTable[v][tadUpdater->vn] == lat->bitTable[i][tadUpdater->vo])
						{
							c = 1;
						}
					}

					if (c == 0)
					{
						dN -= ((Jlp_Valency < hetNeighborhood[lat->bitTable[v][tadUpdater->vn]]+1) ? 0. : 1.);
					}
				}
			}

			return Jlp / 2 * (((Jpl_Valency < spinNeighborhood[tadUpdater->vo]) ? Jpl_Valency : spinNeighborhood[tadUpdater->vo]) - ((Jpl_Valency < spinNeighborhood[tadUpdater->vn]) ? Jpl_Valency : spinNeighborhood[tadUpdater->vn]) + dN);

		}
	}
	
	return 0.;
}
