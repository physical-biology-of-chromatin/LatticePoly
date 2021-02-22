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
		std::vector<std::pair<int, int>> domains;

		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + domainPath);
		
		while ( std::getline(domainFile, line) )
		{
			std::istringstream ss(line);
			
			int d1;
			int d2;
			
			if ( ss >> d1 >> d2 )
			{
				if ( (d1 >= 0) && (d2 >= 0) && (d1 < Nchain) && (d2 < Nchain) )
					domains.push_back((d1 < d2) ? std::make_pair(d1, d2) : std::make_pair(d2, d1));
				else
					throw std::runtime_error("MCHeteroPoly: Found inconsistent domain boundaries '" + line + "' in file " + domainPath);
			}
			
			else
				throw std::runtime_error("MCHeteroPoly: Bad line '" + line + "' in file " + domainPath);
		}

		for ( auto it = domains.begin(); it != domains.end(); ++it )
		{
			for ( int t = it->first; t <= it->second; ++t )
				tadConf[t].type = 1;
		}
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
void MCHeteroPoly::OriginMove()
{}
void MCHeteroPoly::ForkMove()
{}
