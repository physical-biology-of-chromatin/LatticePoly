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


MCHeteroPoly::MCHeteroPoly(MCLattice* _lat): MCPoly(_lat) {};

void MCHeteroPoly::Init(int Ninit)
{
	MCPoly::Init(Ninit);
	
	for ( int vi = 0; vi < Ntot; ++vi )
		tadHetTable[vi] = 0;
	
	if ( RestartFromFile )
	{
		for ( int t = 0; t < Nchain; ++t )
		{
			if ( tadType[t] == 1 )
				++tadHetTable[tadConf[t]];
		}
	}
	
	else
	{
		std::ifstream domainFile(domainPath);
		
		std::string line;
		std::vector<std::pair<int, int>> domains;

		if ( !domainFile.good() )
			throw std::runtime_error("MCHeteroPoly: Couldn't open file " + domainPath);
		
		while ( std::getline(domainFile, line) )
		{
			int d1, d2;
			std::istringstream ss(line);

			if ( ss >> d1 >> d2 )
			{
				if ( (d1 >= 0) && (d2 >= 0) && (d1 < Nchain) && (d2 < Nchain) )
					domains.push_back((d1 < d2) ? std::make_pair(d1, d2) : std::make_pair(d2, d1));
				else
					throw std::runtime_error("MCHeteroPoly: Found inconsistent domain boundaries " + std::to_string(d1) + "-" + std::to_string(d2));
			}
			
			else
				throw std::runtime_error("MCHeteroPoly: Bad line " + line + " in file " + domainPath);
		}

		for ( auto it = domains.begin(); it != domains.end(); ++it )
		{
			for ( int i = it->first; i <= it->second; ++i )
			{
				tadType[i] = 1;
				tadHetTable[tadConf[i]] = 1;
			}
		}
	}
}

void MCHeteroPoly::AcceptMove()
{
	MCPoly::AcceptMove();
	
	if ( tadType[tad->n] == 1 )
	{
		--tadHetTable[tad->vo];
		++tadHetTable[tad->vn];
	}
}

double MCHeteroPoly::GetSpecificEnergy() const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int v = 0; v < 13; ++v )
		{
			if ( v == 0 )
			{
				E1 -= tadHetTable[tad->vo] - 1.;
				E2 -= tadHetTable[tad->vn];
			}
			
			else
			{
				E1 -= tadHetTable[lat->bitTable[v][tad->vo]];
				E2 -= tadHetTable[lat->bitTable[v][tad->vn]];
			}
		}
	}
	
	return Jpp * (E2-E1);
}

double MCHeteroPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tad->vo : lat->bitTable[v][tad->vo];
			int vi2 = (v == 0) ? tad->vn : lat->bitTable[v][tad->vn];
			
			E1 -= spinTable[vi1];
			E2 -= spinTable[vi2];
		}
	}
	
	return Jlp * (E2-E1);
}
