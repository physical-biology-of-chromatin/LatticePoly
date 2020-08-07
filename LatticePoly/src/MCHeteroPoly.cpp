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
	
	for ( int i = 0; i < Ntot; i++ )
		tadHetTable[i] = 0;
	
	if ( RestartFromFile )
	{
		for ( int i = 0; i < Nchain; i++ )
		{
			if ( tadType[i] == 1 )
				tadHetTable[tadConf[i]]++;
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
				if ( (d1 < Nchain) && (d2 < Nchain) )
					domains.push_back(std::make_pair(d1, d2));
				else
					throw std::runtime_error("MCHeteroPoly: Domain " + std::to_string(d1) + "-" + std::to_string(d2) + " incompatible with polymer size");
			}
			
			else
				throw std::runtime_error("MCHeteroPoly: Bad line " + line + " in file " + domainPath);
		}

		for ( auto it = domains.begin(); it != domains.end(); ++it )
		{
			int d1 = std::min(it->first, it->second);
			int d2 = std::max(it->first, it->second);
			
			for ( int i = d1; i <= d2; i++ )
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
		tadHetTable[tad->en]--;
		tadHetTable[tad->v2]++;
	}
}

double MCHeteroPoly::GetSpecificEnergy() const
{
	double E1 = 0.;
	double E2 = 0.;
	
	if ( tadType[tad->n] == 1 )
	{
		for ( int i = 0; i < 13; i++ )
		{
			if ( i == 0 )
			{
				E1 -= tadHetTable[tad->en] - 1.;
				E2 -= tadHetTable[tad->v2];
			}
			
			else
			{
				E1 -= tadHetTable[lat->bitTable[i][tad->en]];
				E2 -= tadHetTable[lat->bitTable[i][tad->v2]];
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
		for ( int i = 0; i < 13; i++ )
		{
			int v1 = (i == 0) ? tad->en : lat->bitTable[i][tad->en];
			int v2 = (i == 0) ? tad->v2 : lat->bitTable[i][tad->v2];
			
			E1 -= spinTable[v1];
			E2 -= spinTable[v2];
		}
	}
	
	return Jlp * (E2-E1);
}
