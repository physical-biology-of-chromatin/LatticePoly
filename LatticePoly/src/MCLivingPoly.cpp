//
//  MCLivingPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 27/04/2021.
//  Copyright © 2021 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <sstream>

#include "MCLivingPoly.hpp"


MCLivingPoly::MCLivingPoly(MCLattice* _lat): MCHeteroPoly_looped(_lat) {}

void MCLivingPoly::Init(int Ninit,int chrom ,int chrom_pos[3])
{
	MCHeteroPoly_looped::Init(Ninit,chrom,chrom_pos);

	if ( !RestartFromFile )
	{
		if ( propagationMode != 2 )
		{
			for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
			{
				if ( tad->type == 1 )
				{
					double rnd = (propagationMode == 0) ? lat->rngDistrib(lat->rngEngine) : 0.;
					
					if ( rnd < inactiveRatio )
					{
						tad->type = 2;

						for ( int v = 0; v < 13; ++v )
						{
							int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
							
							--hetTable[vi];
						}
					}
				}
			}
		}
	}
	
	if ( propagationMode == 1 )
	{
		std::ifstream colorFile(colorPath);

		if ( !colorFile.good() )
			throw std::runtime_error("MCLivingPoly: Couldn't open file " + colorPath);
		
		std::string line;

		while ( std::getline(colorFile, line) )
		{
			int idxTad;
			
			std::vector<int> paintedIds;
			std::istringstream ss(line);
			
			while ( ss >> idxTad )
			{
				if ( idxTad < Nchain )
					paintedIds.push_back(idxTad);
				else
					throw std::runtime_error("MCLivingPoly: Color data in " + colorPath + " incompatible with chain dimensions");
			}
			
			colorData.push_back(paintedIds);
		}
		
		colorFile.close();
	}
}

void MCLivingPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	
	if ( propagationMode == 0 )
	{
		if ( tadTrial->type == 2 )
		{
			double rnd = lat->rngDistrib(lat->rngEngine);
			
			if ( rnd < propRate / ((double) Ninter*Nmeas) )
			{
				tadTrial->type = 1;
				
				for ( int v = 0; v < 13; ++v )
				{
					int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
					
					++hetTable[vi];
				}
			}
		}
	}
}

void MCLivingPoly::UpdateFromFile(int idx)
{
	if ( idx < (int) colorData.size() + Nrelax )
	{
		std::vector<int> paintedIds = colorData[idx-Nrelax];

		for ( auto idxTad = paintedIds.begin(); idxTad != paintedIds.end(); ++idxTad )
		{
			if ( *idxTad >= 0 )
			{
				MCTad* tad = &tadConf[*idxTad];
				
				if ( tad->type == 2 )
				{
					tad->type = 1;
					
					for ( int v = 0; v < 13; ++v )
					{
						int vi = (v == 0) ? tad->pos : lat->bitTable[v][tad->pos];
						
						++hetTable[vi];
					}
				}
			}
		}
	}
}

void MCLivingPoly::ToVTK(int frame,std::string number)
{
	MCPoly::ToVTK(frame,number);
	
	if ( propagationMode == 1 )
		UpdateFromFile(frame);
}
