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
#include <vtkXMLPolyDataWriter.h>

#include "MCLiqLattice.hpp"


void MCLiqLattice::Init(std::mt19937_64& rngEngine)
{
	if ( Qcg != 1 )
		throw std::runtime_error("MCLiqLattice: Must set Qcg to 1 for microscopic liquid simulations");
	
	MCCGLattice::Init(rngEngine);

	for ( int i = 0; i < Ntot; i++ )
	{
		if ( spinTable[i] > 0 )
		{
			spinIdTable[i] = (int) spinConf.size();
			spinConf.push_back(i);
		}
		
		else
			spinIdTable[i] = -1;
	}
}

void MCLiqLattice::TrialMove(std::mt19937_64& rngEngine, double* dE)
{
	int v1 = rngEngine() % nLiq;
	int v2 = rngEngine() % 12;

	idx1 = spinConf[v1];
	idx2 = bitTable[v2+1][idx1];
	
	*dE = GetSpinEnergy();	
}

void MCLiqLattice::AcceptMove()
{
	int spin1 = spinTable[idx1];
	int spin2 = spinTable[idx2];
	
	if ( (spin1 == 1) || (spin2 == 1) )
	{
		DisplaceSpins();

		spinTable[idx1] = spin2;
		spinTable[idx2] = spin1;
						
		if ( (spin1 == 1) && (spin2 == 0) )
		{
			int id1 = spinIdTable[idx1];
		
			spinIdTable[idx1] = -1;
			spinIdTable[idx2] = id1;
			
			spinConf[id1] = idx2;
		}
		
		else if ( (spin1 == 0) && (spin2 == 1) )
		{
			int id2 = spinIdTable[idx2];
			
			spinIdTable[idx1] = id2;
			spinIdTable[idx2] = -1;
				
			spinConf[id2] = idx1;
		}
		
		else
		{
			int id1 = spinIdTable[idx1];
			int id2 = spinIdTable[idx2];
			
			spinIdTable[idx1] = id2;
			spinIdTable[idx2] = id1;
			
			spinConf[id1] = idx2;
			spinConf[id2] = idx1;
		}
	}
}

void MCLiqLattice::DisplaceSpins()
{	
	double dx = xyzTable[0][idx2] - xyzTable[0][idx1];
	double dy = xyzTable[1][idx2] - xyzTable[1][idx1];
	double dz = xyzTable[2][idx2] - xyzTable[2][idx1];
	
	if ( std::abs(dx) > L/2. ) dx -= std::copysign(L, dx);
	if ( std::abs(dy) > L/2. ) dy -= std::copysign(L, dy);
	if ( std::abs(dz) > L/2. ) dz -= std::copysign(L, dz);
	
	if ( spinTable[idx1] == 1 )
	{
		int id1 = spinIdTable[idx1];

		spinDisp[id1].dx += dx;
		spinDisp[id1].dy += dy;
		spinDisp[id1].dz += dz;
	}
	
	if ( spinTable[idx2] == 1 )
	{
		int id2 = spinIdTable[idx2];

		spinDisp[id2].dx -= dx;
		spinDisp[id2].dy -= dy;
		spinDisp[id2].dz -= dz;
	}
}

double MCLiqLattice::GetSpinEnergy() const
{
	double dE = 0.;
	
	if ( spinTable[idx1] != spinTable[idx2] )
	{
		for ( int i = 0; i < 12; i++ )
		{
			if ( bitTable[i+1][idx1] != idx2 )
				dE += spinTable[bitTable[i+1][idx1]];
		
			if ( bitTable[i+1][idx2] != idx1 )
				dE -= spinTable[bitTable[i+1][idx2]];
		}
	
		dE *= Jll * (spinTable[idx1]-spinTable[idx2]);
	}
	
	return dE;
}

double MCLiqLattice::GetCouplingEnergy(const int tadHetTable[Ntot]) const
{
	if ( spinTable[idx1] != spinTable[idx2] )
		return MCCGLattice::GetCouplingEnergy(tadHetTable);
	
	return 0.;
}

double MCLiqLattice::GetSpinDensity(int idx) const
{
	double aveLiq = 0.;

	for ( int j = 0; j < 12; j++ )
		aveLiq += spinTable[bitTable[j+1][idx]] / 12.;
		
	return aveLiq;
}
