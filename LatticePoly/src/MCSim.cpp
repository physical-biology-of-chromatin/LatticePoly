//
//  MCSim.cpp
//  LatticePoly
//
//  Created by mtortora on 29/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCSim.hpp"


template<class lattice, class polymer>
MCSim<lattice, polymer>::MCSim()
{
	lat = new lattice;
	pol = new polymer(lat);
}

template<class lattice, class polymer>
MCSim<lattice, polymer>::~MCSim()
{
	delete lat;
	delete pol;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Init()
{
	step = 0;
	
	acceptCountLiq = 0;
	acceptCountPoly = 0;

	InitRNG();
		
	lat->Init(rngEngine);
	pol->Init(rngEngine);
	
	DumpVTK(0);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::InitRNG()
{
	int seed;
	
	FILE* tmp = fopen("/dev/urandom", "rb");
	
	if ( (tmp != NULL) && (fread((void*) &seed, sizeof(seed), 1, tmp) != 0) )
		std::cout << "Using entropy-harvested random seed: " << seed << std::endl;
	
	else
	{
		seed = (int) time(NULL);
		
		std::cout << "Using system time as RNG seed: " << seed << std::endl;
	}
	
	fclose(tmp);

	rngEngine.seed(seed);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::PrintStats()
{
	unsigned long long nMovePoly = step * Nchain;
	double polyRate = acceptCountPoly / ((long double) nMovePoly);
	
	std::cout << "Polymer acceptance rate: " << 100*polyRate << "% (" << nMovePoly << " trial moves)" << std::endl;
	
	if ( latticeType == "MCLiqLattice" )
	{
		unsigned long long nMoveLiq = step * NliqMC;
		double liqRate = acceptCountLiq / ((long double) nMoveLiq);

		std::cout << "Liquid acceptance rate: " << 100*liqRate << "% (" << nMoveLiq << " trial moves)" << std::endl;
	}
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::DumpVTK(int idx)
{
	lat->ToVTK(idx);
	pol->ToVTK(idx);
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::Run()
{
	for ( int i = 0; i < Nchain; i++ )
		UpdateTAD();
	
	if ( latticeType == "MCLiqLattice" )
	{		
		for ( int i = 0; i < NliqMC; i++ )
			UpdateSpin();
		
		if ( step == Tbleach )
			static_cast<MCLiqLattice*>(lat)->BleachSpins();
	}
	
	step++;
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::UpdateTAD()
{
	bool acceptMove;
	double dEpol;
	
	pol->TrialMoveTAD(rngEngine, &dEpol);
	
	if ( pol->tad->legal )
	{
		acceptMove = MetropolisMove(dEpol);

		if ( acceptMove )
		{
			pol->AcceptMoveTAD();
			acceptCountPoly++;
		}
	}
}

template<class lattice, class polymer>
void MCSim<lattice, polymer>::UpdateSpin() {}

template<class lattice, class polymer>
bool MCSim<lattice, polymer>::MetropolisMove(double dE)
{
	if ( dE > 0. )
	{
		double rnd = rngDistrib(rngEngine);
		
		return (rnd < exp(-dE));
	}
	
	return true;
}

template<class lattice, class polymer>
bool MCSim<lattice, polymer>::ArrheniusMove(double dE, double Ebind)
{
	double rnd = rngDistrib(rngEngine);

	if ( dE > 0. )
		return (rnd < exp(Ebind-dE));
	
	return (rnd < exp(Ebind));
}

template<>
void MCSim<MCLattice, MCHeteroPoly>::UpdateTAD()
{
	bool acceptMove;
	double dEpol, dEspe;
	
	pol->TrialMoveTAD(rngEngine, &dEpol);
	
	if ( pol->tad->legal )
	{
		dEspe = (step < Trel) ? 0. : pol->GetSpecificEnergy();
		acceptMove = MetropolisMove(dEpol+dEspe);
		
		if ( acceptMove )
		{
			pol->AcceptMoveTAD();
			acceptCountPoly++;
		}
	}
}

template<>
void MCSim<MCLiqLattice, MCPoly>::UpdateSpin()
{
	bool acceptMove;
	double dElat;
	
	lat->TrialMoveSpin(rngEngine, &dElat);
	acceptMove = MetropolisMove(dElat);

	if ( acceptMove )
	{
		lat->AcceptMoveSpin();
		acceptCountLiq++;
	}
}

template<>
void MCSim<MCLiqLattice, MCHeteroPoly>::UpdateTAD()
{
	bool acceptMove;
	double dEpol;
	
	pol->TrialMoveSpinTAD(rngEngine, &dEpol);

	if ( pol->tad->legal )
	{
		double dEcpl = (step < Trel) ? 0. : pol->GetCouplingEnergy(lat->spinTable);
		acceptMove = MetropolisMove(dEpol+dEcpl);
		
		if ( acceptMove )
		{
			pol->AcceptMoveSpinTAD();			
			acceptCountPoly++;
		}
	}
}

template<>
void MCSim<MCLiqLattice, MCHeteroPoly>::UpdateSpin()
{
	bool acceptMove;
	double dElat;
	
	lat->TrialMoveSpin(rngEngine, &dElat);
		
	if ( Arrhenius )
	{
		double Ebind = (step < Trel) ? 0. : lat->GetBindingEnergy(pol->tadTable);
		acceptMove = ArrheniusMove(dElat, Ebind);
	}
	
	else
	{
		double dEcpl = (step < Trel) ? 0. : lat->GetCouplingEnergy(pol->tadTable);
		acceptMove = MetropolisMove(dElat+dEcpl);
	}

	if ( acceptMove )
	{
		lat->AcceptMoveSpin();
		acceptCountLiq++;
	}
}

template class MCSim<MCLattice, MCPoly>;
template class MCSim<MCLattice, MCHeteroPoly>;

template class MCSim<MCLiqLattice, MCPoly>;
template class MCSim<MCLiqLattice, MCHeteroPoly>;
