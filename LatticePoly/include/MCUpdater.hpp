//
//  MCUpdater.hpp
//  LatticePoly
//
//  Created by mtortora on 28/02/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCUpdater_hpp
#define MCUpdater_hpp

#include "MCHeteroPoly.hpp"


template<class lattice>
inline bool MetropolisMove(lattice* lat, double dE)
{
	if ( dE > 0. )
	{
		double rnd = lat->rngDistrib(lat->rngEngine);
		return (rnd < exp(-dE));
	}
	
	return true;
}

template<class lattice, class polymer>
inline void UpdateTAD(lattice* lat, polymer* pol,
					  unsigned long long* acceptCountPoly)
{
	double dE;
	
	pol->TrialMove(&dE);

	if ( pol->tad->legal )
	{
		double dEspe = pol->GetSpecificEnergy();
		bool acceptMove = MetropolisMove(lat, dE+dEspe);
	
		if ( acceptMove )
		{
			pol->AcceptMove();
			++(*acceptCountPoly);
		}
	}
}

template<>
inline void UpdateTAD<MCLiqLattice, MCHeteroPoly>(MCLiqLattice* lat, MCHeteroPoly* pol,
												  unsigned long long* acceptCountPoly)
{
	double dE;

	pol->TrialMove(&dE);

	if ( pol->tad->legal )
	{
		double dEcpl = pol->GetCouplingEnergy(lat->spinTable);
		bool acceptMove = MetropolisMove(lat, dE+dEcpl);
		
		if ( acceptMove )
		{
			pol->AcceptMove();
			++(*acceptCountPoly);
		}
	}
}

template<class lattice, class polymer>
inline void UpdateSpin(lattice*, polymer*, unsigned long long*) {}

template<>
inline void UpdateSpin<MCLiqLattice, MCPoly>(MCLiqLattice* lat, MCPoly*,
											 unsigned long long* acceptCountLiq)
{
	double dE;
			
	lat->TrialMove(&dE);
	bool acceptMove = MetropolisMove(lat, dE);

	if ( acceptMove )
	{
		lat->AcceptMove();
		++(*acceptCountLiq);
	}
}

template<>
inline void UpdateSpin<MCLiqLattice, MCHeteroPoly>(MCLiqLattice* lat, MCHeteroPoly* pol,
												   unsigned long long* acceptCountLiq)
{
	double dE;
		
	lat->TrialMove(&dE);
	
	double dEcpl = lat->GetCouplingEnergy(pol->tadHetTable);
	bool acceptMove = MetropolisMove(lat, dE+dEcpl);

	if ( acceptMove )
	{
		lat->AcceptMove();
		++(*acceptCountLiq);
	}
}


#endif /* MCUpdater_hpp */
