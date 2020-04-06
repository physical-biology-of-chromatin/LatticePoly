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


inline bool MetropolisMove(std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
						   double dE)
{
	if ( dE > 0. )
	{
		double rnd = rngDistrib(rngEngine);
		return (rnd < exp(-dE));
	}
	
	return true;
}

inline bool MCCGMove(std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
					 double dE, int spin1, int spin2)
{
	double rnd = rngDistrib(rngEngine);
	
	double c = spin1 * (Qcg - spin2) / (double) Qcg;
	double rate = (dE > 0.) ? c * exp(-dE) : c;
	
	return (rnd < rate);
}

// Update functions with template specializations
template<class lattice, class polymer>
inline void UpdateTAD(lattice*, polymer* pol,
					  std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
					  unsigned long long* acceptCountPoly)
{
	bool acceptMove;
	double dE, dEspe;
	
	pol->TrialMove(rngEngine, &dE);

	if ( pol->tad->legal )
	{
		dEspe = pol->GetSpecificEnergy();
		acceptMove = MetropolisMove(rngEngine, rngDistrib, dE+dEspe);
	
		if ( acceptMove )
		{
			pol->AcceptMove();
			(*acceptCountPoly)++;
		}
	}
}

template<class lattice, class polymer>
inline void UpdateSpin(lattice*, polymer*,
					   std::mt19937_64&, std::uniform_real_distribution<double>&,
					   unsigned long long*) {}

template<>
inline void UpdateTAD<MCCGLattice, MCHeteroPoly>(MCCGLattice* lat, MCHeteroPoly* pol,
												 std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
												 unsigned long long* acceptCountPoly)
{
	double dE;
	bool acceptMove;
		
	pol->TrialMove(rngEngine, &dE);

	if ( pol->tad->legal )
	{
		double dEcpl = pol->GetCouplingEnergy(lat->spinTable);
		acceptMove = MetropolisMove(rngEngine, rngDistrib, dE+dEcpl);
		
		if ( acceptMove )
		{
			pol->AcceptMove();
			(*acceptCountPoly)++;
		}
	}
}

template<>
inline void UpdateSpin<MCCGLattice, MCPoly>(MCCGLattice* lat, MCPoly*,
											std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
											unsigned long long* acceptCountLiq)
{
	bool acceptMove;
	double dE;
			
	lat->TrialMove(rngEngine, &dE);

	if ( lat->legal )
	{
		int spin1 = lat->spinTable[lat->idx1];
		int spin2 = lat->spinTable[lat->idx2];

		acceptMove = MCCGMove(rngEngine, rngDistrib, dE, spin1, spin2);

		if ( acceptMove )
		{
			lat->AcceptMove();
			(*acceptCountLiq)++;
		}
	}
}

template<>
inline void UpdateSpin<MCCGLattice, MCHeteroPoly>(MCCGLattice* lat, MCHeteroPoly* pol,
												  std::mt19937_64& rngEngine, std::uniform_real_distribution<double>& rngDistrib,
												  unsigned long long* acceptCountLiq)
{
	bool acceptMove;
	double dE;
		
	lat->TrialMove(rngEngine, &dE);

	if ( lat->legal )
	{
		int spin1 = lat->spinTable[lat->idx1];
		int spin2 = lat->spinTable[lat->idx2];
		
		double dEcpl = lat->GetCouplingEnergy(pol->tadHetTable);
		
		acceptMove = MCCGMove(rngEngine, rngDistrib, dE+dEcpl, spin1, spin2);

		if ( acceptMove )
		{
			lat->AcceptMove();
			(*acceptCountLiq)++;
		}
	}
}


#endif /* MCUpdater_hpp */
