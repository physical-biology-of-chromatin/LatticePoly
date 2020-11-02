//
//  MCUpdater.hpp
//  LatticePoly
//
//  Created by mtortora on 28/02/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCUpdater_hpp
#define MCUpdater_hpp

#include "MCReplicPoly.hpp"


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


/* struct containers are required to circumvent function template partial specialization */

// UpdateTAD template specialisations
template<class lattice, class polymer>
struct UpdateTADImpl
{
	static inline void _(lattice* lat, polymer* pol, unsigned long long* acceptCount)
	{
		double dE;

		pol->TrialMove(&dE);


		if ( pol->tadUpdater->legal )
		{
			double dEcpl = pol->GetCouplingEnergy(lat->spinTable);
			bool acceptMove = MetropolisMove(lat, dE+dEcpl);
			
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
	}
};

template<class polymer>
struct UpdateTADImpl<MCLattice, polymer>
{
	static inline void _(MCLattice* lat, polymer* pol, unsigned long long* acceptCount)
	{
		double dE;
		
		pol->TrialMove(&dE);

		if ( pol->tadUpdater->legal )
		{
			double dEeff = pol->GetEffectiveEnergy();
			bool acceptMove = MetropolisMove(lat, dE+dEeff);
		
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
	}
};

template<>
struct UpdateTADImpl<MCLattice, MCPoly>
{
	static inline void _(MCLattice* lat, MCPoly* pol, unsigned long long* acceptCount)
	{
		double dE;
		
		pol->TrialMove(&dE);

		if ( pol->tadUpdater->legal )
		{
			bool acceptMove = MetropolisMove(lat, dE);
		
			if ( acceptMove )
			{
				pol->AcceptMove();
				++(*acceptCount);
			}
		}
	}
};


// UpdateSpin template specialisations
template<class lattice, class polymer>
struct UpdateSpinImpl
{
	static inline void _(lattice* lat, polymer* pol, unsigned long long* acceptCount)
	{
		double dE;
			
		lat->TrialMove(&dE);
		
		double dEcpl = lat->GetCouplingEnergy(pol->hetTable);
		bool acceptMove = MetropolisMove(lat, dE+dEcpl);

		if ( acceptMove )
		{
			lat->AcceptMove();
			++(*acceptCount);
		}
	}
};

template<class polymer>
struct UpdateSpinImpl<MCLattice, polymer>
{
	static inline void _(MCLattice*, polymer*, unsigned long long*) {}
};

// UpdateFork template specialisations

template<class lattice,class polymer>
struct UpdateForkImpl
{
	static inline void _(lattice* lat, polymer* pol) {}
};

template<class lattice>
struct UpdateForkImpl<lattice, MCReplicPoly>
{
	static inline void _(lattice* lat, MCReplicPoly* pol)
	{
		
		double rnd2 = lat->rngDistrib(lat->rngEngine);
		if(rnd2<Replicationrate){
			int t = lat->rngEngine() % (int) pol->activeforks.size();
			pol->MoveFork(pol->activeforks[t],t);
		}
	}
};

// createFork template specialisations

template<class lattice,class polymer>
struct CreateForkImpl
{
	static inline void _(lattice* lat, polymer* pol) {}
};

template<class lattice>
struct CreateForkImpl<lattice, MCReplicPoly>
{
	static inline void _(lattice* lat, MCReplicPoly* pol)
	{
		double rnd2 = lat->rngDistrib(lat->rngEngine);
		if(rnd2<Originrate){
			pol->CreateFork();

		}
	}
};

// Wrapper functions
template<class lattice, class polymer>
inline void UpdateTAD(lattice* lat, polymer* pol, unsigned long long* acceptCount)
{
	UpdateTADImpl<lattice, polymer>::_(lat, pol, acceptCount);
}

template<class lattice, class polymer>
inline void UpdateSpin(lattice* lat, polymer* pol, unsigned long long* acceptCount)
{
	UpdateSpinImpl<lattice, polymer>::_(lat, pol, acceptCount);
}

template<class lattice,class polymer>
inline void UpdateFork(lattice* lat, polymer* pol)
{
	UpdateForkImpl<lattice, polymer>::_(lat, pol);
}

template<class lattice,class polymer>
inline void CreateFork(lattice* lat, polymer* pol)
{
	CreateForkImpl<lattice, polymer>::_(lat, pol);
}
#endif /* MCUpdater_hpp */
