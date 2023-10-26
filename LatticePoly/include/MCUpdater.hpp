//
//  MCUpdater.hpp
//  LatticePoly
//
//  Created by mtortora on 28/02/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCUpdater_hpp
#define MCUpdater_hpp

#include "MCLivingPoly.hpp"
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

template<>
struct UpdateTADImpl<MCLiqLattice, MCReplicPoly>
{
	static inline void _(MCLiqLattice* lat, MCReplicPoly* pol, unsigned long long* acceptCount)
	{
		double dE;
		
		pol->TrialMove(&dE);
		
		//double dEcpl = pol->GetCouplingEnergy(lat->spinConf);
		double dEeff = pol->GetEffectiveEnergy();

		if ( pol->tadUpdater->legal )
		{
			bool acceptMove = MetropolisMove(lat, dE+dEeff);
			
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
				if( polyType == "MCLivingPoly")
					static_cast<MCLivingPoly*>(pol)->AcceptMove();
				else if ( polyType == "MCHeteroPoly")
					static_cast<MCHeteroPoly*>(pol)->AcceptMove();
				else if ( polyType == "MCHeteroPoly_looped")
					static_cast<MCHeteroPoly_looped*>(pol)->AcceptMove();
				else if ( polyType == "MCReplicPoly")
					static_cast<MCReplicPoly*>(pol)->AcceptMove();
				else
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
// UpdateSpin template specialisations
template<>
struct UpdateSpinImpl<MCLiqLattice, MCReplicPoly>
{
	static inline void _(MCLiqLattice* lat, MCReplicPoly* pol, unsigned long long* acceptCount)
	{
		

		if(pol->Spin_pos_toDelete.size()>0)
			for (int i=0 ; i < (int) pol->Spin_pos_toDelete.size() ; ++i)
				lat->DeleteSpin(pol->Spin_pos_toDelete.at(i));

		pol->Spin_pos_toDelete.clear();



		double dE;
		
		if(lat->nLiq>0)
		{
			lat->TrialMove(&dE);
			
			if(lat->stop_update==false)
			{
				//double dEcpl = lat->GetCouplingEnergy(pol->hetTable);
				double dEcpl =0;
				bool acceptMove = MetropolisMove(lat, dE+dEcpl);
				
				if ( acceptMove )
				{
					lat->AcceptMove();
					
					++(*acceptCount);
				}
			}
		}
		if(pol->Spin_pos_toCreate.size()>0)
		{
			std::cout << "MERGING  "  << std::endl;
			
			for (int i=0 ; i < (int) pol->Spin_pos_toCreate.size() ; ++i)
				lat->CreateSpin(pol->Spin_pos_toCreate.at(i));
		}
		pol->Spin_pos_toCreate.clear();
	}
};

template<class polymer>
struct UpdateSpinImpl<MCLattice, polymer>
{
	static inline void _(MCLattice*, polymer*, unsigned long long*) {}
};

template<class lattice, class polymer>
struct UpdateReplImpl
{
	static inline void _(lattice* , polymer* ) {}

};

template<>
struct UpdateReplImpl<MCLiqLattice, MCReplicPoly>
{
	static inline void _(MCLiqLattice* lat, MCReplicPoly* pol)
	{
		pol->OriginMove_explicit(lat->spinTable);
		pol->ForkMove();
		
		//extruders moves
		if(pol->Ntad==Nchain*2 or 0==0)
		{
			if(N_extruders!=0)
			{
				pol->unLoadExtruders();
				int active_extruders_count=0;
				for (int i=0 ; i < (int) pol->active_extruders.size() ; ++i)
					active_extruders_count=active_extruders_count+ (int) pol->active_extruders.at(i)->N_loaded_extruders;
				int extruders_moves = N_extruders - active_extruders_count;
				for (int i=0 ; i < extruders_moves ; ++i)
					pol->LoadExtruders();
			}
		}
	}
};

template<>
struct UpdateReplImpl<MCLattice, MCReplicPoly>
{
	static inline void _(MCLattice* , MCReplicPoly* pol)
	{
		pol->OriginMove_implicit();
		pol->ForkMove();
		//extruders moves
		if(pol->Ntad==Nchain*2)
		{
			if(N_extruders!=0)
			{
				pol->unLoadExtruders();
				int active_extruders_count=0;
				for (int i=0 ; i < (int) pol->active_extruders.size() ; ++i)
					active_extruders_count=active_extruders_count+ (int) pol->active_extruders.at(i)->N_loaded_extruders;
				int extruders_moves = N_extruders - active_extruders_count;
				for (int i=0 ; i < extruders_moves ; ++i)
				{
					std::cout <<  " LOADING n " << i << std::endl;
					pol->LoadExtruders();

				}
			}
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

template<class lattice, class polymer>
inline void UpdateRepl(lattice* lat, polymer* pol)
{
	UpdateReplImpl<lattice, polymer>::_(lat, pol);
}
#endif /* MCUpdater_hpp */
