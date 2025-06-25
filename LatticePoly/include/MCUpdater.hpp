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
		if(pol->Ntad==Nchain*2)
		{
			if(pol->individual_N_extruders!=0)
			{
				pol->unLoadExtruders();
				//int active_extruders_count=0;
				//for (int i=0 ; i < (int) pol->active_extruders.size() ; ++i)
				//active_extruders_count=active_extruders_count+ (int) pol->active_extruders.at(i)->N_loaded_extruders;
				int extruders_moves = pol->individual_N_extruders- (pol->individual_N_extruders-1)- pol->active_extruders.size();
				for (int i=0 ; i < extruders_moves ; ++i)
				{
					pol->LoadExtruders();
				}
			}
		}
	}
};

template<>
struct UpdateReplImpl<MCLattice, MCReplicPoly>
{
	static inline void _(MCLattice* , MCReplicPoly* pol)
	{
		//pol->OriginMove_implicit();
		pol->ForkMove();

		
		//extruders moves


		if(n_barriers>-1)
		{
			if(pol->individual_N_extruders!=0)
			{
				for (int i=0 ; i < (int) pol->active_extruders.size() ; ++i)
					if(pol->active_extruders.at(i)->status!=pol->active_extruders.at(i)->binding_site->status)
						throw std::runtime_error("extruder illegal");



				pol->unLoadExtruders();
				// check if any extruder has reached the end, if so, I desattatch and remove from active extruders
				for (int i=0 ; i < (int) pol->active_extruders.size() ; ++i)
				{

					auto LeftAnchor= pol->active_extruders.at(i);
					auto RightAnchor= pol->active_extruders.at(i)->binding_site;

 					if(LeftAnchor->neighbors[0]->isLeftEnd() or RightAnchor->neighbors[1]->isRightEnd())
					{

					//delete info of old extruder
					LeftAnchor->isCohesin=false;
					RightAnchor->isCohesin=false;
					RightAnchor->binding_site=nullptr;
					//LeftAnchor->binding_site=nullptr; gives segmentation error, in the context of 1D model (but also when is not a cohesin in 3D) it doesn't really matter 
					}
				}
				pol->active_extruders.erase(std::remove_if(pol->active_extruders.begin(), pol->active_extruders.end(), [](const MCTad* tad){return !tad->isCohesin;}), pol->active_extruders.end());



				int extruders_moves = pol->individual_N_extruders - pol->active_extruders.size();
				for (int i=0 ; i < extruders_moves ; ++i)
				{
					pol->LoadExtruders();

					if(Instantaneus_Extrusion==1)
					{
						if(pol->active_extruders.size()!=0)
						{
							bool left_stalled= (pol->active_extruders.back()->neighbors[0]->isCAR) or (pol->active_extruders.back()->neighbors[0]->isCohesin) or (pol->active_extruders.back()->neighbors[0]->isLeftEnd());
							bool right_stalled= (pol->active_extruders.back()->binding_site->neighbors[1]->isCAR) or (pol->active_extruders.back()->binding_site->neighbors[1]->isCohesin) or (pol->active_extruders.back()->binding_site->neighbors[1]->isRightEnd());

							while(left_stalled==0 or right_stalled==0)
							{
								pol->Move_Last_Extruders();
								left_stalled= (pol->active_extruders.back()->neighbors[0]->isCAR) or (pol->active_extruders.back()->neighbors[0]->isCohesin) or (pol->active_extruders.back()->neighbors[0]->isLeftEnd());
								right_stalled= (pol->active_extruders.back()->binding_site->neighbors[1]->isCAR) or (pol->active_extruders.back()->binding_site->neighbors[1]->isCohesin) or (pol->active_extruders.back()->binding_site->neighbors[1]->isRightEnd());

							}
						}
					}						
				}
				if(Instantaneus_Extrusion ==false)
					for (int i=0 ; i < (int) pol->active_extruders.size(); ++i)
						pol->Move_Extruders(); 
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
