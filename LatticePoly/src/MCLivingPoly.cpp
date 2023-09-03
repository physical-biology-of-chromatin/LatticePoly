//
//  MCLivingPoly.cpp
//  LatticePoly
//
//  Created by mtortora/amithzafal on 27/04/2021.
//  Copyright © 2021 ENS Lyon. All rights reserved.
//

#include "MCLivingPoly.hpp"

#include <fstream>
#include <sstream>


MCLivingPoly::MCLivingPoly(MCLattice* _lat): MCHeteroPoly(_lat) {}

void MCLivingPoly::Init(int Ninit)
{
	MCHeteroPoly::Init(Ninit);

	for ( int vi = 0; vi < Ntot; ++vi )
	{
		painterTable[vi] = 0.;
		boostTable[vi] = 0.;
	}
    
	if ( !RestartFromFile )
	{
		std::ifstream painterFile(painterPath);
        
		std::string line;
		std::vector<std::pair<int, double>> painters;

		if ( !painterFile.good() )
			throw std::runtime_error("MCLivingPoly: Couldn't open file " + painterPath);
        
		while ( std::getline(painterFile, line) )
		{
			std::istringstream ss(line);

			int t;
			double charge;
            
			if ( ss >> t >> charge )
			{
				if ( t < Nchain )
                    tadConf[t].painter = charge;
				else
					throw std::runtime_error("MCLivingPoly: Found inconsistent values greater than chain length '" + line + "' in file " + painterPath);
			}
            
			else
				throw std::runtime_error("MCLivingPoly: Bad line '" + line + "' in file " + painterPath);
		}
        
		if ( propagationMode == 0 )
		{
			for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
			{
				if ( tad->type == 1 )
				{
					double rnd = lat->rngDistrib(lat->rngEngine);

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

	for ( auto tad = tadConf.begin(); tad != tadConf.end(); ++tad )
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi = (v ==0 ) ? tad->pos : lat->bitTable[v][tad->pos];
				
			painterTable[vi] += tad->painter;
			boostTable[vi] += tad->painter*tad->type;
			
		}
		//std::cout <<"het "<<hetTable[tad->pos] << std::endl;
	}
	//double tothet=0;
	//for ( int v = 0; v < Ntot; ++v )
	//	tothet+=painterTable[v];
	//std::cout <<"het "<<tothet << std::endl;
}

void MCLivingPoly::AcceptMove()
{
	MCHeteroPoly::AcceptMove();
	
	if ( tadTrial->painter != 0 )
	{
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			painterTable[vi1] -= tadTrial->painter;
			painterTable[vi2] += tadTrial->painter;
			
			boostTable[vi1] -= tadTrial->painter * tadTrial->type;
			boostTable[vi2] += tadTrial->painter * tadTrial->type;
		}		
	}
}



void MCLivingPoly::TrialMove(double* dE)
{   
    MCHeteroPoly::TrialMove(dE);

    if ( latticeType == "MCLattice" )  
		{  
       // PropagationMove();
		}

    if ( latticeType == "MCLiqLattice" )    
       {}
	   // LiqPropagationMove();               // to neglect the living part

	//if ( propagationMode == 1 )
	//	PropagationMove();
	
    //if ( propagationMode == 2 )
	//	LiqPropagationMove();
	
    //else
	//{
	//    if ( tadTrial->type == 2 )
	//    {
	//	    double rnd = lat->rngDistrib(lat->rngEngine);
		
	//	    if ( rnd < propRate / ((double) Ninter*Nmeas) )
	//	    {
	//		    tadTrial->type = 1;
			
	//		    for ( int v = 0; v < 13; ++v )
	//		    {
	//			    int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
	//			    ++hetTable[vi];
	//		    }
	//	    }
	//    }
	//}
}

void MCLivingPoly::PropagationMove()
{
	if ( tadTrial->type == 0 )
	{
		int numHet     = hetTable[tadTrial->pos];
		
		double painterCharge = painterTable[tadTrial->pos];
		double boostCharge   = boostTable[tadTrial->pos];

		double Rui = 0.;
		
		if ( painterAct != 0 )
		{
			int numCis = 0;
			
			double cisBoostCharge = 0.;
			double cisPainterCharge = 0.;

			if ( !tadTrial->isLeftEnd() )
			{
				MCTad* nb1 = tadTrial->neighbors[0];
				
				if ( nb1->type == 1 )
					numCis += 1;

				if ( nb1->painter != 0 )
				{
					cisPainterCharge += nb1->painter;
					
					if ( nb1->type == 1 )
						cisBoostCharge += nb1->painter;
				}
			}

			if ( !tadTrial->isRightEnd() )
			{
				MCTad* nb2 = tadTrial->neighbors[1];

				if ( nb2->type == 1 )
					numCis += 1;

				if ( nb2->painter != 0 )
				{
					cisPainterCharge += nb2->painter;

					if ( nb2->type == 1 )
						cisBoostCharge += nb2->painter;
				}
			}

			int numTrans = numHet-numCis;
			double onSite = tadTrial->painter;

			double transPainterCharge = painterCharge - cisPainterCharge;
			double transBoostCharge = boostCharge - cisBoostCharge;
			
			Rui = painterAct * (onSite + cisSpread * (cisPainterCharge + boost*cisBoostCharge) + 
            transSpread * (transPainterCharge + boost*transBoostCharge) + readerWriter * (cisSpread*numCis + transSpread*numTrans));
		}

		double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Rui )     // / ((double) Ninter*Nmeas)
		{
		    tadTrial->type = 1;
			
		    for ( int v = 0; v < 13; ++v )
		    {
			    int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
				boostTable[vi] += tadTrial->painter;
				
				++hetTable[vi];
			}
		}
	}
	
	else  //if ( tadTrial ->type == 1 )
	{
		double Riu = nucleoTurn;
		double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Riu )    // / ((double) Ninter*Nmeas)
		{
			tadTrial->type = 0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
				boostTable[vi] -= tadTrial->painter;

				--hetTable[vi];
			}
		}
	}
}

double MCLivingPoly::GetEffectiveEnergy() const
{
	double E1 = 0.;
	double E2 = 0.;

    double dEpainter    = 0.;  
    double dEcrosspaint = 0.;
    double dEcrosshet   = 0.; 
    double dEcrossInt   = 0.;  

	Jppp = sqrt(Jpppp)*sqrt(Jpp); 
	if ( Jns > 0. || Jpppp > 0. || Jppp > 0. )
	{
        
		for ( int v = 0; v < 13; ++v )
		{
			int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
			int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
			
			E1 += lat->bitTable[0][vi1];  //nonspecific
			E2 += lat->bitTable[0][vi2];
		}	

        if ( tadTrial->painter != 0. )
        {
            dEpainter  = painterTable[tadUpdater->vo] - painterTable[tadUpdater->vn];
            dEcrosshet = hetTable[tadUpdater->vo] - hetTable[tadUpdater->vn];
                
        }

        if ( tadTrial->type != 0. )
                
        {
            dEcrosspaint += painterTable[tadUpdater->vo] - painterTable[tadUpdater->vn];
        }
    }
        
    dEcrossInt = Jppp*(tadTrial->painter*dEcrosshet + tadTrial->type*dEcrosspaint);
    dEpainter  *= Jpppp * tadTrial->painter; 
		
	double dE =  + dEpainter + dEcrossInt + MCHeteroPoly::GetEffectiveEnergy(); //-Jns * (E2-E1)
	return dE;
}

double MCLivingPoly::GetCouplingEnergy(const int spinTable[Ntot]) const
{
    double dE1 = MCHeteroPoly::GetCouplingEnergy(spinTable);
	
	if ( ( Jlpp > 0. ) )
	{
        if ( tadTrial->painter != 0 )
		{
			double dE = 0.;
			for ( int v = 0; v < 13; ++v )
			{
				int vi1 = (v == 0) ? tadUpdater->vo : lat->bitTable[v][tadUpdater->vo];
				int vi2 = (v == 0) ? tadUpdater->vn : lat->bitTable[v][tadUpdater->vn];
				
                dE += spinTable[vi1];
				dE -= spinTable[vi2];
			}

			dE1 += Jlpp * dE * tadTrial->painter;

		}

	}

	if ( ( EV > 0. ) )
	{
		dE1 += EV * (spinTable[tadUpdater->vn]-spinTable[tadUpdater->vo]);
	}

	
	
	return dE1;
}


void MCLivingPoly::LiqPropagationMove()
{
	if ( tadTrial ->type == 0 )
    //lat->spinTable[tadTrial->pos] == 1 
	{
		int numLiqHet = 0;
		for ( int v = 0; v < 12; ++v ) //sum over neighbours
		{
			int vi =  lat->bitTable[v+1][tadTrial->pos];

            if ( lat->spinTable[vi] == 1 )
			    ++numLiqHet;
		}

		double painterCharge = painterTable[tadTrial->pos];
		double boostCharge   = boostTable[tadTrial->pos];

		double Rui = 0.;
		
		if ( painterAct != 0 )
		{
			int numCis = 0;
			
			double cisBoostCharge = 0.;
			double cisPainterCharge = 0.;

			if ( !tadTrial->isLeftEnd() )
			{
				MCTad* nb1 = tadTrial->neighbors[0];
				
				if ( lat->spinTable[nb1->pos] == 1 )
					numCis += 1;
	
				if ( nb1->type == 1 )
					cisBoostCharge += 1;
				
			}

			if ( !tadTrial->isRightEnd() )
			{
				MCTad* nb2 = tadTrial->neighbors[1];

				if ( lat->spinTable[nb2->pos] == 1 )
					numCis += 1;


				if ( nb2->type == 1 )
					cisBoostCharge += 1;
				
			}

			int numTrans  = numLiqHet-numCis;
			double onSite = lat->spinTable[tadTrial->pos]; 

			double transBoostCharge = boostCharge - cisBoostCharge;
			
			Rui = painterAct * (onSite + cisSpread * (numCis + boost*cisBoostCharge) + 
            transSpread * (numTrans + boost*transBoostCharge)); //+ readerWriter * (cisSpread*numCis + transSpread*numTrans));
		}

		double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Rui )     // / ((double) Ninter*Nmeas)
		{
		    tadTrial->type = 1;
			
		    for ( int v = 0; v < 13; ++v )
		    {
			    int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
				boostTable[vi] += tadTrial->painter;
				
				++hetTable[vi];
			}
		}
	}
	
	else  //if ( tadTrial ->type == 1 )
	{
		double Riu = nucleoTurn;
		double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Riu )    // / ((double) Ninter*Nmeas)
		{
			tadTrial->type = 0;
			
			for ( int v = 0; v < 13; ++v )
			{
				int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
				boostTable[vi] -= tadTrial->painter;

				--hetTable[vi];
			}
		}
	}
}
