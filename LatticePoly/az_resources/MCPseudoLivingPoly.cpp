//
//  MCLivingPoly.cpp
//  LatticePoly
//
//  Created by mtortora on 27/04/2021.
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
	painterTable[vi] = 0;
    
    if ( !RestartFromFile )
	{
        std::ifstream painterFile(painterMonomer);
        
        std::string line;
        std::vector<std::pair<int, int>> painters;

        if ( !painterFile.good() )
            throw std::runtime_error("MCLivingPoly: Couldn't open file " + painterMonomer);
        
        while ( std::getline(painterFile, line) )
        {
            std::istringstream ss(line);
            
            int d1;
            int d2;
            
            if ( ss >> d1 >> d2 )
            {
                if ( (d1 >= 0) && (d2 >= 0) && (d1 < Nchain) && (d2 < Nchain) )
                    painters.push_back((d1 < d2) ? std::make_pair(d1, d2) : std::make_pair(d2, d1));
                else
                    throw std::runtime_error("MCLivingPoly: Found inconsistent domain boundaries '" + line + "' in file " + painterMonomer);
            }
            
            else
                throw std::runtime_error("MCLivingPoly: Bad line '" + line + "' in file " + painterMonomer);
        }
        
        for ( auto it = painters.begin(); it != painters.end(); ++it )
        {
            for ( int t = it->first; t <= it->second; ++t )
                tadConf[t].painter = 1.;
        }

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

            if ( tad->painter != 0 )
            {
                for ( int v = 0; v < 13; ++v )
                {
                    int vi = (v ==0 ) ? tad->pos : lat->bitTable[v][tad->pos];

                    painterTable[vi] += tad->painter;
                }
            }    
		}
	}
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
		}
	}
}

void MCLivingPoly::TrialMove(double* dE)
{
	MCHeteroPoly::TrialMove(dE);
	
	if (propagationMode == 1)
		PropagationMove(dE);
	
	else
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

void MCLivingPoly::PropagationMove(double* dE)
{
	if ( tadTrial->type == 0 )
	{
		int numTot      = 0;
		int numHet      = hetTable[tadTrial->pos];
        int numPaint    = painterTable[tadTrial->pos]; 
        int transBoost  = 0; 

        double Rui = 0; 
        if ( painterAct != 0 ) 
        {

            for ( int v = 0; v < 12; ++v )
            {
                int vi = lat->bitTable[v+1][tadTrial->pos];
                transBoost += tadConf[vi].type;
                numTot += lat->bitTable[0][vi];

            }
            
            int numCis = 0;
            //int numTrans = numTot-numHet;
            double cisPaint = 0; 
            int cisBoost  = 0;    
            
            if ( !tadTrial->isLeftEnd() )
            {
                    MCTad* nb1 = tadTrial->neighbors[0];
                    
                    if ( nb1->type == 1 )
                        numCis += 1;

                    if ( tadTrial->neighbors[0]->painter != 0 )
                        if ( tadTrial->neighbors[0]->type == 1 )
                            {   
                                cisPaint += cisSpread*( tadTrial->neighbors[0]->painter ) + boost*cisSpread*( tadTrial->neighbors[0]->painter );
                                cisBoost += 1;
                            }
                        else
                            cisPaint += cisSpread*( tadTrial->neighbors[0]->painter );
            }
                
            if ( !tadTrial->isRightEnd() )
            {
                    MCTad* nb2 = tadTrial->neighbors[1];
                    
                    if ( nb2->type == 1 )
                        numCis += 1;
                
                    if ( tadTrial->neighbors[1]->painter != 0 )
                        if ( tadTrial->neighbors[1]->type == 1 )
                            { 
                                cisPaint += cisSpread*( tadTrial->neighbors[1]->painter ) + boost*cisSpread*( tadTrial->neighbors[0]->painter );
                                cisBoost += 1;
                            }
                        else
                            cisPaint += cisSpread*( tadTrial->neighbors[1]->painter );
                             
            }
            
            int numTrans = numHet-numCis;
            double onsite = tadTrial->painter;

            transBoost = transBoost - cisBoost;
            double transPaint = transSpread*( numPaint - cisPaint - transBoost); 
            double transPaintBoost = transSpread*transBoost*boost;      
            
            Rui = painterAct*( onsite + cisPaint + transPaint + transPaintBoost + readerWriter*( cisSpread*numCis + transSpread*numTrans )); 
        }
        
        double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Rui / ((double) Ninter*Nmeas) )
		{
		    tadTrial->type = 1;
			
		    for ( int v = 0; v < 13; ++v )
		    {
			    int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
                ++hetTable[vi];

			}
        }    

	}
    else  //if ( tadTrial ->type ==1 )
    {
        double Riu = nucleoturn;
        double rnd = lat->rngDistrib(lat->rngEngine);
		
		if ( rnd < Riu / ((double) Ninter*Nmeas) )
		{
		    tadTrial->type = 0;
			
		    for ( int v = 0; v < 13; ++v )
		    {
			    int vi = (v == 0) ? tadTrial->pos : lat->bitTable[v][tadTrial->pos];
				
                ++hetTable[vi];

			}
        }                   
    }
    
}
