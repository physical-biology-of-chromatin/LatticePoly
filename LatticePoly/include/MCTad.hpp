//
//  MCTad.hpp
//  LatticePoly
//
//  Created by mtortora on 02/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCTad_hpp
#define MCTad_hpp

#include "MCLiqLattice.hpp"


struct MCBond
{
	MCBond();

	int id1;
	int id2;
	int dir;
	
	bool isSet;
};


struct MCTad
{
	MCTad();
	MCTad& operator= (MCTad&);
				 
	inline bool isLeftEnd() const {return (links == 1) ? !bonds[0] : false;};
	inline bool isRightEnd() const {return (links == 1) ? !bonds[1] : false;};
	
	inline bool isFork() const {return links == 3;};
	inline bool isLeftFork() const {return isFork() ? (this == neighbors[2]->neighbors[0]) : false;};
	inline bool isRightFork() const {return isFork() ? (this == neighbors[2]->neighbors[1]) : false;};
	
	int pos;
	double CAR_weight;
	int type;
	int domain;
	int insulator_type;
	int links;
	int status;
	int SisterID;
	bool isCohesin;
	bool isCAR;

	int binding_particle;
	bool isCentromere;
	int N_loaded_extruders;
	MCTad* binding_site;
	
	
	MCBond* bonds[3];
	MCTad* neighbors[3];
};


#endif /* MCTad_hpp */
