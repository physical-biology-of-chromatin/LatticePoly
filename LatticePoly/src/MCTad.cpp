//
//  MCTad.cpp
//  LatticePoly
//
//  Created by mtortora on 02/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include "MCTad.hpp"


MCBond::MCBond()
{
	id1 = -1;
	id2 = -1;
	dir = -1;
	
	isSet = false;
}

MCTad::MCTad(): bonds(), neighbors()
{
	pos = -1;
	type = 0;
	domain = -1;
	insulator_type = {};
	links = 0;
	status = 0;
	SisterID = -1;
	isCohesin=false;
	isCAR=false;
	isrDNA=false;

	CAR_weight=0;
	N_loaded_extruders=0;
	insulator_score=0.0;
	

	isCentromere=false;
	binding_particle=-1;
	binding_site=nullptr;
}

MCTad& MCTad::operator= (MCTad& tad)
{
	if ( &tad != this )
	{
		pos = tad.pos;
		type = tad.type;
		domain=tad.domain;
		
		SisterID = tad.SisterID;
		
		status = +1;
		tad.status = -1;

	}
	
	return *this;
}
