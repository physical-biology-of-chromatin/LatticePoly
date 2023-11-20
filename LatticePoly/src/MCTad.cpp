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
	sisterID = -1;
	
	type = 0;
	links = 0;
	status = 0;

	density = 0;
	homdensity = 0;

	isCohesin = 0;
	isBarrier = 0;
	loops     = 0;
	loopDir   = -1;
}

MCTad& MCTad::operator= (MCTad& tad)
{
	if ( &tad != this )
	{
		pos = tad.pos;
		type = tad.type;
		
		sisterID = tad.sisterID;
		
		status = +1;
		tad.status = -1;

		density = tad.density;
		homdensity = tad.homdensity;

		isCohesin = tad.isCohesin;
		isBarrier = tad.isBarrier;
		loops     = tad.loops;
		loopDir   = tad.loopDir;
	}
	
	return *this;
}
