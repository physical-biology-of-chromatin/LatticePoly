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
	
	set = false;
}

MCTad::MCTad(): bonds(), neighbors()
{
	pos = -1;
	
	type = 0;
	links = 0;
	status = 0;
}

MCTad& MCTad::operator= (MCTad& tad)
{
	if ( &tad != this )
	{
		pos = tad.pos;
		type = tad.type;
		
		status = +1;
		tad.status = -1;
	}
	
	return *this;
}
