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


class MCTad
{
public:
	MCTad(MCLattice*);
	
	void RandomMove(const int[Nchain], const int[Nchain]);

	int n;
	int vo;
	int vn;
	int nv1;
	int nv2;
	
	double dE;
	
	bool legal;
	
private:
	void Reset();

	MCLattice* lat;
};


#endif /* MCTad_hpp */
