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
	
	void TrialMovePos(const int[Nchain], const int[Nchain-1]);
	void AcceptMovePos(int[Nchain], int[Nchain-1]) const;

	int n;
	int vo;
	int vn;
	
	double dE;
	
	bool legal;
	
private:
	void Reset();

	int nv1;
	int nv2;
	
	MCLattice* lat;
};


#endif /* MCTad_hpp */
