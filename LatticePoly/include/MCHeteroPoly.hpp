//
//  MCHeteroPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCHeteroPoly_hpp
#define MCHeteroPoly_hpp

#include "MCPoly.hpp"


class MCHeteroPoly: public MCPoly
{
public:
	MCHeteroPoly(MCLattice*);
	
	void Init(int);
	void AcceptMove();
	
	double GetEffectiveEnergy() const;
	double GetCouplingEnergy(const int[Ntot]) const;
	//void MoveFork(int,int);
	std::vector<int> activeforks;

	
	int hetTable[Ntot];
};


#endif /* MCHeteroPoly_hpp */
