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
	
	void Init(std::mt19937_64&);
	void TrialMoveSpinTAD(std::mt19937_64&, double*);
	void AcceptMoveSpinTAD();
	void AcceptMoveTAD();
	
	double GetSpecificEnergy() const;
	double GetCouplingEnergy(const int[Ntot]) const;
	
	int tadHetTable[Ntot];
};


#endif /* MCHeteroPoly_hpp */
