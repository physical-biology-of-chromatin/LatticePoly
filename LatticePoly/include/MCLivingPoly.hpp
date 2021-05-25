//
//  MCLivingPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 27/04/2021.
//  Copyright © 2021 ENS Lyon. All rights reserved.
//

#ifndef MCLivingPoly_hpp
#define MCLivingPoly_hpp

#include "MCHeteroPoly.hpp"


class MCLivingPoly: public MCHeteroPoly
{
public:
	MCLivingPoly(MCLattice*);
	
	void Init(int);	
	void TrialMove(double*);
    void PropagationMove(double*);
    void AcceptMove();

    double painterTable[Ntot];
};


#endif /* MCLivingPoly_hpp */
