//
//  MCReplicPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 02/10/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCReplicPoly_hpp
#define MCReplicPoly_hpp

#include "MCHeteroPoly.hpp"


class MCReplicPoly: public MCHeteroPoly
{
public:
	MCReplicPoly(MCLattice*);
	
	void Init(int);
	void Update();

	void Replicate(int, int);
	void ReplicateTADs(std::vector<MCTad>::iterator, std::vector<MCTad>::iterator);
	void ReplicateBonds(std::vector<MCTad>::iterator, std::vector<MCTad>::iterator);
};


#endif /* MCReplicPoly_hpp */
