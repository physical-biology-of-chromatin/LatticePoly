//
//  MCPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCPoly_hpp
#define MCPoly_hpp

#include "MCTad.hpp"


class MCPoly
{
public:
	MCPoly(MCLattice*);
	~MCPoly();

	void Init(std::mt19937_64&);
	void TrialMoveTAD(std::mt19937_64&, double*);
	void AcceptMoveTAD();
	void ToVTK(int);
	
	int tadType[Nchain];
	int config[2][Nchain];
	
	MCTad* tad;
	MCLattice* lat;
	
private:
	double centreMass[3];
};


#endif /* MCPoly_hpp */
