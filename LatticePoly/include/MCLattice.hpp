//
//  MCLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCLattice_hpp
#define MCLattice_hpp

#include <random>
#include <iostream>

#include "globals.hpp"


class MCLattice
{
public:
    MCLattice();
	
	void Init(int);
	void ToVTK(int);

	int opp[13];
	int nbNN[13][13][13];
	int bitTable[13][Ntot];
		
	double nbXYZ[3][13];
	double cTheta[13][13];
	double xyzTable[3][Ntot];
	
	std::mt19937_64 rngEngine;
	std::uniform_real_distribution<double> rngDistrib{0., 1.};
					
private:
	void InitConstArrays();
	void ReadInputArrays();
};


#endif /* MCLattice_hpp */
