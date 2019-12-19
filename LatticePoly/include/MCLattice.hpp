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
	
	void Init(std::mt19937_64&);
	void ToVTK(int);

	int opp[13];
	int nbNN[13][13][13];
	int bitTable[13][Ntot];
		
	double nbXYZ[3][13];
	double cTheta[13][13];
	double xyzTable[3][Ntot];
		
	std::string dPath;
			
private:
	void InitConstArrays();
	void ReadInputArrays();
};


#endif /* MCLattice_hpp */
