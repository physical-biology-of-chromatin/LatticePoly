//
//  MCLiqLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCLiqLattice_hpp
#define MCLiqLattice_hpp

#include "MCLattice.hpp"


class MCLiqLattice: public MCLattice
{
public:
	void Init(int);
	
	void GenerateRandom();
	void GenerateDroplets();

	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	
	double GetCouplingEnergy(const int[Ntot]) const;

	int spinTable[Ntot];

private:
	void DisplaceSpins();
	double GetSpinEnergy() const;

	int spinIdTable[Ntot];

	int v1;
	int v2;
	
	typedef struct {double dx,dy,dz;} disp;
	
	std::vector<int> spinConf;
	std::vector<disp> spinDisp;
};


#endif /* MCLiqLattice_hpp */
