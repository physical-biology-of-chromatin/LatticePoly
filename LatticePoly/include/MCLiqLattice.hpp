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
	void Init();
	void ToVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	
	double GetCouplingEnergy(const int[Ntot]) const;

	int spinTable[Ntot];

private:
	void DisplaceSpins();
	double GetSpinEnergy() const;

	int nLiq;
	int spinIdTable[Ntot];

	int idx1;
	int idx2;
	
	typedef struct {double dx,dy,dz;} disp;
	
	std::vector<disp> spinDisp;
	std::vector<int> spinConf;	
};


#endif /* MCLiqLattice_hpp */
