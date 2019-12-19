//
//  MCLiqLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCLiqLattice_hpp
#define MCLiqLattice_hpp

#include <vector>

#include "MCLattice.hpp"


class MCLiqLattice: public MCLattice
{
public:
	void Init(std::mt19937_64&);
	void BleachSpins();
	void TrialMoveSpin(std::mt19937_64&, double*);
	void AcceptMoveSpin();
	void ToVTK(int);
	
	double GetSpinEnergy();
	double GetBindingEnergy(const int[Ntot]);
	double GetCouplingEnergy(const int[Ntot]);

	int spinTable[Ntot];
	
	int idxSpin1;
	int idxSpin2;

private:
	typedef struct {double dx,dy,dz;} disp;

	void DisplaceSpins();
	
	int nLiq;
	int spinIdTable[Ntot];
	int spinTypeTable[Ntot];
	
	std::vector<int> spinConf;
	std::vector<disp> spinDisp;
};


#endif /* MCLiqLattice_hpp */
