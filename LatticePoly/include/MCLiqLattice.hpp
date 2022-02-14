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
#include "MCTad.hpp"


class MCLiqLattice: public MCLattice
{
public:
	void Init(int);
	
	void ToVTK(int);
	void FromVTK(int);
	
	void TrialMove(double*);
	void AcceptMove();
	
	double GetCouplingEnergy(const int[Ntot]) const;
	double GetCouplingForkEnergy(const std::vector<int>) const;

	
	int spinTable[Ntot];
	int OriginCheck(std::vector<int>);
	void unLockSpins(std::vector<int>);
	void LockSpin();
	
	int nLiq;
	int n;
	bool stop_update = false;
	std::vector<int> SpinLocked;
	std::vector<std::vector<int>> fork_pos;
	std::vector<int> spinConf;

	
private:
	void GenerateRandom();
	void GenerateDroplets();
	
	void DisplaceSpins();
	
	double GetSpinEnergy() const;
	
	int lookupTable[Ntot];
	
	int v1;
	int v2;

	std::vector<double3> spinDisp;
};


#endif /* MCLiqLattice_hpp */

