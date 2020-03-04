//
//  MCLiqLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCLiqLattice_hpp
#define MCLiqLattice_hpp

#include "MCCGLattice.hpp"


class MCLiqLattice: public MCCGLattice
{
public:
	void Init(std::mt19937_64&);
	
	void TrialMove(std::mt19937_64&, double*);
	void AcceptMove();
	
	double GetCouplingEnergy(const int[Ntot]) const;

protected:
	double GetSpinEnergy() const;
	double GetSpinDensity(int) const;
	
private:
	void DisplaceSpins();
	
	std::vector<int> spinConf;	
};


#endif /* MCLiqLattice_hpp */
