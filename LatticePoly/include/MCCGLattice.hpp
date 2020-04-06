//
//  MCCGLattice.hpp
//  LatticePoly
//
//  Created by mtortora on 29/01/2020.
//  Copyright © 2020 ENS Lyon. All rights reserved.
//

#ifndef MCCGLattice_hpp
#define MCCGLattice_hpp

#include <vector>

#include "MCLattice.hpp"


class MCCGLattice: public MCLattice
{
public:
	virtual ~MCCGLattice() {};
	
	void Init(std::mt19937_64&);
	void ToVTK(int);

	virtual void TrialMove(std::mt19937_64&, double*);
	virtual void AcceptMove();
	
	virtual double GetCouplingEnergy(const int[Ntot]) const;

	bool legal;
	
	int idx1;
	int idx2;
	
	int spinTable[Ntot];
		
protected:
	double GetSpinEnergy() const;
	virtual double GetSpinDensity(int) const;
	
	int nLiq;

	int spinIdTable[Ntot];

	typedef struct {double dx,dy,dz;} disp;
	std::vector<disp> spinDisp;
};


#endif /* MCCGLattice_hpp */
