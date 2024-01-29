//
//  MCHeteroPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCHeteroPoly_looped_hpp
#define MCHeteroPoly_looped_hpp

#include "MCHeteroPoly.hpp"


class MCHeteroPoly_looped: public MCHeteroPoly
{
public:
	MCHeteroPoly_looped(MCLattice*);
	
	void Init(int);	
	void AcceptMove();
	void BuildLoopTable();

	double GetEffectiveEnergy() const;
	double GetCouplingEnergy(const int[Ntot]) const;

	double hetTable_insulator[9][Ntot];

protected:
	virtual vtkSmartPointer<vtkPolyData> GetVTKData();
	virtual void SetVTKData(const vtkSmartPointer<vtkPolyData>);
};


#endif /* MCHeteroPoly_hpp */
