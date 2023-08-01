//
//  MCHeteroPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 03/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCHeteroPoly_hpp
#define MCHeteroPoly_hpp

#include "MCPoly.hpp"


class MCHeteroPoly: public MCPoly
{
public:
	MCHeteroPoly(MCLattice*);
	
	void Init(int);	
	void AcceptMove();
	void BuildHetTable();

	double GetEffectiveEnergy() const;
	double GetCouplingEnergy(const int[Ntot]) const;
	int hetTable[Ntot];
	int hetTable_tads[51][Ntot];
	int hetTable_domain[2][Ntot];
	int hetTable_insulator[3][Ntot];


protected:
	virtual vtkSmartPointer<vtkPolyData> GetVTKData();
	virtual void SetVTKData(const vtkSmartPointer<vtkPolyData>);
};


#endif /* MCHeteroPoly_hpp */
