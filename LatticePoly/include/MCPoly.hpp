//
//  MCPoly.hpp
//  LatticePoly
//
//  Created by mtortora on 30/11/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef MCPoly_hpp
#define MCPoly_hpp

#include "MCTadUpdater.hpp"


class MCPoly
{
public:
	MCPoly(MCLattice*);
	virtual ~MCPoly();

	void Init(int);
	void GenerateHedgehog(int);
	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	
	int Ntad;
	int Nbond;
	
	MCLattice* lat;
	MCTadUpdater* tadUpdater;
		
protected:
	MCTad* tadTrial;
	
	std::vector<MCTad> tadConf;
	std::vector<MCBond> tadTopo;
	
	void SetBond(MCBond&);
		
	virtual vtkSmartPointer<vtkPolyData> GetVTKData();
	virtual void SetVTKData(const vtkSmartPointer<vtkPolyData>);
	
private:
	std::array<double3, 2> centerMass;

	void SetPBCCenterMass(std::vector<double3>&);
	void BuildPBCPair(std::vector<MCTad*>&, std::vector<double3>&, MCTad*, MCTad*);

	std::vector<double3> BuildUnfoldedConf();
};


#endif /* MCPoly_hpp */
