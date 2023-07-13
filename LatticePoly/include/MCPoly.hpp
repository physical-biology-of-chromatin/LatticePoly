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
#include <vtkMultiBlockDataSet.h>


class MCPoly
{
public:
	MCPoly(MCLattice*);
	virtual ~MCPoly();

	void Init(int);
	void GenerateHedgehog(int);
	void GenerateRing(int);
	void Update_rcms_before_separation();

	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	bool PrintCohesins();

	
	int Ntad;
	int Nbond;

	//vtkSmartPointer<vtkMultiBlockDataSet> mbds;


	std::vector<double3> BuildUnfoldedConf();

	


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

	void FixPBCPair(std::vector<double3>&, MCTad*, MCTad*);
	void FixPBCCenterMass(std::vector<double3>&);

};


#endif /* MCPoly_hpp */
