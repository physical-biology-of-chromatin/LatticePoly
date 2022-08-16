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
	void GenerateRabl(int);

	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	void OriginMove(MCTad*);
	void ForkMove();
	void GenerateCAR();
	
	int Ntad;
	int Nbond;
	std::vector<MCTad*> activeForks;
	int NbindedForks;
	int NbindedCohesin;


	std::vector<MCTad*> activeOrigins;
	std::vector<int> MergedForkPos;

	


	


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

	std::vector<double3> BuildUnfoldedConf();
};


#endif /* MCPoly_hpp */
