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
	void OriginMove(MCTad*);
	void ForkMove();
	void GenerateCAR();
	
	int Ntad;
	int Nbond;
	std::vector<MCTad*> activeForks;
	std::vector<MCTad*> activeOrigins;
	
	std::vector<MCBond> IntraBonds;
	std::vector<MCTad*> CAR;
	std::vector<MCTad*> interCAR;

	


	MCLattice* lat;
	MCTadUpdater* tadUpdater;

		
protected:
	MCTad* tadTrial;
	
	std::vector<MCTad> tadConf;
	std::vector<MCBond> tadTopo;
<<<<<<< HEAD

	void CreateBond(MCBond&);
	void FixPBCPair(std::vector<double3>&, int, int);
	
	virtual std::vector<double3> GetPBCConf();
	double3 GetPBCCenterMass(std::vector<double3>::iterator, std::vector<double3>::iterator);
=======
	
	void SetBond(MCBond&);
		
	virtual vtkSmartPointer<vtkPolyData> GetVTKData();
	virtual void SetVTKData(const vtkSmartPointer<vtkPolyData>);
>>>>>>> origin/master
	
private:
	std::array<double3, 2> centerMass;

	void FixPBCPair(std::vector<double3>&, MCTad*, MCTad*);
	void FixPBCCenterMass(std::vector<double3>&);

	std::vector<double3> BuildUnfoldedConf();
};


#endif /* MCPoly_hpp */
