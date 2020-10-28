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
	~MCPoly();

	void Init(int);
	void GenerateRandom(int);
	
	void ToVTK(int);
	void FromVTK(int);

	void TrialMove(double*);
	void AcceptMove();
	//void MoveFork(int,int);
	
	int Ntad;
	int Nbond;
	
	MCLattice* lat;
	MCTadUpdater* tadUpdater;
	
	std::vector<int> activeforks;

protected:
	MCTad* tadTrial;
	
	std::vector<MCTad> tadConf;
	std::vector<MCBond> tadTopo;
	
	void CreateBond(MCBond&);

private:
	double3 centreMass;
};


#endif /* MCPoly_hpp */
